# Helper: build per-frame displacement Vector{CartesianIndex} from a 2×T raw matrix
function _dx_vec(raw::AbstractMatrix{<:Integer}, lengthT)
    [CartesianIndex(raw[1,t], raw[2,t]) for t in 1:lengthT]
end

# Helper: run compositeimages and return plain Float64 array (nrow × ncol × t)
function _composite_float(imgsz, template, centers, intensity; dx=nothing)
    result = compositeimages(imgsz, template, centers, intensity; dx)
    Float64.(Array(result.images))
end

function gaussian1D(ncells, lengthT, nevents, bias=0.1, SNR=10) # should lengthT > 50
    S₁ = gaussiantemplate(Float64, (2.0, 2.0))   # 2-D template; extract middle column
    mid = (size(S₁, 2) + 1) ÷ 2
    S₁ = S₁[:, mid]                              # 1-D gaussian profile
    T₀ = sparseevents((lengthT, ncells), round(Int, nevents))
    T₀[37, 1:min(2, ncells)] = [0.2, 0.7][1:min(2, ncells)]
    if ncells > 3
        T₀[50, 3:4] = [0.2, 0.7]
    elseif ncells > 2
        T₀[20, 2:3] = [0.2, 0.7]
    end
    S₀ = zeros(eltype(S₁), 2*ncells*mid, ncells)
    for i = 1:ncells
        S₀[(i-1)*length(S₁)+1:(i-1)*length(S₁)+length(S₁), i] = S₁
    end
    img = S₀ * T₀'
    signalpwr = sum(img.^2) / length(img)
    bg = randn(size(img)) .* sqrt(signalpwr / 10^(SNR/10)) .+ maximum(img)*bias
    img = img + bg
    img, S₀, T₀', bg
end

# bias is a percentage of maximum(img₂) before noise
function gaussian2D(sigma, imgsz::NTuple{2}, lengthT, revent=10; fovsz::NTuple{2}=imgsz, jitter=0,
        drift=0, bias=0.1, SNR=10, useCalciumT=false, orthogonal=true, overlaplevel=1,
        inhibitindices=0, gtincludebg=false) # should lengthT > 50
    template = gaussiantemplate(Float64, (sigma, sigma))
    centers  = spreadcells(fovsz, sigma)          # Vector{CartesianIndex{2}}
    imagecenter = fovsz .÷ 2
    overlap_r = imagecenter[1] - round(Int, min(overlaplevel, imagecenter[1]) / 1.5)
    overlap_c = imagecenter[2] + min(overlaplevel, imagecenter[2])
    !orthogonal && push!(centers, CartesianIndex(overlap_r, overlap_c))
    ncells  = length(centers)
    nevents = round(Int, revent * ncells * lengthT / 100)

    if useCalciumT
        fr = firingrate((lengthT, ncells); lambda=0.01)
        T₂ = calciumtransient(fr, [0.85]; noisesigma=0.0, alpha=1, baseline=0.0)
    else
        T₂ = sparseevents((lengthT, ncells), nevents)
    end

    dx_raw = rand(-jitter:jitter, 2, lengthT)
    dx_raw = map(x -> round(Int, x), dx_raw .+ range(0, stop=drift, length=lengthT)')
    dx_ci  = _dx_vec(dx_raw, lengthT)

    if !isempty(inhibitindices)
        for inhibitidx in inhibitindices
            if inhibitidx > 0 && inhibitidx <= ncells
                T₂[:, inhibitidx] .*= -0.5    # scale and invert
            end
        end
        img_wobias = _composite_float(imgsz, template, centers, T₂; dx=dx_ci)
        signalpwr  = sum(img_wobias.^2)
        for inhibitidx in inhibitindices
            if inhibitidx > 0 && inhibitidx <= ncells
                T₂[:, inhibitidx] .+= -minimum(T₂[:, inhibitidx])   # shift to non-negative
            end
        end
    end

    img₂      = _composite_float(imgsz, template, centers, T₂; dx=dx_ci)
    signalpwr  = isempty(inhibitindices) ? sum(img₂.^2) : signalpwr
    bg         = maximum(img₂) * bias
    noise      = randn(size(img₂)); noisepwr0 = sum(noise.^2)
    noise    .*= sqrt(signalpwr / 10^(SNR/10) / noisepwr0)
    img₂       = img₂ .+ bg + noise
    img₂a      = AxisArray(img₂, :x, :y, :time)
    gtW, gtH, gtWimgc = makegt(template, centers, T₂, imgsz, ncells, bg; gtincludebg)
    gtbg  = copy(bg)
    imgrs = Matrix(reshape(img₂a, prod(imgsz), nimages(img₂a)))
    ncells, imgrs, img₂, gtW, gtH, gtWimgc, gtbg
end

function gaussian2D_AmplitudeSNR(sigma, imgsz::NTuple{2}, lengthT, revent=10; fovsz::NTuple{2}=imgsz,
        jitter=0, drift=0, bias=0.1, SNR=10, useCalciumT=false, orthogonal=true, overlaplevel=1,
        inhibitindices=0, gtincludebg=false) # should lengthT > 50
    template = gaussiantemplate(Float64, (sigma, sigma))
    centers  = spreadcells(fovsz, sigma)
    imagecenter = fovsz .÷ 2
    overlap_r = imagecenter[1] - round(Int, min(overlaplevel, imagecenter[1]) / 1.5)
    overlap_c = imagecenter[2] + min(overlaplevel, imagecenter[2])
    !orthogonal && push!(centers, CartesianIndex(overlap_r, overlap_c))
    ncells  = length(centers)
    nevents = round(Int, revent * ncells * lengthT / 100)

    if useCalciumT
        fr = firingrate((lengthT, ncells); lambda=0.01)
        T₂ = calciumtransient(fr, [0.85]; noisesigma=0.0, alpha=1, baseline=0.0)
    else
        T₂ = sparseevents((lengthT, ncells), nevents)
    end

    dx_raw = rand(-jitter:jitter, 2, lengthT)
    dx_raw = map(x -> round(Int, x), dx_raw .+ range(0, stop=drift, length=lengthT)')
    dx_ci  = _dx_vec(dx_raw, lengthT)

    for inhibitidx in inhibitindices
        if inhibitidx > 0 && inhibitidx <= ncells
            T₂[:, inhibitidx] .*= -0.5
        end
    end

    img₂  = _composite_float(imgsz, template, centers, T₂; dx=dx_ci)
    bg    = maximum(img₂) * bias
    img₂  = img₂ .+ bg + 10^(-SNR/20) .* randn(size(img₂))
    img₂a = AxisArray(img₂, :x, :y, :time)
    gtW, gtH, gtWimgc = makegt(template, centers, T₂, imgsz, ncells, bg; gtincludebg)
    gtbg  = copy(bg)
    imgrs = Matrix(reshape(img₂a, prod(imgsz), nimages(img₂a)))
    ncells, imgrs, img₂, gtW, gtH, gtWimgc, gtbg
end

function gaussian_two_objs(sigma, imgsz::NTuple{2}, lengthT, distance, overlap_rate=0.5;
        fovsz::NTuple{2}=imgsz, jitter=0, drift=0, bias=0.1, SNR=10, gtincludebg=false)
    template = gaussiantemplate(Float64, (sigma, sigma))
    rcenter  = fovsz[1] ÷ 2
    ccenter  = (fovsz[2] - distance) ÷ 2
    centers  = [CartesianIndex(rcenter, ccenter), CartesianIndex(rcenter, fovsz[2]-ccenter)]
    ol = lengthT * overlap_rate
    l1 = round(Int, (lengthT + ol) / 2)
    l2 = lengthT - l1
    T₂ = reshape(Float64[rand(l1)..., zeros(l2)..., zeros(l2)..., rand(l1)...], lengthT, 2)

    dx_raw = rand(-jitter:jitter, 2, lengthT)
    dx_raw = map(x -> round(Int, x), dx_raw .+ range(0, stop=drift, length=lengthT)')
    dx_ci  = _dx_vec(dx_raw, lengthT)

    img₂      = _composite_float(imgsz, template, centers, T₂; dx=dx_ci)
    signalpwr = sum(img₂.^2) / length(img₂)
    bg        = maximum(img₂) * bias
    noise     = randn(size(img₂)) .* sqrt(signalpwr / 10^(SNR/10))
    img₂      = img₂ .+ bg + noise
    img₂a     = AxisArray(img₂, :x, :y, :time)
    gtW, gtH, gtWimgc = makegt(template, centers, T₂, imgsz, 2, bg; gtincludebg)
    gtbg  = copy(bg)
    imgrs = Matrix(reshape(img₂a, prod(imgsz), nimages(img₂a)))
    2, imgrs, img₂, gtW, gtH, gtWimgc, gtbg
end

function makegt(template::AbstractArray{Ts}, centers::Vector{CartesianIndex{2}},
                Ht, imgsz::NTuple{2}, ncells, bg; gtincludebg=true) where Ts
    numgtcomp = gtincludebg ? ncells + 1 : ncells
    gtimg = Array{Ts}(undef, imgsz..., numgtcomp)
    gtH   = Array{Ts}(undef, numgtcomp, size(Ht, 1))
    for i = 1:ncells
        intensity = zeros(Ts, 1, ncells)
        intensity[1, i] = one(Ts)
        img0 = _composite_float(imgsz, template, centers, intensity)[:,:,1]
        nrm  = sqrt(sum(img0.^2))
        img0 ./= nrm
        gtH[i, :]    = Ht[:, i] .* nrm
        gtimg[:,:,i] = img0
    end
    if gtincludebg
        fillval = bg == 0 ? zero(Ts) : one(Ts) / sqrt(prod(imgsz))
        fill!(view(gtimg, :, :, numgtcomp), fillval)
        fill!(view(gtH, numgtcomp, :), sqrt(prod(imgsz)) * bg)
    end
    gtimgrs = reshape(gtimg, prod(imgsz), numgtcomp)
    mxabs   = maximum(abs, gtimg)
    fsc     = scalesigned(mxabs)
    fcol    = colorsigned()
    gtimgrs, gtH, mappedarray(fcol ∘ fsc, gtimgrs)
end

nothing
