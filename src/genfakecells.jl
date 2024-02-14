function gaussian1D(ncells, lengthT, nevents, bias=0.1, SNR = 10) # should lengthT > 50
    S₁ = gaussiantemplate(Float64,2)
    mid = (size(S₁,2) + 1) ÷ 2
    S₁ = S₁[:,mid]  # 1d
    T₀ = makeintensity(lengthT, ncells, nevents)
    T₀[37,1:2] = [0.2, 0.7]   # to ensure at least one frame when both cells are active (→ SVD has negative components)
    if ncells > 3
        T₀[50,3:4] = [0.2, 0.7]
    elseif ncells > 2
        T₀[20,2:3] = [0.2, 0.7]
    end
    S₀ = zeros(eltype(S₁), 2*ncells*mid, ncells)
    for i = 1:ncells
        S₀[(i-1)*length(S₁)+1:(i-1)*length(S₁)+length(S₁), i] = S₁
    end
    img = S₀*T₀'
    signalpwr = sum(img.^2)/length(img)
    bg = randn(size(img)).*sqrt(signalpwr/10^(SNR/10)) .+ maximum(img₂)*bias
    img = img + bg# add noise and bg
    img, S₀, T₀', bg
end

# bias is a percentage for the maximum(img₂)
function gaussian2D(sigma, imgsz::NTuple{2}, lengthT, revent=10; fovsz::NTuple{2}=imgsz, jitter=0,
        drift = 0, bias=0.1, SNR = 10, useCalciumT=false, orthogonal=true, overlaplevel=1,
        inhibitindices=0, gtincludebg=false) # should lengthT > 50
    S₂ = gaussiantemplate(Float64, sigma)
    centers = spreadcells(fovsz, sigma)
    imagecenter = fovsz.÷2
    overlap_cell_center = [imagecenter[1]-min(overlaplevel,imagecenter[1])÷1.5,imagecenter[2]+min(overlaplevel,imagecenter[2])] # to add an overlapped cell
    !orthogonal && (centers = [centers overlap_cell_center])
    ncells = size(centers, 2)
    nevents = revent*ncells*lengthT/100
    if useCalciumT
        fr = makefiringrate(lengthT, ncells, lambda=0.01) # rate parameter λ in the Poisson distribution
        # (firingrate::Array{T1,2}, arcoeff::Vector{T2}, noisesigma, alpha, baseline
        T₂ = calcium_transient(fr, [0.85], 0., 1, 0.) # add noise and baseline later
    else
        T₂ = makeintensity(lengthT, ncells, nevents) # nevents=10*ncells
    end
    dx = rand(-jitter:jitter,2,lengthT)   # simulate some registration jitter
    dx = map(x->round(Int, x), dx .+ range(0, stop=drift, length=lengthT)')
    if !isempty(inhibitindices)
        for inhibitidx in inhibitindices
            if inhibitidx > 0 && inhibitidx <= ncells
                T₂[:,inhibitidx] .*= 0.5
                T₂[:,inhibitidx] .*= -1 # inverted
            end
        end
        img_wobias = makeimages(imgsz..., S₂, centers, T₂, dx)
        signalpwr = sum(img_wobias.^2)
        for inhibitidx in inhibitindices
            if inhibitidx > 0 && inhibitidx <= ncells
                biasinhibit = -minimum(T₂[:,inhibitidx])
                T₂[:,inhibitidx] .+= biasinhibit # add individual bias to inverted T₂
            end
        end
    end
    img₂ = makeimages(imgsz..., S₂, centers, T₂, dx)
    signalpwr = isempty(inhibitindices) ? sum(img₂.^2) : signalpwr
    bg = maximum(img₂)*bias
    noise = randn(size(img₂)); noisepwr0 = sum(noise.^2)
    noiseamp = sqrt(signalpwr/10^(SNR/10)/noisepwr0); noise .*= noiseamp
    img₂ = img₂ .+ bg + noise # add noise and bg
    img₂a = AxisArray(img₂, :x, :y, :time)
    gtW, gtH, gtWimgc = makegt(S₂,centers,T₂,imgsz,ncells,bg; gtincludebg=gtincludebg)
    gtbg = copy(bg)
    imgrs = Matrix(reshape(img₂a, prod(imgsz), nimages(img₂a)))
    ncells, imgrs, img₂, gtW, gtH, gtWimgc, gtbg
end

function gaussian2D_AmplitudeSNR(sigma, imgsz::NTuple{2}, lengthT, revent=10; fovsz::NTuple{2}=imgsz, jitter=0,
        drift = 0, bias=0.1, SNR = 10, useCalciumT=false, orthogonal=true, overlaplevel=1,
        inhibitindices=0, gtincludebg=false) # should lengthT > 50
    S₂ = gaussiantemplate(Float64, sigma)
    centers = spreadcells(fovsz, sigma)
    imagecenter = fovsz.÷2
    overlap_cell_center = [imagecenter[1]-min(overlaplevel,imagecenter[1])÷1.5,imagecenter[2]+min(overlaplevel,imagecenter[2])] # to add an overlapped cell
    !orthogonal && (centers = [centers overlap_cell_center])
    ncells = size(centers, 2)
    nevents = revent*ncells*lengthT/100
    if useCalciumT
        fr = makefiringrate(lengthT, ncells, lambda=0.01) # rate parameter λ in the Poisson distribution
        # (firingrate::Array{T1,2}, arcoeff::Vector{T2}, noisesigma, alpha, baseline
        T₂ = calcium_transient(fr, [0.85], 0., 1, 0.) # add noise and baseline later
    else
        T₂ = makeintensity(lengthT, ncells, nevents) # nevents=10*ncells
    end
    dx = rand(-jitter:jitter,2,lengthT)   # simulate some registration jitter
    dx = map(x->round(Int, x), dx .+ range(0, stop=drift, length=lengthT)')
    if !isempty(inhibitindices)
        for inhibitidx in inhibitindices
            if inhibitidx > 0 && inhibitidx <= ncells
                T₂[:,inhibitidx] .*= 0.5
                T₂[:,inhibitidx] .*= -1
            end
        end
    end
    img₂ = makeimages(imgsz..., S₂, centers, T₂, dx)
    bg = maximum(img₂)*bias
    noise = 10^(-SNR/20)*randn(size(img₂))
    img₂ = img₂ .+ bg + noise # add noise and bg
    img₂a = AxisArray(img₂, :x, :y, :time)
    gtW, gtH, gtWimgc = makegt(S₂,centers,T₂,imgsz,ncells,bg; gtincludebg=gtincludebg)
    gtbg = copy(bg)
    imgrs = Matrix(reshape(img₂a, prod(imgsz), nimages(img₂a)))
    ncells, imgrs, img₂, gtW, gtH, gtWimgc, gtbg
end

function gaussian_two_objs(sigma, imgsz::NTuple{2}, lengthT, distance, overlap_rate = 0.5; fovsz::NTuple{2}=imgsz,
        jitter=0, drift = 0, bias=0.1, SNR = 10, gtincludebg=false) # should lengthT > 50
    S₂ = gaussiantemplate(Float64, sigma)
    rcenter = fovsz[1]÷2; ccenter = (fovsz[2]-distance)÷2
    centers = Float64[rcenter ccenter; rcenter fovsz[2]-ccenter]
    ol = lengthT*overlap_rate; l1 = Int(round((lengthT+ol)/2)); l2 = lengthT-l1
    T₂ = reshape(Float64[rand(l1)..., zeros(l2)..., zeros(l2)..., rand(l1)...],lengthT,2)
    dx = rand(-jitter:jitter,2,lengthT)   # simulate some registration jitter
    dx = map(x->round(Int, x), dx .+ range(0, stop=drift, length=lengthT)')
    img₂ = makeimages(imgsz..., S₂, centers, T₂, dx)
    signalpwr = sum(img₂.^2)/length(img₂)
    bg = maximum(img₂)*bias
    noise = randn(size(img₂)).*sqrt(signalpwr/10^(SNR/10))
    img₂ = img₂ .+ bg + noise# add noise and bg
    img₂a = AxisArray(img₂, :x, :y, :time)
    gtW, gtH, gtWimgc = makegt(S₂,centers,T₂,imgsz,2,bg,gtincludebg=gtincludebg)
    gtbg = copy(bg)
    imgrs = Matrix(reshape(img₂a, prod(imgsz), nimages(img₂a)))
    2, imgrs, img₂, gtW, gtH, gtWimgc, gtbg
end

# function gaussian2D_calcium_transient(sigma, imgsz::NTuple{2}, lengthT, revent=10; fovsz::NTuple{2}=imgsz,
#         jitter=0, drift = 0, lambda=0.01, decayconst=0.85, bias=0.1, SNR = 10, overlaplevel=1) # should lengthT > 50
#     S₂ = gaussiantemplate(Float64, sigma)
#     centers = spreadcells(fovsz, sigma)
#     imagecenter = fovsz.÷2
#     overlap_cell_center = [imagecenter[1]-min(overlaplevel,imagecenter[1])÷1.5,imagecenter[2]+min(overlaplevel,imagecenter[2])] # to add an overlapped cell
#     centers = [centers overlap_cell_center]
#     ncells = size(centers, 2)
#     nevents = revent*ncells*lengthT/100
#     if lambda == 0
#         firingrate = zeros(UInt, lengthT, ncells)
#         offset = 10
#         step = max(offset*ncells, Int(round(lengthT/nevents*ncells)))
#         for (i,ofst) in enumerate(1:offset:offset*ncells)
#             for j in ofst:step:lengthT
#                 firingrate[j,i] = 1.0
#             end
#         end
#     else
#         firingrate = makefiringrate(lengthT, ncells, lambda = lambda)
#     end
#     T₂  = calcium_transient(firingrate, [decayconst], 0., 1, 0.)
#     dx = rand(-jitter:jitter,2,lengthT)   # simulate some registration jitter
#     dx = map(x->round(Int, x), dx .+ range(0, stop=drift, length=lengthT)')
#     img₂ = makeimages(imgsz..., S₂, centers, T₂, dx)
#     signalpwr = sum(img₂.^2)/length(img₂)
#     bg = randn(size(img₂)).*sqrt(signalpwr/10^(SNR/10)) .+ maximum(img₂)*bias
#     img₂ = img₂ + bg# add noise and bg
#     img₂a = AxisArray(img₂, :x, :y, :time)
#     gtW, gtH, gtWimgc = makegt(S₂,centers,T₂,imgsz,ncells)
#     gtbg = copy(bg)
#     imgrs = Matrix(reshape(img₂a, prod(imgsz), nimages(img₂a)))
#     ncells, imgrs, img₂, gtW, gtH, gtWimgc, gtbg
# end

function makegt(S₂::AbstractArray{Ts},centers,H,imgsz,ncells,bg; gtincludebg=true) where Ts
    numgtcomp = gtincludebg ? ncells+1 : ncells
    gtimg = Array{Ts}(undef,imgsz...,numgtcomp) # ncells + bg
    gtH = Array{Ts}(undef,size(H,1),numgtcomp)
    for i = 1:ncells
        T = zeros(1,ncells) # total ncells timepoints
        T[i] = 1.0
        img0 = makeimages(imgsz..., S₂, centers, T) # ith cell is activated at timepoint 1.
                                                    # This means length 1 movie X is the ith cell shape
        nrm = sqrt(sum(img0.^2))
        img0 ./= nrm # normalized
        gtH[:,i] = H[:,i] .* nrm                    # power is all stored to H component
        gtimg[:,:,i] = dropdims(img0,dims=3)        # each gtimg[:,:,i] holds the ith cell shape
    end
    gtincludebg && begin
        fillval = bg == 0 ? zero(Ts) : 1/sqrt(*(imgsz...))
        fill!(view(gtimg,:,:,numgtcomp), fillval)
        fill!(view(gtH,:,numgtcomp), sqrt(*(imgsz...))*bg)
    end
    gtimgrs = reshape(gtimg  , prod(imgsz), size(gtimg)[end])
    mxabs = maximum(abs, gtimg)
    fsc = scalesigned(mxabs)
    fcol = colorsigned()
    gtimgrs, gtH, mappedarray(fcol ∘ fsc, reshape(gtimg, Val(2)))
end

nothing