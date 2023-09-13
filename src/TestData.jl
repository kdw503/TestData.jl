module TestData

using MAT, GLMakie, Colors, JLD2, Printf, ImageView, Images
using FakeCells, AxisArrays, ImageCore, MappedArrays, DataStructures
using ImageAxes # avoid using ImageCore.nimages for AxisArray type array
include("genfakecells.jl")

export load_data, imsave_data, plot_convergence, plotWH_data, plotH_data

if Sys.iswindows()
    datapath="C:\\Users\\kdw76\\WUSTL\\Work\\Data"
elseif Sys.isunix()
    datapath=ENV["MYSTORAGE"]*"/work/Data"
end

set_datapath(path) = (datapath = path)

function load_cbcl()
    dirpath = joinpath(datapath,"MIT-CBCL-face","face.train","train","face")
    nRow = 19; nCol = 19; nFace = 2429; imgsz = (nRow, nCol)
    X = zeros(nRow*nCol,nFace)
    for i in 1:nFace
        fname = "face"*@sprintf("%05d",i)*".pgm"
        img = load(joinpath(dirpath,fname))
        X[:,i] = vec(img)
    end
    ncells = 49
    X, imgsz, nFace, ncells, 0, Dict()
end

 # :orlface : 112X92, lengthT = 400, ncomponents = 49; 1iter:10min(convex)
function load_orl()
    dirpath = joinpath(datapath,"ORL_Face")
    nRow = 112; nCol = 92; nPeople = 40; nFace = 10; imgsz = (nRow, nCol); lengthT = nPeople*nFace
    X = zeros(nRow*nCol,lengthT)
    for i in 1:nPeople
        subdirpath = "s$i"
        for j in 1:nFace
            fname = "$j.bmp"
            img = load(joinpath(dirpath,subdirpath,fname))
            X[:,i] = vec(Gray.(img))
        end
    end
    ncells = 25
    X, imgsz, lengthT, ncells, 0, Dict()
end

# :natural : 12X12, lengthT=100000, ncomponents = 72
function load_natural(;num_trials=1000, batch_size=100, sz=12, BUFF=4) # BUFF(border margin)
    dirpath = joinpath(datapath,"natural")
    dd = matread(joinpath(dirpath,"IMAGES.mat")); img = dd["IMAGES"]
    
    lengthT = num_trials*batch_size
    imgsz = (sz, sz); l=sz^2
    num_images=size(img,3); image_size=512
    X=zeros(l,lengthT);
    for t=1:num_trials
        # choose an image for this batch
        i=Int(ceil(num_images*rand()));
        this_image=img[:,:,i] #reshape(IMAGES(:,i),image_size,image_size)';
        # extract subimages at random from this image to make data vector X(64 X batch_size)
        for i=1:batch_size
            r=Int(BUFF+ceil((image_size-sz-2*BUFF)*rand()));
            c=Int(BUFF+ceil((image_size-sz-2*BUFF)*rand()));
            patchp = this_image[r:r+sz-1,c:c+sz-1]
            X[:,(t-1)*batch_size+i]=reshape(patchp,l,1)
        end
    end
    ncells=72
    X, imgsz, lengthT, ncells, 0, Dict()
end

# :natural : 12X12, lengthT=100000, ncomponents = 72
function load_onoffnatural(;num_trials=1000, batch_size=100, sz=12, BUFF=4) # BUFF(border margin)
    dirpath = joinpath(datapath,"natural")
    dd = matread(joinpath(dirpath,"IMAGES.mat")); img = dd["IMAGES"]

    lengthT = num_trials*batch_size
    imgsz = (sz, sz); l=sz^2
    num_images=size(img,3); image_size=512
    X=zeros(2*l,lengthT);
    for t=1:num_trials
        # choose an image for this batch
        i=Int(ceil(num_images*rand()));
        this_image=img[:,:,i] #reshape(IMAGES(:,i),image_size,image_size)';
        # extract subimages at random from this image to make data vector X(64 X batch_size)
        for i=1:batch_size
            r=Int(BUFF+ceil((image_size-sz-2*BUFF)*rand()));
            c=Int(BUFF+ceil((image_size-sz-2*BUFF)*rand()));
            patchp = this_image[r:r+sz-1,c:c+sz-1]
            patchm = -this_image[r:r+sz-1,c:c+sz-1]
            patchp[patchp.<0].=0; patchm[patchm.<0].=0
            X[:,(t-1)*batch_size+i]=[reshape(patchp,l,1); reshape(patchm,l,1)]
        end
    end
    ncells=72
    X, imgsz, lengthT, ncells, 0, Dict()
end

function load_neurofinder()
    roi = (71:170, 31:130)
    dirpath = joinpath(datapath,"neurofinder")
    fullorgimg = load(joinpath(dirpath,"neurofinder.02.00.cut100250_sqrt.tif"))
    orgimg = fullorgimg[roi...,:]
    # save("neurofinder.02.00.cut.gif",orgimg.*2)
    imgsz = size(orgimg)[1:end-1]; lengthT = size(orgimg)[end]
    X = Matrix{Float64}(reshape(orgimg, prod(imgsz), lengthT))
    ncells = 40; gtncells = 32
    X, imgsz, lengthT, ncells, gtncells, Dict("imgrs"=>orgimg)
end

function load_neurofinder_small()
    dirpath = joinpath(datapath,"neurofinder")
    orgimg = load(joinpath(dirpath,"neurofinder.02.00.cut100250_small.tif"))
    dcells = load(joinpath(dirpath,"neurofinder.02.00.cut100250_small_gt.jld2"))
    # save("neurofinder.02.00.cut.gif",orgimg.*2)
    imgsz = size(orgimg)[1:end-1]; lengthT = size(orgimg)[end]
    X = Matrix{Float64}(reshape(orgimg, prod(imgsz), lengthT))
    ncells = 30; gtncells = 16
    X, imgsz, lengthT, ncells, gtncells, Dict("imgrs"=>orgimg, "cells"=>dcells["cells"])
end

function load_urban()
    dirpath = joinpath(datapath,"Hyperspectral-Urban")
    vars = matread(joinpath(dirpath,"Urban_R162.mat"))
    Y=vars["Y"]; nRow = Int(vars["nRow"]); nCol = Int(vars["nCol"]); nBand = 162; imgsz = (nRow, nCol)
    maxvalue = maximum(Y)
    X = Array(Y')./maxvalue;  # reinterpret(N0f8,UInt8(255))=1.0N0f8 but vars["Y"] has maximum 1000
    # img = reshape(X, (nRow, nCol, nBand));
    gtncells = 7; ncells = gtncells + 5
    X, imgsz, nBand, ncells, gtncells, Dict()
end

function load_audio()
    dirpath = joinpath(datapath,"audio")
    vars = load(joinpath(dirpath,"Maryhadalittlelamb.jld2")) # 257 frequency bins X 295 frames
    X=vars["Maryhadalittlelamb"]; lengthT = size(X)[end]; imgsz = size(X)
    # img = reshape(X, (nRow, nCol, nBand));
    gtncells = 3; ncells = gtncells
    X, imgsz, lengthT, ncells, gtncells, Dict()
end

function load_inhibit_real()
    dirpath = joinpath(datapath,"inhibit")

    orgimg = load(joinpath(dirpath,"Slice22_Fish201712127_ELO_ERO_Both_HS.tif")) # (540, 640, 420)
    st = (161,181); ed = (270,320) # inhibit neuron location (253,282) -> in the cut image (93,101)
    orgimgcut = orgimg[st[1]:ed[1], st[2]:ed[2], :]; sz = size(orgimgcut)
    fixed_index=200; margin=10; mxshift = (margin,margin); mxrot=(0.5,)
    cimg = orgimgcut[margin+1:sz[1]-margin,margin+1:sz[2]-margin,:]
    # rimg_trans = register(orgimgcut, fixed_index, mxshift, mxrot, margin; method=:rigid)
    # img2rigid = reshape(vcat(reshape(cimg,10800,420),reshape(rimg_rigid,10800,420)),90,240,420)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut.tif"),cimg)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_rigid.tif"),rimg_rigid)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_rigid_compare.tif"),img2rigid)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_rigid_mean.tif"),mean(img2rigid,dims=3))
    # # img2rigid = load(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_rigid_compare.tif"))
    # rimg_affine = register(orgimgcut, fixed_index, mxshift, mxrot, margin; method=:affine)
    # img2affine = reshape(vcat(reshape(cimg,10800,420),reshape(rimg_affine,10800,420)),90,240,420)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_affine.tif"),rimg_affine)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_affine_compare.tif"),img2affine)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_affine_mean.tif"),mean(img2affine,dims=3))
    # img2both = reshape(vcat(reshape(rimg_rigid,10800,420),reshape(rimg_affine,10800,420)),90,240,420)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_both_compare.tif"),img2both)
    # rimg_trans = register(orgimgcut, fixed_index, mxshift, mxrot, margin; method=:translate)
    # img2trans = reshape(vcat(reshape(cimg,10800,420),reshape(rimg_trans,10800,420)),90,240,420)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_trans.tif"),rimg_trans)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_trans_compare.tif"),img2trans)
    # Images.save(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_trans_mean.tif"),mean(img2trans,dims=3))

    img = load(joinpath(datapath,"inhibit","Slice22_Fish201712127_ELO_ERO_Both_HS_cut_trans.tif"))[:,:,51:end]
    (nRow, nCol, lengthT) = size(img); imgsz = (nRow, nCol) # activated neuron (38,21) inhibited neuron location (82,90)
    X=Float32.(reshape(img,nRow*nCol,lengthT))
    # img = reshape(X, (nRow, nCol, nBand));
    gtncells = 150; ncells = gtncells
    X, imgsz, lengthT, ncells, gtncells, Dict("activated_loc"=>(38,21),"inhibited_loc"=>(82,90))
end

function register(imgs,fixed_index,mxshift, mxrot, margin; method=:rigid, presmoothed=false, SD=I, initial_tfm=RegisterQD.IdentityTransformation(), kwargs...)
    slicedim = ndims(imgs); sz = size(imgs)
    fixed = imgs[:,:,fixed_index]
    minwidth_rot = fill(0.002, 3) # specifies the lower limit of resolution for the rotation
    rimgs = Array{eltype(imgs)}(undef,sz[1]-2*margin,sz[2]-2*margin,sz[3])
    for (i,moving) in enumerate(eachslice(imgs,dims=slicedim))
        @show i
        if method == :rigid
            tform, mm = qd_rigid(fixed, moving, mxshift, mxrot; presmoothed=presmoothed, SD=SD,
                    minwidth_rot=minwidth_rot, initial_tfm=initial_tfm, maxevals=1000, rtol=0, fvalue=0.0002, kwargs...)
        elseif method == :affine
            tform, mm = qd_affine(fixed, moving, mxshift; SD = SD, maxevals=1000, rtol=0, fvalue=0.0002)
        elseif method == :translate
            @show "translate"
            tform, mm = qd_translate(fixed, moving, mxshift; maxevals=1000, rtol=0, fvalue=0.0003)
        end
        wmoving = warp(moving,tform)
        # inds = intersect.(Base.axes(fixed[margin+1:end-margin]), Base.axes(moving))
        rimgs[:,:,i] .= wmoving[margin+1:sz[1]-margin,margin+1:sz[2]-margin]
    end
    rimgs
end

function load_fakecells(;SNR=10, user_ncells=0, imgsz=(40,20), fovsz=imgsz, lengthT=1000, bias=0.1, useCalciumT=false,
        jitter=0, inhibitindices=0, gtincludebg=false, only2cells=false, issave=true, isload=true, save_maxSNR_X=false,
        save_X=false, save_gtimg=false)
    dirpath = joinpath(datapath,"fakecells")
    calciumstr = useCalciumT ? "_calcium" : ""
    fprefix = "fakecells$(inhibitindices)$(calciumstr)_sz$(imgsz)_lengthT$(lengthT)_J$(jitter)_SNR$(SNR)_bias$(bias)"
    dfprefix = joinpath(dirpath,fprefix)
    X, imgsz, fakecells_dic, img_nl, maxSNR_X = loadfakecell(Float64, dfprefix*".jld2", only2cells=only2cells,
        fovsz=fovsz, imgsz=imgsz, lengthT=lengthT, bias=bias, useCalciumT=useCalciumT, jitter=jitter, SNR=SNR,
        inhibitindices=inhibitindices, gtincludebg=gtincludebg, issave=issave, isload=isload);
    gtncells = fakecells_dic["gt_ncells"]
    if save_gtimg
        gtW = fakecells_dic["gtW"]; gtH = fakecells_dic["gtH"]
        W3,H3 = copy(gtW), copy(gtH')
        imsaveW(dfprefix*"_GT_W.png", W3, imgsz, borderwidth=1)
        imsaveH(dfprefix*"_GT_H.png", H3, 100, colors=g1wm())
    end
    if save_maxSNR_X
        imsaveW(dfprefix*"_maxSNR_W.png", maxSNR_X, imgsz, borderwidth=1,colors=bbw())
    end
    if save_X
        options = (crf=23, preset="medium")
        clamp_level=1.0; X_max = maximum(abs,X)*clamp_level; Xnor = X./X_max;  X_clamped = clamp.(Xnor,0.,1.)
        Xuint8 = UInt8.(round.(map(clamp01nan, X_clamped)*255))
        VideoIO.save(dfprefix*".mp4", reshape.(eachcol(Xuint8),imgsz...), framerate=30, encoder_options=options)
    end
    ncells = user_ncells==0 ? gtncells + 8 : user_ncells
    X, imgsz, lengthT, ncells, gtncells, fakecells_dic
end

# TODO: load and save sometimes cause error. Do I need to change to below fucntions?
# jldopen("mydata.jld", "w") do file
#     write(file, "A", A)  # alternatively, say "@write file A"
# end
# c = jldopen("mydata.jld", "r") do file
#     read(file, "A")
# end
function loadfakecell(T::Type, fname; lengthT=100, imgsz=(40,20), fovsz=imgsz, SNR=10, bias=0.1,
        useCalciumT=false, jitter=0, distance = 10, only2cells=false, overlap_rate = 0.3,
        inhibitindices=0, gtincludebg=true, issave=true, isload=true)
   if isload && isfile(fname)
        fakecells_dic = JLD2.load(fname)
#        JLD.load(fname)
        gt_ncells = fakecells_dic["gt_ncells"]
        imgrs = fakecells_dic["imgrs"]
        img_nl = fakecells_dic["img_nl"]
        gtW = fakecells_dic["gtW"]
        gtH = fakecells_dic["gtH"]
        gtWimgc = fakecells_dic["gtWimgc"]
        gtbg = fakecells_dic["gtbg"]
        imgsz = fakecells_dic["imgsz"]
        SNR = fakecells_dic["SNR"]
    else
        @show fname
        isload && @warn "$fname not found. Generating fakecell data..."
        sigma = 5; imgsz = imgsz; lengthT = lengthT
        if only2cells
            gt_ncells, imgrs, img_nl, gtW, gtH, gtWimgc, gtbg = gaussian_two_objs(sigma, imgsz, lengthT,
                                    distance, overlap_rate; jitter=jitter, fovsz=fovsz, SNR = SNR)
        else
            revent = 10
            gt_ncells, imgrs, img_nl, gtW, gtH, gtWimgc, gtbg = gaussian2D(sigma, imgsz, lengthT, revent,
                bias=bias, useCalciumT=useCalciumT,jitter=jitter, fovsz=fovsz, SNR=SNR, orthogonal=false,
                inhibitindices=inhibitindices,gtincludebg=gtincludebg)
        end
        fakecells_dic = Dict()
        fakecells_dic["gt_ncells"] = gt_ncells
        fakecells_dic["imgrs"] = imgrs
        fakecells_dic["img_nl"] = img_nl
        fakecells_dic["gtW"] = gtW
        fakecells_dic["gtH"] = gtH
        fakecells_dic["gtWimgc"] = Array(gtWimgc)
        fakecells_dic["gtbg"] = gtbg
        fakecells_dic["imgsz"] = imgsz
        fakecells_dic["SNR"] = SNR
        if issave
            JLD2.save(fname,fakecells_dic)
        end
    end
    X = T.(imgrs)
    maxindices = argmax.(eachcol(gtH))
    maxSNR_X = X[:,[maxindices...]]

    return X, imgsz, fakecells_dic, img_nl, maxSNR_X
end

function load_data(dataset; SNR=10, user_ncells=0, imgsz=(40,20), fovsz=imgsz, lengthT=1000, useCalciumT=false,
        jitter=0, bias=0.1, inhibitindices=0, gtincludebg=false, issave=true, isload=true, save_maxSNR_X=false,
        save_X=false, save_gtimg=false)
    if dataset == :cbclface
        println("loading CBCL face dataset")
        load_cbcl()
    elseif dataset == :orlface
        println("loading ORL face dataset")
        load_orl()
    elseif dataset == :natural
        println("loading natural image")
        load_natural()
    elseif dataset == :onoffnatural
        println("loading On/OFF-contrast filtered natural image")
        load_onoffnatural()
    elseif dataset == :urban
        println("loading hyperspectral urban dataset")
        load_urban()
    elseif dataset == :audio
        println("loading maryhadalittlelamb audio dataset")
        load_audio()
    elseif dataset == :neurofinder
        println("loading Neurofinder dataset")
        load_neurofinder()
    elseif dataset == :neurofinder_small
        println("loading Neurofinder small dataset")
        load_neurofinder_small()
    elseif dataset == :inhibit_real
        load_inhibit_real()
    elseif dataset == :fakecells
        println((isload ? "Loading" : "Generating") * " image of fakecells")
        load_fakecells(only2cells=false, SNR=SNR, user_ncells=user_ncells, imgsz=imgsz, fovsz=fovsz, lengthT=lengthT,
            bias=bias, useCalciumT=useCalciumT, jitter=jitter, inhibitindices=inhibitindices, gtincludebg=gtincludebg,
            issave=issave, isload=isload, save_maxSNR_X=save_maxSNR_X, save_X = save_X, save_gtimg=save_gtimg)
    elseif dataset == :fakecellsmall
        println((isload ? "Loading" : "Generating") * "image of fakecells with only two cells")
        load_fakecells(only2cells=true, SNR=SNR, user_ncells=6, imgsz=(20,30), lengthT=lengthT, bias=bias,
            gtincludebg=gtincludebg, issave=issave, isload=isload, useCalciumT=useCalciumT,jitter=jitter,
            save_maxSNR_X=save_maxSNR_X, save_gtimg=save_gtimg)
    else
        error("Not supported dataset")
    end
end

#======== Image Save ==========================================================#

function imsave_data(dataset,fprefix,W,H,imgsz,lengthT; mssdwstr="", mssdhstr="",
        signedcolors=nothing,gridcols=nothing,saveH=true)
    if dataset == :cbclface
        println("Saving image of CBCL face dataset")
        imsave_cbcl(fprefix,W,H,imgsz,lengthT; signedcolors=signedcolors)
    elseif dataset == :orlface
        println("Saving image of ORL face dataset")
        imsave_orl(fprefix,W,H,imgsz,lengthT; signedcolors=signedcolors)
    elseif dataset == :natural
        println("Saving image of natural image")
        imsave_natural(fprefix,W,H,imgsz,lengthT; signedcolors=signedcolors)
    elseif dataset == :onoffnatural
        println("Saving image of On/OFF-contrast filtered natural image")
        imsave_onoffnatural(fprefix,W,H,imgsz,lengthT; signedcolors=signedcolors)
    elseif dataset == :urban
        println("Saving hyperspectral urban dataset")
        imsave_urban(fprefix,W[:,2:7],H[2:7,:],imgsz,lengthT; signedcolors=signedcolors)
    elseif dataset == :audio
        @warn "$(dataset) dataset doesn't have a image save method"
    elseif dataset ∈ [:neurofinder, :neurofinder_small, :inhibit_real]
        println("Saving image of Neurofinder dataset")
        imsave_neurofinder(fprefix,W,H,imgsz,lengthT, signedcolors=signedcolors, gridcols=gridcols, saveH=saveH)
    elseif dataset == :fakecells
        println("Saving image of fakecells")
        imsave_fakecell(fprefix,W,H,imgsz,lengthT; mssdwstr=mssdwstr, mssdhstr=mssdhstr,
                        signedcolors=signedcolors, saveH=saveH)
    elseif dataset == :fakecellsmall
        println("loading fakecells with only two cells")
        imsave_fakecell(fprefix,W,H,imgsz,lengthT; mssdwstr=mssdwstr, mssdhstr=mssdhstr,
                        signedcolors=signedcolors, saveH=saveH)
    end
end

dgwm() = (colorant"darkgreen", colorant"white", colorant"magenta")
dgwdm() = (colorant"darkgreen", colorant"white", colorant"darkmagenta")
bwm() = (colorant"blue", colorant"white", colorant"magenta")
g1wm() = (colorant"green1", colorant"white", colorant"magenta")
g1bw() = (colorant"green1", colorant"black", colorant"white")
bbw() = (colorant"black", colorant"black", colorant"white")
bgw() = (colorant"black", colorant"gray", colorant"white")
wwb() = (colorant"white", colorant"white", colorant"black") # cdcl-face NMF

function interpolate(img,n)
    imginter = zeros(size(img).*(n+1))
    offset = CartesianIndex(1,1)
    for i in CartesianIndices(img)
        ii = (i-offset)*n+i
        for j in CartesianIndices((n+1,n+1))
            imginter[ii+j-offset] = img[i]
        end
    end
    imginter
end

function imsave_cbcl(fprefix,W,H,imgsz,tlength; gridcols=Int(ceil(sqrt(size(W,2)))), borderwidth=1,
        signedcolors=nothing)
    signedcolors = signedcolors === nothing ? g1wm() : signedcolors
    # clamp_level=0.5; W3_max = maximum(abs,W3)*clamp_level; W3_clamped = clamp.(W3,0.,W3_max)
    # signedcolors = (colorant"green1", colorant"white", colorant"magenta")
    # imsaveW(fprefix*"_iter$(iter)_rt$(rt2).png", W3, imgsz, gridcols=7, colors=signedcolors,
    #       borderval=0.5, borderwidth=1)
    imsaveW(fprefix*"_W.png", W, imgsz; gridcols=gridcols, borderwidth=borderwidth,
        colors=signedcolors)
end

function imsave_reconstruct(fprefix,X,W,H,imgsz; index=100, gridcols=7, clamp_level=1.0) # 0.2 for urban
    # TODO: save all face pictures
    # imsave("cbcl_faces.png", X[:,1:49], imgsz;; gridcols=7, borderwidth=1, colors=bbw())
    intplimg = interpolate(reshape(H[:,index],gridcols,gridcols),3)
    imsaveW(fprefix*"_H$(index).png",reshape(intplimg,length(intplimg),1),size(intplimg))
    reconimg = W*H[:,index]
    W_max = maximum(abs,reconimg)*clamp_level; W_clamped = clamp.(reconimg,0.,W_max)
    #save(fprefix*"_Org$index.png", map(clamp01nan, reshape(X[:,index],imgsz...)))
    save(fprefix*"_recon$index.png", map(clamp01nan, reshape(W_clamped,imgsz...)))
end

function imsave_orl(fprefix,W,H,imgsz,tlength; gridcols=Int(ceil(sqrt(size(W,2)))),
        borderwidth=1, signedcolors=nothing)
    signedcolors = signedcolors === nothing ? g1wm() : signedcolors
    imsaveW(fprefix*"_W.png", W, imgsz; gridcols=gridcols, borderwidth=borderwidth,
        colors=signedcolors)
end

function imsave_natural(fprefix,W,H,imgsz,tlength; gridcols=12, borderwidth=1,
        signedcolors=nothing)
    signedcolors = signedcolors === nothing ? bgw() : signedcolors
    imsaveW(fprefix*"_W.png", W, imgsz; gridcols=gridcols, borderwidth=borderwidth,
            colors=signedcolors)
 end
 
function imsave_onoffnatural(fprefix,W,H,imgsz,tlength; gridcols=12, borderwidth=1,
        signedcolors=nothing)
    signedcolors = signedcolors === nothing ? bgw() : signedcolors
    l = imgsz[1]^2; W = W[1:l,:]-W[l+1:2*l,:]
    imsaveW(fprefix*"_W.png", W, imgsz; gridcols=gridcols, borderwidth=borderwidth,
            colors=signedcolors)
end

function imsave_urban(fprefix,W,H,imgsz,tlength; gridcols=2, borderwidth=5, signedcolors=nothing,
        clamp_level=0.2)
    signedcolors = signedcolors === nothing ? bbw() : signedcolors
    W_max = maximum(abs,W)*clamp_level; W_clamped = clamp.(W,-W_max,W_max)
    imsaveW(fprefix*"_W.png", W_clamped, imgsz; gridcols=gridcols,
            borderwidth=borderwidth, colors=signedcolors)
end

function imsave_reconstruct_urban(fprefix,X,W,H,imgsz; index=100, clamp_level=0.2)
    reconimg = W*H[:,index]
    W_max = maximum(abs,reconimg)*clamp_level; W_clamped = clamp.(reconimg,0.,W_max)
    save(fprefix*"_Org$index.png", map(clamp01nan, reshape(X[:,index],imgsz...)))
    save(fprefix*"_CG$index.png", map(clamp01nan, reshape(W_clamped,imgsz...)))
    #save(fprefix*"_CG_clamped$index.png", map(clamp01nan, reshape(clamp.(W,0.,Inf)*clamp.(H[:,index],0.,Inf),imgsz...)))
end

function imsave_neurofinder(fprefix,W,H,imgsz,tlength; gridcols=nothing, borderwidth=1,
        signedcolors=nothing, saveH=true)
    #@show ncells, gridcols
    signedcolors = signedcolors === nothing ? g1wm() : signedcolors
    gridcols = gridcols === nothing ? 5 : gridcols
    ncells = size(W,2)
    ncells < gridcols && @warn "ncells($(ncells)) should not be smaller than gridcols($(gridcols))"
    for i = 1:ncells÷gridcols
        imsaveW(fprefix*"_W_$(i).png",W[:,((i-1)*gridcols+1):(i*gridcols)],imgsz;
                gridcols=gridcols,colors=signedcolors,borderwidth=borderwidth)
    end
    saveH && imsaveH(fprefix*"_H.png", H, tlength; colors=signedcolors)
end

function imsave_fakecell(fprefix,W,H,imgsz,tlength; mssdwstr="", mssdhstr="",
        gridcols=size(W,2), borderwidth=1, signedcolors=nothing, saveH=true)
    signedcolors = signedcolors === nothing ? g1wm() : signedcolors
    imsaveW(fprefix*"_W$(mssdwstr).png", W, imgsz; gridcols=gridcols,
            borderwidth=borderwidth, colors=signedcolors)
    saveH && imsaveH(fprefix*"_H$(mssdhstr).png", H, tlength; colors=signedcolors)
end

function match_order_to_gt(W,H,gtW,gtH)
    mssd, ml, ssds = matchedWnssda(gtW,W); mssdH = ssdH(ml, gtH,H')
    neworder = zeros(Int,length(ml))
    for (gti, i) in ml
        neworder[gti]=i
    end
    p = size(W,2)
    for i in 1:p
        i ∉ neworder && push!(neworder,i)
    end
    copy(W[:,neworder]), copy(H[neworder,:]), mssd, mssdH
end

function imsave_data_gt(dataset,fprefix,W,H,gtW,gtH,imgsz,lengthT; signedcolors=nothing, saveH=true)
    signedcolors = signedcolors === nothing ? g1wm() : signedcolors
    Wno, Hno, mssd, mssdH = match_order_to_gt(W,H,gtW,gtH)
    imsave_data(dataset,fprefix,Wno,Hno,imgsz,lengthT; mssdwstr="_MSE"*@sprintf("%1.4f",mssd),
            mssdhstr="_MSE"*@sprintf("%1.4f",mssdH), signedcolors=signedcolors, saveH=saveH)
end
#=========== Plot Convergence ====================================================#
function plot_convergence(fprefix, x_abss, xw_abss, xh_abss, f_xs; title="")
    xlbl = "iteration"; ylbl = "log10(penalty)"
    lstrs = ["log10(f_x)", "log10(x_abs)", "log10(xw_abs)", "log10(xh_abs)"]
    logfx = log10.(f_xs); logxabs = log10.(x_abss);
    logxw = log10.(xw_abss); logxh = log10.(xh_abss)
    ax = plotW([logfx logxabs logxw logxh], fprefix*"_plot_all.png"; title=title,
            xlbl=xlbl, ylbl="", legendloc=1, arrange=:combine, legendstrs=lstrs)
            #,axis=[480,1000,-0.32,-0.28])
    ax = plotW(logfx, fprefix*"_plot_fx.png"; title=title, xlbl=xlbl, ylbl=ylbl,
            legendstrs = [lstrs[1]], legendloc=1)
    ax = plotW([logxabs logxw logxh], fprefix*"_plot_xabs.png"; title=title, xlbl=xlbl,
            ylbl="", legendloc=1, arrange=:combine, legendstrs = lstrs[2:4])
    ax
end

function plot_convergence(fprefix, x_abss, f_xs, f_x_abss=[]; title="")
    xlbl = "iteration"; ylbl = "log10(penalty)"
    lstrs = ["log10(f_x)", "log10(f_x_abs)", "log10(x_abs)"]
    logfx = log10.(f_xs); logxabs = log10.(x_abss);
    ax = plotW([logfx logxabs], fprefix*"_plot_all.png"; title=title,
            xlbl=xlbl, ylbl="", legendloc=1, arrange=:combine, legendstrs=lstrs)
            #,axis=[480,1000,-0.32,-0.28])
    if length(f_x_abss) > 0
        logfxabs = log10.(f_x_abss)
        ax = plotW([logfx logfxabs], fprefix*"_plot_fxfxabs.png"; title=title, xlbl=xlbl, ylbl=ylbl,
            legendstrs = [lstrs[1],lstrs[2]], legendloc=1)
    else
        ax = plotW(logfx, fprefix*"_plot_fx.png"; title=title, xlbl=xlbl, ylbl=ylbl,
            legendstrs = [lstrs[1]], legendloc=1)
    end
    ax = plotW(logxabs, fprefix*"_plot_xabs.png"; title=title, xlbl=xlbl, ylbl="",
            legendstrs = [lstrs[3]], legendloc=1)
    ax
end
#=========== Plot W and H ====================================================#
using GLMakie

function plotWH_data(dataset,fprefix,W,H; resolution = (800,400), space=1.0, issave=true,
            colors=distinguishable_colors(size(W,2); lchoices=range(0, stop=50, length=5)))
    if dataset == :audio
        plotWH_audio_data(fprefix,W,H; resolution=resolution, space=space, colors=colors, issave=issave)
    else
        error("Not supported for $(dataset) dataset")
    end
end

function plotWH_audio_data(fprefix,W,H; resolution = (800,400), space=1.0, title="",issave=true,
        colors=distinguishable_colors(size(W,2); lchoices=range(0, stop=50, length=5)))
    fig = GLMakie.Figure(resolution = resolution)
    fn = fprefix*"_plot_WH.png"
    ax1 = GLMakie.Axis(fig[1, 1], xlabel = "W column", ylabel = "Frequency (kHz)", xgridvisible=false,
                xticksvisible=false, xtickformat = "", ygridvisible=false, title = title)
    plotW_data(ax1, W, colors, rng=0:size(W,1)-1, scale=:linear, rotate=true, space=space)
    ax2 = GLMakie.Axis(fig[1, 2], xlabel = "Time[s]", ylabel = "Activations", xgridvisible=false,
                ygridvisible=false, yticksvisible=false, ytickformat="", title = title)
    plotW_data(ax2, H', colors, rng=0:size(H,2)-1, scale=:linear, rotate=false, space=space)
    issave && save(fn,fig,px_per_unit=2)
    fig
end

function plotW_data(ax::Makie.Axis, W::AbstractArray, colors::AbstractVector; rng=0:size(W,1)-1,
        scale=:linear, rotate=false, space=1.0)
    for (i,w) in enumerate(eachcol(W))
        spacei = (i-1)*space
        rotate ? lines!(ax,w.+spacei,rng,color=colors[i]) : lines!(ax,rng,w.+spacei,color=colors[i])
    end
end

function plotH_data(fprefix, H; resolution = (800,400), space=1.0, title="",issave=true,
        colors=distinguishable_colors(size(H,1); lchoices=range(0, stop=50, length=5)))
    fig = GLMakie.Figure(resolution = resolution)
    fn = fprefix*"_plot_H.png"
    ax2 = GLMakie.Axis(fig[1, 1], xlabel = "Time index", ylabel = "Intensity", xgridvisible=true,
                ygridvisible=true, yticksvisible=true, ytickformat="", title = title)
    plotW_data(ax2, H', colors, rng=0:size(H,2)-1, scale=:linear, rotate=false, space=space)
    issave && save(fn,fig,px_per_unit=2)
    fig
end

#=========== Plot H ====================================================#
function plotH_data(dataset, fprefix, H)
    if dataset == :urban
        f = plotH_urban(H[2:7,:]; titles = ["2","3","4","5","6","7"])
    else
        plotH_data_old(fprefix,H; arrange=:combine)
    end
    save(joinpath(subworkpath,fprefix*"_H.png"),f)
end

function plotH_urban(H; title="", titles=fill("",size(H,1)))
    n = size(H,2); rng = 0:n-1
    black=RGB{N0f8}(0.0,0.0,0.0)
    f = Figure(resolution = (900,1500))
    ax11=GLMakie.Axis(f[1,1],title=titles[1], titlesize=25)
    ax12=GLMakie.Axis(f[1,2],title=titles[2], titlesize=25)
    ax21=GLMakie.Axis(f[2,1],title=titles[3], titlesize=25)
    ax22=GLMakie.Axis(f[2,2],title=titles[4], titlesize=25)
    ax31=GLMakie.Axis(f[3,1],title=titles[5], titlesize=25)
    ax32=GLMakie.Axis(f[3,2],title=titles[6], titlesize=25)
    ln1 = lines!(ax11, rng, H[1,:], color=black)
    ln2 = lines!(ax12, rng, H[2,:], color=black)
    ln3 = lines!(ax21, rng, H[3,:], color=black)
    ln4 = lines!(ax22, rng, H[4,:], color=black)
    ln5 = lines!(ax31, rng, H[5,:], color=black)
    ln6 = lines!(ax32, rng, H[6,:], color=black)
    f
end

#========= Image show and save ===========================================#
is_ImageView_available = true
try
    Sys.islinux() && run(`ls /usr/bin/x11vnc`) # check if this is noVNC graphical platform
    using ImageView
    using Gtk.ShortNames
catch # not a graphical platform
    @warn("Not a RIS noVNC graphical platform")
    global is_ImageView_available = false
end

function mkimgW(W::Matrix{T},imgsz; gridcols=size(W,2), borderwidth=1, borderval=0.7, scalemtd=:maxwhole,
        colors=(colorant"green1", colorant"white", colorant"magenta")) where T
    ncells = size(W,2)
    if scalemtd == :maxwhole
        mxabs = max(eps(eltype(W)),maximum(abs, W))
        fsc = scalesigned(mxabs)
        fcol = colorsigned(colors...)
        Wcolor = Array(mappedarray(fcol ∘ fsc, reshape(W, Val(2))))
    elseif scalemtd == :maxcol
        fsc = (x) -> (mxabs=max(eps(eltype(x)),maximum(abs, x)); x./mxabs)
        sW = similar(W)
        for (i,x) in enumerate(eachcol(W)) sW[:,i] = fsc(x) end
        fcol = colorsigned(colors...)
        Wcolor = Array(mappedarray(fcol, reshape(sW, Val(2))))
    elseif scalemtd == :maxgridrow
        fsc = (x) -> (mxabs=max(eps(eltype(x)),maximum(abs, x)); x./mxabs)
        sW = similar(W)
        for i in 1:gridcols:ncells
            x = view(W,:,i:min(i+4,ncells))
            sW[:,i:min(i+4,ncells)]=fsc(x)
        end
        fcol = colorsigned(colors...)
        Wcolor = Array(mappedarray(fcol, reshape(sW, Val(2))))
    elseif scalemtd == :avgwhole
        avgabs = sum(abs,W)/length(W)*18; sW = W./avgabs; clamp!(sW,-1,1)
        fcol = colorsigned(colors...)
        Wcolor = Array(mappedarray(fcol, reshape(sW, Val(2))))
    end
    gridsz = ((ncells-1)÷gridcols+1,gridcols)
    add_dim_sz = ntuple(i->1,Val(length(imgsz)-length(gridsz)))
    bordersz = ntuple(i->borderwidth,Val(2))
    bimgsz = imgsz.+(bordersz..., add_dim_sz...)
    gimgsz = bimgsz.*(gridsz..., add_dim_sz...).+bordersz
    fill_val = eltype(Wcolor)(borderval)
    Wrs = fill(fill_val, gimgsz...)
    for i in 1:ncells
        gi = (i-1)÷gridsz[2]+1
        gj = i-(gi-1)*gridsz[2]
        gindices = (gi-1,gj-1, (add_dim_sz.-add_dim_sz)...)
        offset = gindices.*bimgsz .+ bordersz
        rngs = ntuple(i->offset[i]+1:offset[i]+imgsz[i], length(imgsz))
        Wrs[rngs...] = reshape(Wcolor[:,i], imgsz...)
    end
    Wrs
end

function mkimgH(H::Matrix{T}, tlength=size(H,2); colors=(colorant"green1", colorant"white", colorant"magenta")) where T
    mxabs = maximum(abs, H)
    fsc = scalesigned(mxabs)
    fcol = colorsigned(colors...)
    Array(mappedarray(fcol ∘ fsc, reshape(H[:,1:tlength], Val(2))))
end

"""
imshowW(W,imgsz; title="", gridcols=size(W,2), borderwidth=1, borderval=0.7, scalemtd=:maxwhole,
        colors=(colorant"green1", colorant"white", colorant"magenta"))
"""
function imshowW(W,imgsz; title="", kwargs...)
    if true #is_ImageView_available
        wimg = mkimgW(W,imgsz; kwargs...)
        gui_dict = ImageView.imshow(wimg)
        #set_gtk_property!(gui_dict["gui"]["window"], :title, title)
        gui_dict
    else
        @warn("ImageView is not available!")
        Dict()
    end
end

function imsaveW(fname,W,imgsz; kwargs...)
    wimg = mkimgW(W,imgsz; kwargs...)
    Images.save(fname, wimg)
    nothing
end

function imshowH(H,tlength=size(H,2); title="", kwargs...)
    if is_ImageView_available
        himg = mkimgH(H,tlength; kwargs...)
        gui_dict = ImageView.imshow(himg)
        set_gtk_property!(gui_dict["gui"]["window"], :title, title)
        gui_dict
    else
        @warn("ImageView is not available!")
        Dict()
    end
end

function imsaveH(fname,H,tlength=size(H,2); kwargs...)
    himg = mkimgH(H,tlength; kwargs...)
    Images.save(fname, himg)
    nothing
end


"""
    matchlist, ssds = matchcomponents(GT, W, errorfn::Function)
Matched the W columns with with those of GT.
GT: ground truch matrix
W: matrix to Compare
errorfn: function used to calculate the error of two vectors
matchlist: list of pair = (column index of GT, column index of W)
ssds: list of ssd for the matched pair columns
"""
function matchWcomponents(GT, W, errorfn::Function) # M X r form
    pq = PriorityQueue{Tuple{Int,Int,Bool}, Float64}(Base.Order.Forward) # Forward(low->high)
    gtcolnum = size(GT,2); wcolnum = size(W,2)
    for i = 1:gtcolnum
        gti = GT[:,i]
        for j = 1:wcolnum
            wj = W[:,j]
            dist, invert = errorfn(gti,wj)
            enqueue!(pq,(i,j,invert),dist)
        end
    end
    matchlist = Tuple{Int,Int,Bool}[]
    errs = Float64[]
    while !isempty(pq)
        p = peek(pq)
        dequeue!(pq)
        found = false
        mllength = length(matchlist)
        for i = 1:mllength
            if p[1][1] == matchlist[i][1] || p[1][2] == matchlist[i][2]
                found = true
                break
            end
        end
        if !found
            push!(matchlist,(p[1][1],p[1][2],p[1][3]))
            push!(errs,p[2])
        end
    end
    matchlist, errs
end

function matchcomponents(GTW::AbstractArray{T}, GTH::AbstractArray{T}, W::AbstractArray{T}, H::AbstractArray{T}; clamp=false) where T
    pq = PriorityQueue{Tuple{Int,Int,Bool}, T}(Base.Order.Forward) # Forward(low->high)
    gtcolnum = size(GTW,2); wcolnum = size(W,2)
    for i = 1:gtcolnum
        gtwi = GTW[:,i]; gthi = GTH[:,i]; gtxi = gtwi*gthi'
        for j = 1:wcolnum
            wj = W[:,j]; hj = H[j,:]; xj = wj*hj'
            clamp && (xj[xj.<0].=0)
            mnssd, invert = nssd(gtxi,xj)
            enqueue!(pq,(i,j,invert),mnssd)
        end
    end
    matchlist = Tuple{Int,Int,Bool}[]; ml = Int[]
    mnssds = T[]
    while !isempty(pq)
        p = peek(pq)
        dequeue!(pq)
        found = false
        mllength = length(matchlist)
        for i = 1:mllength
            if p[1][1] == matchlist[i][1] || p[1][2] == matchlist[i][2]
                found = true
                break
            end
        end
        if !found
            push!(matchlist,(p[1][1],p[1][2],p[1][3]))
            push!(ml,p[1][2])
            push!(mnssds,p[2])
        end
    end
    # Calculate unmatched power
    unmatchlist = collect(1:wcolnum)
    filter!(a->a ∉ ml,unmatchlist)
    gtxi = zeros(T,size(W,1),size(H,2))
    rerrs = T[]
    for j in unmatchlist
        wj = W[:,j]; hj = H[j,:]; xj = wj*hj'
        rerr = sum(abs2,xj)
        push!(rerrs,rerr)
    end
    matchlist, mnssds, rerrs
end

function fitcomponents(GTW::AbstractArray{T}, GTH::AbstractArray{T}, W::AbstractArray{T}, H::AbstractArray{T};
            clamp=false, iscalunmatched=false) where T
    pq = PriorityQueue{Tuple{Int,Int,Bool}, T}(Base.Order.Reverse) # Reverse(high->low)
    gtcolnum = size(GTW,2); wcolnum = size(W,2)
    for i = 1:gtcolnum
        gtwi = GTW[:,i]; gthi = GTH[:,i]; gtxi = gtwi*gthi'; ngtxi = norm(gtxi)
        for j = 1:wcolnum
            wj = W[:,j]; hj = H[j,:]; xj = wj*hj'
            clamp && (xj[xj.<0].=0)
            fitval = fitd(gtxi,xj,ngtxi)
            enqueue!(pq,(i,j,false),fitval)
        end
    end
    matchlist = Tuple{Int,Int,Bool}[]; ml = Int[]
    fitvals = T[]
    while !isempty(pq)
        p = peek(pq)
        dequeue!(pq)
        found = false
        mllength = length(matchlist)
        for i = 1:mllength
            if p[1][1] == matchlist[i][1] || p[1][2] == matchlist[i][2]
                found = true
                break
            end
        end
        if !found
            push!(matchlist,(p[1][1],p[1][2],p[1][3]))
            push!(ml,p[1][2])
            push!(fitvals,p[2])
        end
    end
    # Calculate unmatched power
    rerrs = T[]
    if iscalunmatched
        unmatchlist = collect(1:wcolnum)
        filter!(a->a ∉ ml,unmatchlist)
        gtxi = zeros(T,size(W,1),size(H,2))
        for j in unmatchlist
            wj = W[:,j]; hj = H[j,:]; xj = wj*hj'
            rerr = sum(abs2,xj)
            push!(rerrs,rerr)
        end
    end
    matchlist, fitvals, rerrs
end

function fitcomponents(X::AbstractArray, GTX::AbstractVector, W::AbstractArray{T}, H::AbstractArray{T};
            clamp=false) where T
    pq = PriorityQueue{Tuple{Int,Int,Bool}, T}(Base.Order.Reverse) # Reverse(high->low)
    gtcolnum = length(GTX); wcolnum = size(W,2)
    for i = 1:gtcolnum
        # read mask from GTX[i][2] and do the masking to X. vecs is a vector of vectors
        vecs = collect.(GTX[i][2]); idxs = eltype(vecs[1])[]; foreach(v->append!(idxs,v),vecs)
        gtxi = X[idxs,:]; foreach(c->c.=c.*GTX[i][3],eachcol(gtxi)); ngtxi = norm(gtxi)
        # calcuate fit value  with each W[:,j]*H[j,:]'
        for j = 1:wcolnum
            wj = W[idxs,j]; hj = H[j,:]; xj = wj*hj'
            clamp && (xj[xj.<0].=0)
            fitval = fitd(gtxi,xj,ngtxi)
            enqueue!(pq,(i,j,false),fitval)
        end
    end
    # find best matched pair (i,j)
    matchlist = Tuple{Int,Int,Bool}[]; ml = Int[]
    fitvals = T[]
    while !isempty(pq)
        p = peek(pq)
        dequeue!(pq)
        found = false
        mllength = length(matchlist)
        for i = 1:mllength
            if p[1][1] == matchlist[i][1] || p[1][2] == matchlist[i][2]
                found = true
                break
            end
        end
        if !found
            push!(matchlist,(p[1][1],p[1][2],p[1][3]))
            push!(ml,p[1][2])
            push!(fitvals,p[2])
        end
    end
    rerrs = T[]
    matchlist, fitvals, rerrs
end

fitx(a,b) = (m=sum(a)/length(a); denom=sum(abs2,a.-m); fitx(a,b,denom))
fitx(a,b,denom) = (1-sum(abs2,a-b)/denom)
fitd(a,b) = (na=norm(a); nb=norm(b); fitd(a,b,na,nb))
fitd(a,b,na) = (nb=norm(b); fitd(a,b,na,nb))
fitd(a,b,na,nb) = (denom=na^2+nb^2+2na*nb; 1-sum(abs2,a-b)/denom)
calfit(a,b) = (dval=fitd(a,b); aval=fitd(a,-b); dval > aval ? (dval, false) : (aval, true))
calfit(a,b,na) = (dval=fitd(a,b,na); aval=fitd(a,-b,na); dval > aval ? (dval, false) : (aval, true))
ssd(a,b) = sum(abs2,a-b)
nssd(a,b) = (ssd(a,b)/(norm(a)*norm(b)), false)
nssda(a,b) = (ssdval=ssd(a,b); ssaval=ssd(a,-b); nab=(norm(a)*norm(b));
                ssdval < ssaval ? (ssdval/nab, false) : (ssaval/nab, true))
function fiterr(a,b)
    init_x = eltype(a)[1, 0]
    f(x) = norm((x[1].*b.+x[2]).-a)^2
    rst = optimize(f,init_x)
    rst.minimum/length(a), false
end

matchedWnssd(GT,W) = ((ml, nssds) = matchWcomponents(GT, W, nssd); (sum(nssds)/length(nssds), ml, nssds))
matchedWnssda(GT,W) = ((ml, nssdas) = matchWcomponents(GT, W, nssda); (sum(nssdas)/length(nssdas), ml, nssdas))
matchedfitval(GTW, GTH, W, H; clamp=false) = ((ml, fitvals, rerrs) = fitcomponents(GTW, GTH, W, H; clamp=clamp);
                    (sum(fitvals)/size(GTW,2), ml, fitvals, rerrs))
matchedfitval(X, GTX::AbstractVector, W, H; clamp=false) = ((ml, fitvals, rerrs) = fitcomponents(X, GTX, W, H; clamp=clamp);
                    (sum(fitvals)/length(GTX), ml, fitvals, rerrs))
matchednssd(GTW, GTH, W, H; clamp=false) = ((ml, mnssds, rerrs) = matchcomponents(GTW, GTH, W, H; clamp=clamp);
                    (sum(mnssds)/length(mnssds), ml, mnssds, rerrs))
function matchedimg(W, matchlist)
    Wmimg = zeros(size(W,1),length(matchlist))
    for mp in matchlist
        Wmimg[:,mp[1]] = W[:,mp[2]]
    end
    Wmimg
end

function ssdH(ml,gtH,H)
    ssd = 0.
    for (gti, i, invert) in ml
        if i>size(H,2) # no match found
            ssd += sum((gtH[:,gti]).^2)
        else
            ssd += invert ? sum((gtH[:,gti]+H[:,i]).^2) : sum((gtH[:,gti]-H[:,i]).^2)
        end
    end
    ssd
end

# match order with gtW, then W[:,neworder] is same order with gtW
function matchedorder(ml,ncells)
    neworder = zeros(Int,length(ml))
    for (gti, i) in ml
        neworder[gti]=i
    end
    for i in 1:ncells
        i ∉ neworder && push!(neworder,i)
    end
    neworder
end

end