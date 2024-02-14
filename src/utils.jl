"""
X = noisefilter(filter,X)
filter : :medT, :meanT, :medS
"""
function noisefilter(filter,X)
    if filter == :medT
        X = mapwindow(median!, X, (1,3)) # just for each row
    elseif filter == :meanT
        X = mapwindow(mean, X, (1,3)) # just for each row
    elseif filter == :medS
        rsimg = reshape(X,imgsz...,lengthT)
        rsimgm = mapwindow(median!, rsimg, (3,3,1))
        X = reshape(rsimgm,*(imgsz...),lengthT)
    end
    X
end

function flip2makepos!(W,H)
    p = size(W,2)
    for i in 1:p
        (w,h) = view(W,:,i), view(H,i,:)
        psum = sum(w[w.>0]); nsum = -sum(w[w.<0])
        psum < nsum && (w .*= -1; h .*= -1) # just '*=' doesn't work
    end
end

function flip2makepos!(W,H,Mw,Mh)
    p = size(W,2)
    for i in 1:p
        (w,h) = view(W,:,i), view(H,i,:); (mw,mh) = view(Mw,:,i), view(Mh,i,:)
        psum = sum(w[w.>0]); nsum = -sum(w[w.<0])
        psum < nsum && (w .*= -1; mw .*=-1; h .*= -1; mh .*= -1) # just '*=' doesn't work
    end
end
