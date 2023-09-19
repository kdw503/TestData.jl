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
