using TestData, Test, Colors, Images

# ── Helpers ──────────────────────────────────────────────────────────────────

has_display = TestData.is_ImageView_available

# ── Color utilities ───────────────────────────────────────────────────────────

@testset "color utilities" begin
    @test length(TestData.dgwm())  == 3
    @test length(TestData.bwm())   == 3
    @test length(TestData.g1wm())  == 3
    @test length(TestData.g1bw())  == 3
    @test length(TestData.bbw())   == 3
    @test length(TestData.bgw())   == 3
    @test length(TestData.wwb())   == 3
    @test all(c -> c isa Colorant, TestData.g1wm())
end

# ── interpolate ───────────────────────────────────────────────────────────────

@testset "interpolate" begin
    img = [1.0 2.0; 3.0 4.0]
    out = TestData.interpolate(img, 1)   # n=1 → each pixel becomes 2×2 tile, with 1-pixel border
    @test size(out) == size(img) .* 2
    @test out[1, 1] == img[1, 1]
end

# ── mkimgW / mkimgH ──────────────────────────────────────────────────────────

@testset "mkimgW" begin
    W = rand(16, 4)
    imgsz = (4, 4)
    img = TestData.mkimgW(W, imgsz)
    @test ndims(img) == 2
    @test eltype(img) <: Colorant

    for sm in [:maxwhole, :maxcol, :maxgridrow, :avgwhole, :avgcol]
        img2 = TestData.mkimgW(W, imgsz; scalemtd=sm)
        @test eltype(img2) <: Colorant
    end
end

@testset "mkimgH" begin
    H = rand(3, 20)
    img = TestData.mkimgH(H)
    @test size(img, 2) == 20
    @test eltype(img) <: Colorant
end

# ── imsaveW / imsaveH (file I/O, no display needed) ──────────────────────────

@testset "imsaveW / imsaveH" begin
    mktempdir() do d
        W = rand(16, 6)
        imgsz = (4, 4)
        fW = joinpath(d, "W.png")
        TestData.imsaveW(fW, W, imgsz)
        @test isfile(fW)
        img = load(fW)
        @test size(img, 1) > 0

        H = rand(6, 30)
        fH = joinpath(d, "H.png")
        TestData.imsaveH(fH, H)
        @test isfile(fH)
    end
end

# ── distance / fit metrics ────────────────────────────────────────────────────

@testset "ssd / nssd / fitd / fitx" begin
    a = [1.0, 2.0, 3.0]
    b = [1.0, 2.0, 3.0]
    @test TestData.ssd(a, b) ≈ 0.0

    c = [2.0, 3.0, 4.0]
    @test TestData.ssd(a, c) ≈ 3.0

    val, inv = TestData.nssd(a, b)
    @test val ≈ 0.0
    @test inv == false

    @test TestData.fitd(a, b) ≈ 1.0
    @test TestData.fitd(a, -a) ≈ 0.0

    @test TestData.fitx(a, b) ≈ 1.0
end

@testset "nssda" begin
    a = [1.0, 0.0, -1.0]
    val, inv = TestData.nssda(a, -a)
    @test inv == true       # -a is the better match
    @test val ≈ 0.0
end

# ── matchWcomponents / matchedWnssda ─────────────────────────────────────────

@testset "matchWcomponents" begin
    GT = [1.0 0.0; 0.0 1.0; 0.0 0.0]
    W  = [0.0 1.0; 1.0 0.0; 0.0 0.0]
    ml, errs = TestData.matchWcomponents(GT, W, TestData.nssd)
    @test length(ml) == 2
    @test length(errs) == 2
end

@testset "matchedWnssda" begin
    GT = [1.0 0.0; 0.0 1.0]
    W  = [0.9 0.1; 0.1 0.9]
    mssd, ml, nssds = TestData.matchedWnssda(GT, W)
    @test mssd isa Float64
    @test length(ml) == 2
end

# ── ssdH / matchedorder / matchedimg ─────────────────────────────────────────

@testset "ssdH / matchedorder / matchedimg" begin
    ml = [(1, 2, false), (2, 1, false)]
    gtH = [1.0 0.5; 0.5 1.0]
    H   = [0.9 0.4; 0.4 0.9]
    val = TestData.ssdH(ml, gtH, H)
    @test val isa Float64 && val >= 0.0

    neworder = TestData.matchedorder(ml, 2)
    @test length(neworder) == 2

    W = rand(4, 3)
    ml2 = [(1, 2, false), (2, 3, false)]
    Wm = TestData.matchedimg(W, ml2)
    @test size(Wm) == (4, 2)
end

# ── flip2makepos! ─────────────────────────────────────────────────────────────

@testset "flip2makepos!" begin
    W = [-1.0 2.0; -1.0 2.0; -1.0 2.0]
    H = [1.0 1.0 1.0; 1.0 1.0 1.0]
    TestData.flip2makepos!(W, H)
    @test all(sum(w for w in W[:, i] if w > 0) >= -sum(w for w in W[:, i] if w < 0) for i in 1:2)
end

# ── gaussian2D (synthetic data generation) ────────────────────────────────────

@testset "gaussian2D" begin
    ncells, imgrs, img2, gtW, gtH, gtWimgc, gtbg = TestData.gaussian2D(3.0, (20, 20), 100)
    @test ncells > 0
    @test size(imgrs, 2) == 100
    @test size(gtW, 1) == 20 * 20
    @test size(gtH, 2) == 100
end

@testset "gaussian_two_objs" begin
    n, imgrs, _, gtW, gtH, _, _ = TestData.gaussian_two_objs(3.0, (20, 20), 100, 8)
    @test n == 2
    @test size(imgrs) == (400, 100)
end

# ── loadfakecell (generate + save + reload) ───────────────────────────────────

@testset "loadfakecell generate" begin
    mktempdir() do d
        fname = joinpath(d, "test_fakecells.jld2")
        X, imgsz, dic, img_nl, maxSNR_X = TestData.loadfakecell(
            Float64, fname;
            sigma=3.0, lengthT=50, imgsz=(20, 20), SNR=5,
            issave=true, isload=false
        )
        @test size(X, 2) == 50
        @test imgsz == (20, 20)
        @test haskey(dic, "gt_ncells")
        @test isfile(fname)

        X2, imgsz2, dic2, _, _ = TestData.loadfakecell(
            Float64, fname;
            sigma=3.0, lengthT=50, imgsz=(20, 20), SNR=5,
            issave=false, isload=true
        )
        @test size(X2) == size(X)
    end
end

# ── Makie-based save functions (CairoMakie fallback — no window needed) ───────

@testset "heatmapWH" begin
    mktempdir() do d
        fprefix = joinpath(d, "test")
        W = rand(10, 3)
        H = rand(3, 20)
        TestData.heatmapWH(fprefix, W, H; Wdivision=1, Hdivision=1)
        @test isfile(fprefix * "_heatmap_W1.png")
        @test isfile(fprefix * "_heatmap_H1.png")
    end
end

@testset "plotH_data" begin
    mktempdir() do d
        fprefix = joinpath(d, "test")
        H = rand(3, 50)
        fig = TestData.plotH_data(fprefix, H)
        @test isfile(fprefix * "_plot_H.png")
    end
end

@testset "imsave_fakecell" begin
    mktempdir() do d
        fprefix = joinpath(d, "fakecell")
        W = rand(400, 4)
        H = rand(4, 50)
        imgsz = (20, 20)
        TestData.imsave_fakecell(fprefix, W, H, imgsz, 50)
        @test isfile(fprefix * "_W.png")
        @test isfile(fprefix * "_H.png")
    end
end

@testset "imsave_neurofinder" begin
    mktempdir() do d
        fprefix = joinpath(d, "neuro")
        W = rand(400, 5)
        H = rand(5, 50)
        TestData.imsave_neurofinder(fprefix, W, H, (20, 20), 50)
        @test isfile(fprefix * "_W.png")
        @test isfile(fprefix * "_H.png")
    end
end

# ── Display / window tests (skipped when no X11 / ImageView) ─────────────────

@testset "imshowW (requires display)" begin
    if !has_display
        @test_skip "No display available — skipping imshowW"
    else
        W = rand(16, 4)
        d = TestData.imshowW(W, (4, 4))
        @test d isa Dict
    end
end

@testset "imshowH (requires display)" begin
    if !has_display
        @test_skip "No display available — skipping imshowH"
    else
        H = rand(4, 30)
        d = TestData.imshowH(H)
        @test d isa Dict
    end
end
