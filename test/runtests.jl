using SP2TExtra
using Test

# The result of f() and what it printed to stdout.
function captured(f)
    mktemp() do path, io
        result = redirect_stdout(f, io)
        flush(io)
        return result, read(path, String)
    end
end

@testset "XML tests" begin
    gttracks = [0.15367499430949305 0.22010250918040386 0.0;
        0.3161351020914415 0.27006612153286086 0.0;
        0.2711288587786091 0.27340129071031466 0.0;
        0.2703490014548096 0.2945822935280891 0.0;
        0.30439380732843224 0.3992083130600652 0.0;
        0.30925668022438907 0.4331489344514949 0.0;
        0.23547467805007016 0.45027805553041783 0.0;
        0.2625778731233346 0.5101664492168924 0.0;
        0.18825167441294663 0.4513309081261772 0.0;
        0.2716093990068874 0.529427986757348 0.0;
        0.3154063193003772 0.614647224629777 0.0;
        0.33022318733972866 0.4094173064084635 0.0;
        0.30588299411579484 0.4438255315805149 0.0;
        0.26892598888856545 0.4464031867929209 0.0;
        0.35138360313662426 0.4423236660738443 0.0;
        0.3863323462070774 0.34835328919284425 0.0;
        0.407679507824874 0.3068153660528501 0.0;
        0.37691638407786715 0.2587751980232748 0.0;
        0.3079657004069656 0.23887797333907246 0.0;
        0.44513955403705946 0.29784996879972314 0.0;
        0.354213834410956 0.24915317053513877 0.0;
        0.38883342690732947 0.25525323046170495 0.0;
        0.3835967456696548 0.29643673264425047 0.0;
        0.4606706600583911 0.3146070755425922 0.0;
        0.4992902758199351 0.5318347298915443 0.0;
        0.5718895296147325 0.3986045417262149 0.0;
        0.6457022521443803 0.4281646139001585 0.0;
        0.685187887834596 0.28938489274827994 0.0;
        0.7060370721393979 0.17069668011566905 0.0;
        0.669630541491682 0.21899871222352615 0.0;
        0.7032447453147959 0.27189352036659087 0.0;
        0.7360280328825375 0.2501445757922171 0.0;
        0.7251632119133569 0.23299726272920693 0.0;
        0.7047287916318168 0.2099248946689253 0.0]
    @test isapprox(xml2tracks("./testtracks.xml", nframes=34, batchsize=1), gttracks)

    # A track with a single detection, and a file with no tracks. Detection t (0-based) fills
    # frames t·batchsize .+ (1:batchsize); frames without a detection are NaN.
    trackmate(particles...) = """
        <?xml version="1.0" encoding="UTF-8"?>
        <Tracks nTracks="$(length(particles))" spaceUnits="micron" frameInterval="1.0" timeUnits="frame">
        $(join(particles, "\n"))
        </Tracks>
        """
    particle(detections...) = """
        <particle nSpots="$(length(detections))">
        $(join(("<detection t=\"$t\" x=\"$x\" y=\"$y\" z=\"$z\" />" for (t, x, y, z) in detections), "\n"))
        </particle>"""
    mktempdir() do dir
        path = joinpath(dir, "tracks.xml")

        write(path, trackmate(particle((0, 1.0, 2.0, 0.0), (1, 1.5, 2.5, 0.0)), particle((2, 3.0, 4.0, 0.5))))
        expected = fill(NaN, 3, 3, 2)
        expected[1:2, :, 1] = [1.0 2.0 0.0; 1.5 2.5 0.0]
        expected[3, :, 2] = [3.0, 4.0, 0.5]
        @test isequal(xml2tracks(path; nframes=3), expected)

        expected = fill(NaN, 6, 3, 2)
        expected[1:4, :, 1] = [1.0 2.0 0.0; 1.0 2.0 0.0; 1.5 2.5 0.0; 1.5 2.5 0.0]
        expected[5:6, :, 2] = [3.0 4.0 0.5; 3.0 4.0 0.5]
        @test isequal(xml2tracks(path; nframes=6, batchsize=2), expected)

        write(path, trackmate(particle((1, 3.0, 4.0, 0.5))))
        expected = fill(NaN, 3, 3, 1)
        expected[2, :, 1] = [3.0, 4.0, 0.5]
        @test isequal(xml2tracks(path; nframes=3), expected)

        write(path, trackmate())
        @test size(xml2tracks(path; nframes=3)) == (3, 3, 0)
    end
end

# The binary readers are documented (a string inside a function body documents nothing).
@testset "Docstring tests" begin
    for name in (:readbin, :getframes, :readbin_dense_directory, :crop_frames,
        :write_frames_to_tiff, :_readbin, :countsignals, :readsignals!, :signal2indices!,
        :popdelimiters!, :shiftindices!, :_getframes)
        @test Base.Docs.hasdoc(SP2TExtra, name)
    end
end

# A click at column i, row j of frame k has the linear index (k-1)·w·h + (j-1)·w + i.
# Extracting a region of clicks must equal slicing the frames, also when the region is not
# square. Every region below has a click in the last frame.
@testset "ROI tests" begin
    w, h = 5, 4
    frames = zeros(UInt16, w, h, 3)
    for (i, j, k) in ((2, 1, 1), (4, 3, 1), (3, 2, 2), (5, 4, 2), (2, 2, 3), (4, 3, 3))
        frames[i, j, k] = 1
    end
    indices = findall(!iszero, vec(frames))
    @test getframes(indices, width=w, height=h, batchsize=1) == frames
    for (origin, range) in (((2, 1, 1), (4, 3, 3)), ((2, 2, 1), (3, 3, 3)), ((2, 1, 1), (3, 2, 3)))
        expected = frames[origin[1]:origin[1]+range[1]-1, origin[2]:origin[2]+range[2]-1, :]
        ROIindices = extractROI(indices, (w, h), origin, range)
        @test getframes(ROIindices, width=range[1], height=range[2], batchsize=1) == expected
    end
    # the single-click and matrix methods agree with the definition
    cartesian = [2 1 1; 4 3 1; 3 2 2; 5 4 2; 2 2 3; 4 3 3]
    @test SP2TExtra.linearize(cartesian, w, h) == indices
    @test [SP2TExtra.linearize(c, w, h) for c in eachrow(cartesian)] == indices
    @test SP2TExtra.cartesianize(indices, w, w * h) == cartesian
end

# The sparse format lists each frame's clicks (0-based pixel indices, Int32) followed by a
# delimiter (the number of pixels). A recording that stops mid-frame has no delimiter after
# its last clicks: they are kept as a last frame, with a warning.
@testset "Sparse binary tests" begin
    w, h = 5, 4
    frames = zeros(UInt16, w, h, 3)
    for (i, j, k) in ((2, 1, 1), (4, 3, 1), (3, 2, 2), (5, 4, 2), (2, 2, 3), (4, 3, 3))
        frames[i, j, k] = 1
    end
    signals = Int32[]
    for k in 1:3
        append!(signals, findall(!iszero, vec(frames[:, :, k])) .- 1)
        push!(signals, w * h)
    end
    toframes(indices) = getframes(indices, width=w, height=h, batchsize=1)
    toindices(signals) = SP2TExtra.signal2indices!(signals .+ Int32(1), w * h)

    @test toframes(@test_logs toindices(signals)) == frames
    warning = (:warn, r"not a frame delimiter")
    @test toframes(@test_logs warning toindices(signals[1:end-1])) == frames
    # no delimiter at all: one frame
    @test toframes(@test_logs warning toindices(signals[1:2])) == frames[:, :, 1:1]
    @test isempty(@test_logs toindices(Int32[]))

    mktempdir() do dir
        write(joinpath(dir, "terminated.bin"), signals)
        write(joinpath(dir, "unterminated.bin"), signals[1:end-1])
        # readbin prints nothing; it logs one line
        info = (:info, r"Reading 1 \.bin file\(s\) from .*terminated\.bin")
        indices, output = captured() do
            @test_logs info readbin(joinpath(dir, "terminated.bin"), width=w, height=h)
        end
        @test toframes(indices) == frames
        @test isempty(output)
        @test toframes(@test_logs info warning readbin(joinpath(dir, "unterminated.bin"), width=w, height=h)) == frames

        # a recording split across files, in name order, reads as one stream, also when a
        # frame starts in one file and ends in the next (signals[4] is frame 2's first click)
        split = mkdir(joinpath(dir, "split"))
        write(joinpath(split, "a.bin"), signals[1:4])
        write(joinpath(split, "b.bin"), signals[5:end])
        info = (:info, r"Reading 2 \.bin file\(s\)")
        @test toframes(@test_logs info readbin(split, width=w, height=h)) == frames
    end
end

# getframes with batches: a batch sums batchsize raw frames, and a last partial batch is kept
# (with a warning). Without nframes the frames end at the last click; with it, trailing
# frames without clicks are kept.
@testset "Batching tests" begin
    w, h = 3, 2
    frames = zeros(UInt16, w, h, 5)
    for (i, j, k) in ((1, 1, 1), (3, 2, 1), (2, 1, 2), (1, 1, 3), (3, 1, 3))
        frames[i, j, k] = 1
    end # frames 4 and 5 have no clicks
    indices = findall(!iszero, vec(frames))
    batched(b) = cat((sum(frames[:, :, k:min(k + b - 1, 5)], dims=3) for k in 1:b:5)..., dims=3)

    @test getframes(indices, width=w, height=h, batchsize=1) == frames[:, :, 1:3]
    @test getframes(indices, width=w, height=h, batchsize=1, nframes=5) == frames
    @test getframes(indices, width=w, height=h, batchsize=5, nframes=5) == batched(5)
    partial = (:warn, r"last batch sums 1 frames")
    @test (@test_logs partial getframes(indices, width=w, height=h, batchsize=2, nframes=5)) == batched(2)
    @test (@test_logs partial getframes(indices, width=w, height=h, batchsize=2)) == batched(2)[:, :, 1:2]
    @test getframes(Int[], width=w, height=h, batchsize=1, nframes=2) == zeros(UInt16, w, h, 2)
    @test_throws ArgumentError getframes(Int[], width=w, height=h, batchsize=1)
    @test_throws ArgumentError getframes(indices, width=w, height=h, batchsize=1, nframes=2)

    # a region whose last frame has no clicks keeps all its frames with nframes
    ROIindices = extractROI(indices, (w, h), (1, 1, 1), (2, 2, 5))
    @test getframes(ROIindices, width=2, height=2, batchsize=1, nframes=5) == frames[1:2, :, :]
end

# The dense format stores one bit per pixel, most significant bit first, frame by frame in
# column-major order. Reading sums batches of frames (keeping a partial last batch, with a
# warning) and can crop; any frame size works.
@testset "Dense binary tests" begin
    w, h = 8, 4
    bits = BitArray(isodd(i * j + k) && (i + 2j + 3k) % 5 != 0 for i in 1:w, j in 1:h, k in 1:10)
    pack(b) = [reduce(|, UInt8(b[8(n-1)+m]) << (8 - m) for m in 1:8) for n in 1:length(b)÷8]
    summed(b, s) = cat((sum(UInt16.(b[:, :, k:min(k + s - 1, end)]), dims=3) for k in 1:s:size(b, 3))..., dims=3)
    mktempdir() do dir
        write(joinpath(dir, "1.bin"), pack(bits[:, :, 1:6]))
        write(joinpath(dir, "2.bin"), pack(bits[:, :, 7:10]))
        read(; kwargs...) = readbin(dir; format=:dense, width=w, height=h, kwargs...)
        info = (:info, r"Reading 2 \.bin file\(s\)")
        binned, output = captured(() -> @test_logs info read(batchsize=1))
        @test binned == summed(bits, 1)
        @test isempty(output)
        @test (@test_logs info read(batchsize=5)) == summed(bits, 5)
        @test (@test_logs info (:warn, r"only contains 2 frames") read(batchsize=4)) == summed(bits, 4)
        @test (@test_logs info read(batchsize=5, crop_origin=(3, 2), crop_width=4, crop_height=2)) ==
              summed(bits[3:6, 2:3, :], 5)
        @test_logs info @test_throws ArgumentError read(batchsize=5, crop_origin=(6, 1), crop_width=4)
        @test_logs (:warn, r"Expected 6 frames") SP2TExtra.readbin_dense_directory(
            dir; width=w, height=h, batchsize=5, frames_per_file=6)
        # an invalid call fails before logging anything
        @test_logs @test_throws ArgumentError readbin(dir; width=w, height=h, crop_width=4)
    end
end

# The MSD is per coordinate, as in SP2T (2·D·period): particle p steps ±s[p] along every
# coordinate in every frame, so its one-frame MSD is s[p]², whatever the number of
# coordinates. With steps of constant sign, the displacement over L frames is L·s[p], so the
# time-averaged MSD at lag L is L²·s[p]². On irregular tracks, tamsd matches a direct
# average over all pairs of frames L apart.
@testset "MSD tests" begin
    s = [0.5, 1.0, 2.0]
    for ndims in (2, 3)
        steps = [(-1)^(t + d) * s[p] for t in 1:9, d in 1:ndims, p in 1:3]
        x = cumsum(cat(zeros(1, ndims, 3), steps; dims=1); dims=1)
        @test vec(tamsd(x; lags=1)) ≈ s .^ 2
        @test msd(x) ≈ sum(s .^ 2) / 3

        steps = [(-1)^d * s[p] for t in 1:9, d in 1:ndims, p in 1:3]
        x = cumsum(cat(zeros(1, ndims, 3), steps; dims=1); dims=1)
        @test tamsd(x) ≈ [L^2 * s[p]^2 for L in 1:9, p in 1:3]
        @test tamsd(x; lags=[5, 2]) ≈ [L^2 * s[p]^2 for L in [5, 2], p in 1:3]
    end

    x = [sin(t * d + p) + t^2 / 7 for t in 1:12, d in 1:2, p in 1:3]
    direct(L, p) = sum((x[t+L, d, p] - x[t, d, p])^2 for t in 1:12-L, d in 1:2) / (2 * (12 - L))
    @test tamsd(x) ≈ [direct(L, p) for L in 1:11, p in 1:3]
    @test msd(x) ≈ sum(tamsd(x; lags=1)) / 3
    @test_throws ArgumentError tamsd(x; lags=0:2)
    @test_throws ArgumentError tamsd(x; lags=[12])
end

# The localization error is the RMS over frames and coordinates of the per-coordinate
# posterior standard deviation. Samples at μ ± s have the variance s² (over the samples
# themselves), so the error is sqrt(mean(s²)). A sample far away is skipped as burn-in;
# samples with other than one emitter are rejected.
@testset "Localization error tests" begin
    S = SP2TExtra.SP2T
    sample(x) = S.Sample(x, 1.0, 1.0, 0, 1.0, NaN, NaN)
    chain(xs...) = S.Chain([sample(x) for x in xs], 10, S.ConstantAnnealing{Float64}(1))
    μ = reshape([1.0 2.0; 1.5 2.5; 2.0 3.0], 3, 2, 1)
    s = reshape([0.1 0.2; 0.3 0.1; 0.2 0.4], 3, 2, 1)
    far = fill(50.0, 3, 2, 1)
    @test localization_error(chain(far, μ .+ s, μ .- s); burn_in=1) ≈ sqrt(sum(s .^ 2) / 6)
    @test localization_error(chain(μ .+ s, μ .- s, μ .+ s, μ .- s)) ≈ sqrt(sum(s .^ 2) / 6)
    @test_throws ArgumentError localization_error(chain(μ .+ s, cat(μ, μ .- s; dims=3)))
    @test_throws ArgumentError localization_error(chain(μ .+ s, μ .- s); burn_in=2)
end

# tracks(chain; burn_in) concatenates, in order, the emitters of every sample after the
# first burn_in, whatever their number.
@testset "Chain tracks tests" begin
    S = SP2TExtra.SP2T
    sample(x) = S.Sample(x, 1.0, 1.0, 0, 1.0, NaN, NaN)
    chain(xs...) = S.Chain([sample(x) for x in xs], 10, S.ConstantAnnealing{Float64}(1))
    a, b, c = (reshape(Float64.(n .+ (1:6)), 3, 2, 1) for n in (0, 10, 20))
    d = reshape(Float64.(31:36), 3, 2, 1)
    @test tracks(chain(a, cat(b, c; dims=3), d)) == cat(a, b, c, d; dims=3)
    @test tracks(chain(a, cat(b, c; dims=3), d); burn_in=1) == cat(b, c, d; dims=3)
    @test size(tracks(chain(a, b); burn_in=2)) == (3, 2, 0)
end

# credible1D and credible2D give, for each frame, the fraction of sampled positions in each
# bin: for one particle (i) or pooled over all. Two identical samples, two frames: particle 1
# is at (0.5, 1.5) then (1.5, 2.5), particle 2 at (2.5, 0.5) in both; bins 0:1:3 on each
# axis. credible2D's rows are x bins and its columns y bins.
@testset "Credible region tests" begin
    S = SP2TExtra.SP2T
    positions = cat([0.5 1.5; 1.5 2.5], [2.5 0.5; 2.5 0.5]; dims=3) # (frame, x/y, particle)
    samples = [S.Sample(positions, 1.0, 1.0, 0, 1.0, NaN, NaN) for _ in 1:2]
    edges = 0:1:3

    x, y = credible1D(samples, 1, edges, edges)
    @test x == [1 0 0; 0 1 0]
    @test y == [0 1 0; 0 0 1]
    x, y = credible1D(samples, edges, edges)
    @test x == [0.5 0 0.5; 0 0.5 0.5]
    @test y == [0.5 0.5 0; 0.5 0 0.5]

    onehot(i, j) = (h = zeros(3, 3); h[i, j] = 1; h)
    c = credible2D(samples, 1, edges, edges)
    @test c[:, :, 1] == onehot(1, 2) # x bin 1, y bin 2
    @test c[:, :, 2] == onehot(2, 3)
    c = credible2D(samples, edges, edges)
    @test c[:, :, 1] == (onehot(1, 2) + onehot(3, 1)) / 2
    @test c[:, :, 2] == (onehot(2, 3) + onehot(3, 1)) / 2

    # factor and shift: positions are shifted, then scaled, before binning
    x, ~ = credible1D(samples, 1, 0:2:6, edges; factor=2, xshift=0.5)
    @test x == [1 0 0; 0 1 0] # 2·(0.5 + 0.5) = 2 and 2·(1.5 + 0.5) = 4
end

# A TIFF written by writetiff reads back the same: frames, period, unit and pixel size. TIFF
# stores the resolution (pixels per unit, rounded by writetiff to 3 decimals), so the pixel
# size comes back exactly for 0.1 and within that rounding for 0.065.
@testset "TIFF tests" begin
    frames = UInt16[(i + 3j + 7k) % 9 for i in 1:6, j in 1:4, k in 1:3]
    mktempdir() do dir
        path = joinpath(dir, "frames.tif")
        writetiff(path, frames; px_size=0.1, unit="nm", period=0.01)
        readouts, metadata = readtiff(path)
        @test readouts == frames
        @test metadata["pixel size"] == 0.1
        @test metadata["period"] == 0.01
        @test metadata["unit"] == "nm"

        writetiff(path, frames; px_size=0.065)
        @test isapprox(readtiff(path)[2]["pixel size"], 0.065; rtol=1e-4)

        # write_frames_to_tiff converts to UInt16 and prints nothing
        ~, output = captured(() -> SP2TExtra.write_frames_to_tiff(path, Int.(frames); pixel_size=0.2))
        @test isempty(output)
        readouts, metadata = readtiff(path)
        @test readouts == frames
        @test metadata["pixel size"] == 0.2
    end
end

# read_stage_positions reads the Timestamp_ms, X_microns and Y_microns columns (in any order,
# among others) of a stage CSV. By default times are seconds and positions microns, both
# relative to the first row; blank lines are skipped.
@testset "Stage position tests" begin
    mktempdir() do dir
        path = joinpath(dir, "stage.csv")
        write(path, """
            Y_microns, Other, Timestamp_ms, X_microns
            20.0, a, 1000, 5.0
            21.5, b, 1250, 4.0

            19.0, c, 1500, 6.5
            """)
        t, x, y = SP2TExtra.read_stage_positions(path)
        @test t ≈ [0, 0.25, 0.5]
        @test x ≈ [0, -1, 1.5]
        @test y ≈ [0, 1.5, -1]
        t, x, y = SP2TExtra.read_stage_positions(path; relative_time=false,
            relative_position=false, time_unit=:ms)
        @test t == [1000, 1250, 1500]
        @test x == [5, 4, 6.5]
        @test y == [20, 21.5, 19]
        @test_throws ErrorException SP2TExtra.read_stage_positions(path; time_unit=:min)

        write(path, "Timestamp_ms,X_microns\n0,1\n")
        @test_throws ErrorException SP2TExtra.read_stage_positions(path) # no Y_microns
        write(path, "Timestamp_ms,X_microns,Y_microns\n0,1\n")
        @test_throws ErrorException SP2TExtra.read_stage_positions(path) # short row
    end
end

# binframes sums and bintracks averages each batch of frames, a partial last batch included,
# so the two stay aligned frame for frame. With x[t] = t + 10d + 100p, the binned position is
# the mean frame number of the batch. Only binframes warns about a partial batch.
@testset "Binning tests" begin
    frames = UInt8[(i * j + k) % 2 for i in 1:4, j in 1:3, k in 1:10]
    x = [t + 10.0d + 100p for t in 1:10, d in 1:2, p in 1:3]
    for (batchsize, batches) in ((5, [1:5, 6:10]), (3, [1:3, 4:6, 7:9, 10:10]))
        expected = cat((sum(frames[:, :, b]; dims=3) for b in batches)...; dims=3)
        binned = if length(batches[end]) < batchsize
            @test_logs (:warn, r"last batch sums 1 frames") binframes(frames, batchsize)
        else
            @test_logs binframes(frames, batchsize)
        end
        @test binned == expected
        means = [sum(b) / length(b) for b in batches]
        @test (@test_logs bintracks(x, batchsize)) ≈
              [m + 10d + 100p for m in means, d in 1:2, p in 1:3]
        @test size(bintracks(x, batchsize), 1) == size(binned, 3)
    end

    # binframes returns UInt16, SP2T's readout type, whatever the input type; batch sums
    # must fit in it: 65535 passes, 65536 and negative sums throw.
    for input in (trues(2, 2, 600), ones(UInt8, 2, 2, 600), ones(Int, 2, 2, 600))
        binned = binframes(input, 300)
        @test eltype(binned) == UInt16
        @test binned == fill(300, 2, 2, 2)
    end
    @test binframes(fill(UInt8(65), 1, 1, 1008), 1008) == fill(65520, 1, 1, 1)
    @test binframes([fill(UInt16(65535), 1, 1, 1);;; zeros(UInt16, 1, 1, 2)], 3) ==
          fill(65535, 1, 1, 1)
    @test_throws ArgumentError binframes([fill(UInt16(65535), 1, 1, 1);;; ones(UInt16, 1, 1, 2)], 3)
    @test_throws ArgumentError binframes(fill(-1, 1, 1, 2), 2)
end

# A pixel that clicks in k of N binary frames has the dark count -log(1 - k/N): 4 frames of
# batchsize 3 give N = 12, and k = 0, 3, 6, 12 give 0, log(4/3), log(2), Inf. A pixel that
# always clicks gives Inf with a warning, and SP2T filters it out. Readouts above the
# batchsize mean a wrong batchsize.
@testset "Dark count tests" begin
    darkframes = zeros(UInt8, 2, 2, 4)
    darkframes[2, 1, :] = [1, 0, 2, 0] # k = 3
    darkframes[1, 2, :] = [3, 0, 3, 0] # k = 6
    darkframes[2, 2, :] .= 3 # k = 12
    darkcounts = @test_logs (:warn, r"1 pixel\(s\) clicked in every binary frame") getdarkcounts(
        darkframes; batchsize=3)
    @test darkcounts ≈ [0 log(2); log(4 / 3) Inf]
    detector = SP2TExtra.SP2T.SPAD{Float64}(period=1, pixel_size=0.1, darkcounts=darkcounts,
        cutoffs=(0, Inf), readouts=zeros(UInt16, 2, 2, 1), batchsize=3)
    @test detector.filter == [1 1; 1 0]

    # binary frames: k = 0, 2, 3 of N = 4 give 0, log(2), log(4)
    binary = zeros(UInt8, 2, 2, 4)
    binary[2, 1, :] = [1, 0, 1, 0]
    binary[1, 2, :] = [1, 1, 1, 0]
    @test (@test_logs getdarkcounts(binary)) ≈ [0 log(4); log(2) 0]

    @test_throws ArgumentError getdarkcounts(darkframes)
end

@testset "Viewer tests" begin
    @test isnothing(Base.get_extension(SP2TExtra, :SP2TExtraMakieExt))
    @eval using CairoMakie
    CairoMakie.activate!(visible = false) # display renders without opening a viewer
    @test !isnothing(Base.get_extension(SP2TExtra, :SP2TExtraMakieExt))

    # SP2T's pixel (i, j) covers [(i - 1)·pxsize, i·pxsize] × [(j - 1)·pxsize, j·pxsize]
    width, height, pxsize = 12, 6, 0.1
    frames = zeros(UInt16, width, height, 3)
    frames[10, 2, :] .= 5
    # from the center of pixel (10, 2) to that of pixel (2, 6)
    tracks = reshape(([10 2; 6 4; 2 6] .- 1 / 2) .* pxsize, 3, 2, 1)
    @test isnothing(first(captured(() -> viewframes(frames; tracks, pxsize))))
    fig = current_figure()
    ax = only(c for c in fig.content if c isa Axis)
    plotsof(T) = [p for p in ax.scene.plots if p isa T]
    @test ax.yreversed[]
    heatmap = only(plotsof(Heatmap))
    xedges, yedges = heatmap[1][], heatmap[2][]
    @test xedges == 0:width && yedges == 0:height # pixel edges, in pixels
    @test heatmap[3][] == frames[:, :, 1] # neither rotated nor transposed
    frameslider = first(only(c for c in fig.content if c isa SliderGrid).sliders)
    set_close_to!(frameslider, 3)
    track = only(plotsof(Lines))[1][]
    @test length(track) == 3
    @test xedges[10] < track[1][1] < xedges[11] && yedges[2] < track[1][2] < yedges[3]

    regionsliders = [c for c in fig.content if c isa IntervalSlider]
    xslider = only(s for s in regionsliders if s.horizontal[])
    yslider = only(s for s in regionsliders if !s.horizontal[])
    region = only(plotsof(Poly))[1]
    @test region[] == Rect(0, 0, width, height)
    set_close_to!(xslider, 2, 8) # x edges 2 to 8
    set_close_to!(yslider, 1, 5) # counted from the bottom: y edges 6 - 5 = 1 to 6 - 1 = 5
    @test region[] == Rect(2, 1, 6, 4)
    @test sort([only(p.text[]) for p in plotsof(CairoMakie.Text)]) == ["(3, 2)", "(6, 4)"] # first pixel, size
end
