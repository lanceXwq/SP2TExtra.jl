using SP2TExtra
using Test

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
        redirect_stdout(devnull) do
            @test toframes(readbin(joinpath(dir, "terminated.bin"), width=w, height=h)) == frames
            @test toframes(@test_logs warning readbin(joinpath(dir, "unterminated.bin"), width=w, height=h)) == frames
        end
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
        redirect_stdout(devnull) do
            read(; kwargs...) = readbin(dir; format=:dense, width=w, height=h, kwargs...)
            @test read(batchsize=1) == summed(bits, 1)
            @test read(batchsize=5) == summed(bits, 5)
            @test (@test_logs (:warn, r"only contains 2 frames") read(batchsize=4)) == summed(bits, 4)
            @test read(batchsize=5, crop_origin=(3, 2), crop_width=4, crop_height=2) ==
                  summed(bits[3:6, 2:3, :], 5)
            @test_throws ArgumentError read(batchsize=5, crop_origin=(6, 1), crop_width=4)
            @test_logs (:warn, r"Expected 6 frames") SP2TExtra.readbin_dense_directory(
                dir; width=w, height=h, batchsize=5, frames_per_file=6)
            @test_throws ArgumentError readbin(dir; width=w, height=h, crop_width=4)
        end
    end
end

# @testset "MSD tests" begin
#     x = randn(10, 5, 3)
#     @test isapprox(msd(x), sum(diff(x, dims=1) .^ 2) / (3 * (10 - 1)))
# end