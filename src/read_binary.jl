# read_binary.jl
#
# Unified interface for reading binary image data in dense and sparse formats.
# - Dense: All pixels stored as bits
# - Sparse: Only activated pixel indices stored as 32-bit integers

# ============================================================================
# SHARED UTILITIES - Coordinate Conversion
# ============================================================================

function cartesianize(
    indices::AbstractVector{T},
    width::Integer,
    framesize::Integer,
) where {T<:Integer}
    cartesianindices = Matrix{T}(undef, length(indices), 3)
    return cartesianize!(cartesianindices, indices, width, framesize)
end

function cartesianize!(
    cartesianindices::AbstractMatrix{<:Integer},
    indices::AbstractVector{<:Integer},
    width::Integer,
    framesize::Integer,
)
    cartesianindices[:, 1] .= mod1.(indices, width)
    cartesianindices[:, 2] .= fld1.(mod1.(indices, framesize), width)
    cartesianindices[:, 3] .= fld1.(indices, framesize)
    return cartesianindices
end

function cartesianize!(
    cartesianindex::AbstractVector{<:Integer},
    index::Integer,
    width::Integer,
    framesize::Integer,
)
    cartesianindex[1] = mod1(index, width)
    cartesianindex[2] = fld1(mod1(index, framesize), width)
    cartesianindex[3] = fld1(index, framesize)
    return cartesianindex
end

linearize(cartesianindices::AbstractVector{<:Integer}, width::Integer, height::Integer) =
    (cartesianindices[3] - 1) * height * width +
    (cartesianindices[2] - 1) * width +
    cartesianindices[1]

function linearize(
    cartesianindices::AbstractMatrix{<:Integer},
    width::Integer,
    height::Integer,
)
    indices = similar(cartesianindices, axes(cartesianindices, 1))
    return linearize!(indices, cartesianindices, width, height)
end

function linearize!(
    indices::AbstractVector{<:Integer},
    cartesianindices::AbstractMatrix{<:Integer},
    width::Integer,
    height::Integer,
)
    @views @. indices =
        (cartesianindices[:, 3] - 1) * height * width +
        (cartesianindices[:, 2] - 1) * width +
        cartesianindices[:, 1]
    return indices
end

# ============================================================================
# SHARED UTILITIES - ROI Extraction
# ============================================================================

function extractROI!(cartesianindices, ROIbounds)
    newupperbounds = view(ROIbounds, 1:1, :) .- 1
    cartesianindices .-= newupperbounds
    newupperbounds .= view(ROIbounds, 2:2, :) .- newupperbounds
    inROI = vec(all(0 .< cartesianindices .≤ newupperbounds, dims=2))
    return cartesianindices[inROI, :]
end

function extractROI(
    indices::AbstractVector{<:Integer},
    framesize::NTuple{2,<:Integer},
    ROIorigin::NTuple{3,<:Integer},
    ROIrange::NTuple{3,<:Integer},
)
    cartesianindex = similar(indices, 3)
    newindices = similar(indices, 0)
    for index in indices
        cartesianize!(cartesianindex, index, framesize[1], framesize[1] * framesize[2])
        cartesianindex .-= ROIorigin .- 1
        if cartesianindex[3] > ROIrange[3]
            break
        elseif 0 < cartesianindex[1] <= ROIrange[1] &&
               0 < cartesianindex[2] <= ROIrange[2] &&
               0 < cartesianindex[3]
            push!(newindices, linearize(cartesianindex, ROIrange[1], ROIrange[2]))
        end
    end
    return newindices
end

# ============================================================================
# SHARED UTILITIES - File Path Handling
# ============================================================================

isbinfile(path) = isfile(path) && endswith(path, ".bin")

function getbinfiles(dir)
    isdir(dir) || return String[]
    filter(f -> endswith(f, ".bin"), readdir(dir, join=true, sort=true))
end

function checkpath(path)
    isbinfile(path) && return [path]
    files = getbinfiles(path)
    isempty(files) && throw(ArgumentError("No .bin files found in path: $path"))
    return files
end

# ============================================================================
# DENSE FORMAT - Bit Unpacking & Reading
# ============================================================================

const BYTE_TO_BITS_MSB = [
    ((UInt8(byte) >> UInt8(bitpos)) & 0x01) == 0x01
    for bitpos in 7:-1:0, byte in 0:255
]

function readbin_dense_one_file(file_path; width=512, height=512)
    reshape(BYTE_TO_BITS_MSB[:, read(file_path).+1], width, height, :)
end

"""
    crop_frames(frames, crop_width, crop_height, crop_origin)

Crop frames to a region of interest (ROI).

crop_origin is 1-based, e.g. (1, 1) means top-left pixel.
"""
function crop_frames(frames, crop_width, crop_height, crop_origin)
    x0 = crop_origin[1]
    y0 = crop_origin[2]
    x1 = x0 + crop_width - 1
    y1 = y0 + crop_height - 1
    return frames[x0:x1, y0:y1, :]
end

"""
    readbin_dense_directory(path; width=512, height=512, frames_per_file=nothing, batchsize=1000, crop_width=width, crop_height=height, crop_origin=(1, 1))

Read dense binary files from a directory with batching and optional cropping.

The crop defaults to the whole frame. Given `frames_per_file`, warn about files with a
different number of frames.
"""
function readbin_dense_directory(
    path;
    width=512,
    height=512,
    frames_per_file=nothing,
    batchsize=1000,
    crop_width=width,
    crop_height=height,
    crop_origin=(1, 1),
)
    files = checkpath(path)
    file_sizes_in_bits = filesize.(files) * 8
    n_frames_per_file = file_sizes_in_bits .÷ (width * height)
    n_bits_left_per_file = rem.(file_sizes_in_bits, width * height)
    any(n_bits_left_per_file .!= 0) && throw(ArgumentError("Some files have extra bits that do not form a complete $width × $height frame. File sizes in bits: $file_sizes_in_bits, bits left per file: $n_bits_left_per_file"))

    n_frames_batched = cld(sum(n_frames_per_file), batchsize)
    n_remaining_frames = rem(sum(n_frames_per_file), batchsize)

    n_remaining_frames != 0 && @warn "Last batched frame only contains $n_remaining_frames frames instead of $batchsize"

    any(crop_origin .< 1) && throw(ArgumentError("Crop_origin must be 1-based, e.g. (1, 1) for the top-left pixel"))
    any(crop_origin .+ (crop_width, crop_height) .- 1 .> (width, height)) && throw(ArgumentError("Selected ROI exceeds FOV bounds"))

    batched_frames = Array{UInt16}(undef, crop_width, crop_height, n_frames_batched)

    batch_idx = 1
    frames_in_batch = 0
    batch_accumulator = zeros(UInt16, crop_width, crop_height)

    @showprogress dt = 3 desc = "Reading binary image data" for file in files
        frames = readbin_dense_one_file(file; width=width, height=height)

        !isnothing(frames_per_file) && size(frames, 3) != frames_per_file &&
            @warn "Expected $frames_per_file frames in $file, found $(size(frames, 3))"

        if crop_width < width || crop_height < height
            frames = crop_frames(frames, crop_width, crop_height, crop_origin)
        end

        for frame in eachslice(frames, dims=3)
            batch_accumulator .+= frame
            frames_in_batch += 1

            if frames_in_batch == batchsize
                copyto!(view(batched_frames, :, :, batch_idx), batch_accumulator)
                batch_idx += 1
                frames_in_batch = 0
                fill!(batch_accumulator, 0)
            end
        end
    end

    n_remaining_frames > 0 && copyto!(view(batched_frames, :, :, batch_idx), batch_accumulator)

    return batched_frames
end

# ============================================================================
# SPARSE FORMAT - Index Reading & Frame Generation
# ============================================================================

"""
    _readbin(files)

Internal: Read raw binary data from files as 32-bit integers.
"""
function _readbin(files::AbstractVector)
    nsignals = countsignals(files)
    signals = Vector{Int32}(undef, sum(nsignals))
    readsignals!(files, signals, nsignals)
    return signals
end

"""
    countsignals(files)

Count signals (32-bit integers) in each file.
"""
function countsignals(files::AbstractVector)
    sizes = filesize.(files)
    corrupted = iscorrupted.(sizes)
    any(corrupted) && throw(
        DomainError(
            sizes[corrupted],
            "File size is not divisible by 4 in $(files[corrupted])",
        ),
    )
    return sizes .÷ 4
end

iscorrupted(sizeinbyte) = sizeinbyte % 4 != 0

"""
    readsignals!(files, signals, nsignals)

Read raw signals from files and increment by 1 (0-indexing to 1-indexing).
"""
function readsignals!(
    files::AbstractVector{String},
    signals::Vector{<:Integer},
    nsignals::AbstractVector{<:Integer},
)
    start = 1
    @inbounds for (file, number) in zip(files, nsignals)
        read!(file, view(signals, range(start, length=number)))
        start += number
    end
    return signals .+= 1
end

"""
    signal2indices!(indices, framesize)

Convert raw signals to pixel indices by removing frame delimiters.
"""
function signal2indices!(indices::AbstractVector, framesize)
    delimiters = popdelimiters!(indices, framesize)
    indices = convert(Vector{Int}, indices)
    return shiftindices!(indices, delimiters, framesize)
end

"""
    popdelimiters!(signals, framesize)

Remove delimiter values from signal array. A delimiter is framesize in the file, so
framesize + 1 after `readsignals!` adds 1: every signal above framesize.
"""
function popdelimiters!(signals::AbstractVector, framesize)
    delimiters = findall(>(framesize), signals)
    # a recording that stops mid-frame has clicks after its last delimiter
    unterminated = !isempty(signals) && (isempty(delimiters) || delimiters[end] != length(signals))
    deleteat!(signals, delimiters)
    delimiters .-= 1:length(delimiters)
    if unterminated
        @warn "The last signal is not a frame delimiter: the clicks after the last delimiter are kept as a last frame, which may be incomplete."
        push!(delimiters, length(signals)) # close the last frame
    end
    return delimiters
end

"""
    shiftindices!(indices, delimiters, framesize)

Add appropriate framesize offsets to indices based on frame boundaries.
"""
function shiftindices!(indices, delimiters, framesize)
    @inbounds for i = 1:length(delimiters)-1
        @views indices[delimiters[i]+1:delimiters[i+1]] .+= i * framesize
    end
    return indices
end

"""
    getframes(indices; width, height, batchsize, nframes=nothing, threshold=typemax(Int))

Convert sparse indices to batched frame arrays.

`nframes` is the number of raw frames. Without it, the frames end with the last click,
so trailing frames without clicks are lost (as is common after `extractROI`). A last
batch of fewer than `batchsize` raw frames is kept, with a warning: it sums fewer frames
than the others.
"""
function getframes(
    indices::AbstractVector{<:Integer};
    width::Integer,
    height::Integer,
    batchsize::Integer,
    nframes::Union{Integer,Nothing}=nothing,
    threshold::Integer=typemax(Int),
)
    framesize = width * height
    if isnothing(nframes)
        isempty(indices) &&
            throw(ArgumentError("There are no clicks to count the frames from; pass `nframes`."))
        nframes = fld1(maximum(indices), framesize)
    elseif !isempty(indices) && maximum(indices) > nframes * framesize
        throw(ArgumentError("There are clicks after frame $nframes (`nframes`)."))
    end
    nbatches, framesleft = countframes(nframes, batchsize)
    framesleft != 0 &&
        @warn "The last batch sums $framesleft frames, fewer than the batchsize ($batchsize)."
    frames = _getframes(indices, framesize, batchsize, nbatches)
    for i in eachindex(frames)
        frames[i] .> threshold && (frames[i] = 0)
    end
    return reshape(frames, width, height, nbatches)
end

# the number of batches, the last one partial when the remainder is not zero
countframes(nframes::Integer, batchsize::Integer) =
    (cld(nframes, batchsize), rem(nframes, batchsize))

"""
    _getframes(indices, framesize, batchsize, count)

Internal: Generate frame array from indices using histogram counting.
"""
function _getframes(
    indices::AbstractVector{<:Integer},
    framesize::Integer,
    batchsize::Integer,
    count::Integer,
)
    npixels = framesize * count
    frames = zeros(UInt16, npixels)
    addcounts!(frames, batchindices(indices, framesize, batchsize), 1:npixels)
    return frames
end

batchindices(indices, framesize, batchsize) =
    @. (indices - 1) ÷ (framesize * batchsize) * framesize + mod1(indices, framesize)

# ============================================================================
# UNIFIED API - Main Entry Points
# ============================================================================

"""
    readbin(path; format=:sparse, width=512, height=512, batchsize=1000, crop_origin, crop_width, crop_height)

Read binary image data in dense or sparse format.

# Arguments
- `path`: File or directory path
- `format`: `:dense` or `:sparse` (default)
- `width`, `height`: Image dimensions
- `batchsize`: Frames per batch (dense only)
- `crop_origin`, `crop_width`, `crop_height`: Region to keep while reading (dense only;
  default: the whole frame). `crop_origin` is the 1-based `(x, y)` of its first pixel. For
  sparse files, crop the indices afterwards with `extractROI`.

# Returns
- Dense format: Batched frames array (width × height × nbatches)
- Sparse format: Pixel indices array
"""
function readbin(
    path;
    format=:sparse,
    width::Integer=512,
    height::Integer=512,
    batchsize::Integer=1000,
    crop_origin::Union{NTuple{2,Integer},Nothing}=nothing,
    crop_width::Union{Integer,Nothing}=nothing,
    crop_height::Union{Integer,Nothing}=nothing,
)
    format in (:dense, :sparse) || throw(ArgumentError("format must be :dense or :sparse"))
    format == :sparse && !all(isnothing, (crop_origin, crop_width, crop_height)) && throw(
        ArgumentError("Cropping while reading applies to dense files only; crop sparse indices with `extractROI`."),
    )
    files = checkpath(path)
    @info "Reading $(length(files)) .bin file(s) from $path"

    if format == :dense
        return readbin_dense_directory(
            path;
            width=width,
            height=height,
            batchsize=batchsize,
            crop_origin=something(crop_origin, (1, 1)),
            crop_width=something(crop_width, width),
            crop_height=something(crop_height, height),
        )
    else
        return signal2indices!(_readbin(files), width * height)
    end
end

# ============================================================================
# VISUALIZATION & I/O UTILITIES
# ============================================================================

# function show_frames(frames; pixel_size=0.1, units="um", threshold=nothing, fig_width_px=400, colormap=:grays, title="")
#     """Display frame(s) as a heatmap with spatial coordinates."""
#     img = ndims(frames) == 3 ? dropdims(sum(frames, dims=3), dims=3) : frames
#     H, W = size(img)

#     x = (0:W-1) .* pixel_size .+ pixel_size / 2
#     y = (0:H-1) .* pixel_size .+ pixel_size / 2

#     fig = Figure(size=(fig_width_px, fig_width_px))
#     ax = Axis(fig[1, 1],
#         xlabel="x ($units)",
#         ylabel="y ($units)",
#         aspect=DataAspect(),
#         title=title
#     )
#     if threshold == nothing
#         threshold = maximum(img)
#     end
#     hm = heatmap!(ax, x, y, img, colormap=colormap, colorrange=(0, threshold))
#     Colorbar(fig[1, 2], hm, label="Intensity")
#     return fig
# end

"""
    write_frames_to_tiff(fname, frames; pixel_size=0.1, unit="um", period=1e-5)

Write frame array to TIFF file with metadata.
"""
function write_frames_to_tiff(fname, frames; pixel_size=0.1, unit="um", period=1e-5)
    frames = convert(Array{UInt16,3}, frames)
    writetiff(fname, frames, px_size=pixel_size, unit=unit, period=period)
end

"""
    read_stage_positions(file_path; relative_time=true, relative_position=true, time_unit=:s)

Read stage position CSV and return `(t, x, y)` arrays.

By default `t` is seconds from the first timestamp, and `x`, `y` are microns 
from the first position. Pass `relative_time=false` or `relative_position=false`
to keep absolute values, or `time_unit=:ms` to return milliseconds.
"""
function read_stage_positions(file_path; relative_time=true, relative_position=true, time_unit=:s)
    lines = readlines(file_path)
    isempty(lines) && error("No rows found in $file_path")

    header = strip.(split(strip(first(lines)), ','))
    rows = [strip.(split(strip(line), ',')) for line in lines[2:end] if !isempty(strip(line))]
    isempty(rows) && error("No data rows found in $file_path")

    function column_index(name)
        idx = findfirst(==(name), header)
        isnothing(idx) && error("Missing required column '$name' in $file_path")
        return idx
    end

    t_col = column_index("Timestamp_ms")
    x_col = column_index("X_microns")
    y_col = column_index("Y_microns")
    max_col = max(t_col, x_col, y_col)

    for (row_number, row) in enumerate(rows)
        length(row) >= max_col || error("Row $(row_number + 1) has $(length(row)) columns; expected at least $max_col")
    end

    t = [parse(Float64, row[t_col]) for row in rows]
    x = [parse(Float64, row[x_col]) for row in rows]
    y = [parse(Float64, row[y_col]) for row in rows]

    if relative_time
        t .-= first(t)
    end

    if relative_position
        x .-= first(x)
        y .-= first(y)
    end

    if time_unit in (:s, :seconds)
        t ./= 1000
    elseif !(time_unit in (:ms, :milliseconds))
        error("time_unit must be :s, :seconds, :ms, or :milliseconds")
    end

    return t, x, y
end
