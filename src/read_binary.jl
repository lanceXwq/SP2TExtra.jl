# read_binary.jl
#
# Reading SPAD recordings from .bin files, in dense and sparse formats, into batched frames.
# - Dense: one bit per pixel
# - Sparse: the indices of the pixels that clicked, as 32-bit integers

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
# DENSE FORMAT - Reading bits into frames in one pass
# ============================================================================

"""
    readdense(files, width, height, batchsize, roiorigin, roisize; chunkbytes = 2^20)

Internal: the batched frames of the region of dense `files` and its number of binary
frames (`roisize[3]`, or `nothing` for the rest of the recording), in one pass. Each file
is read about `chunkbytes` at a time, in blocks of whole frames that are also whole bytes,
and each region pixel's bit (most significant first) is added to its batch. When the
region's rows are whole bytes (`width`, `x - 1` and the region's width multiples of 8),
whole bytes are added at once, 8 pixels per operation. Files and frames outside the region
are skipped.
"""
function readdense(files, width, height, batchsize, roiorigin, roisize; chunkbytes=2^20)
    framesize = width * height
    file_sizes_in_bits = filesize.(files) * 8
    n_frames_per_file = file_sizes_in_bits .÷ framesize
    n_bits_left_per_file = rem.(file_sizes_in_bits, framesize)
    any(n_bits_left_per_file .!= 0) && throw(ArgumentError("Some files have extra bits that do not form a complete $width × $height frame. File sizes in bits: $file_sizes_in_bits, bits left per file: $n_bits_left_per_file"))

    x0, y0, f0 = roiorigin
    nframes = something(roisize[3], sum(n_frames_per_file) - f0 + 1)
    lastframe = f0 + nframes - 1
    (nframes >= 1 && lastframe <= sum(n_frames_per_file)) || throw(ArgumentError(
        "The region's frames ($f0 to $lastframe) are not all in the recording's $(sum(n_frames_per_file)) frames."))

    frames = zeros(UInt16, roisize[1], roisize[2], cld(nframes, batchsize))
    unit = 8 ÷ gcd(framesize, 8) # frames in the smallest block of whole bytes
    buffer = Vector{UInt8}(undef, unit * max(1, chunkbytes * 8 ÷ (framesize * unit)) * framesize ÷ 8)
    wholebytes = width % 8 == 0 && (x0 - 1) % 8 == 0 && roisize[1] % 8 == 0
    counts = zeros(UInt64, wholebytes ? roisize[1] * roisize[2] ÷ 8 : 0)
    before = 0 # frames in the files already passed
    @showprogress dt = 3 desc = "Reading binary image data" for (file, n) in zip(files, n_frames_per_file)
        first, last = max(1, f0 - before), min(n, lastframe - before) # this file's frames in the region
        first <= last && adddense!(frames, counts, wholebytes, buffer, file, first, last,
            before - f0 + 1, width, framesize, x0, y0, batchsize, unit)
        before += n
    end
    return frames, nframes
end

# Add the bits of frames first to last of a dense file to their batches of frames (frame f
# of the file is frame f + shift of the region), reading blocks of frames into buffer. With
# wholebytes, a frame's bits go first to counts, 8 one-byte counters per UInt64, which are
# added to the batch every 255 frames, at the batch's end and at the file's end.
function adddense!(frames, counts, wholebytes, buffer, file, first, last, shift, width,
    framesize, x0, y0, batchsize, unit)
    perchunk = length(buffer) * 8 ÷ framesize # a multiple of unit
    frame = (first - 1) ÷ unit * unit # frames skipped, a whole number of bytes
    incounts = 0 # frames in counts
    io = open(file)
    try
        seek(io, frame * framesize ÷ 8)
        while frame < last
            n = min(perchunk, last - frame)
            read!(io, view(buffer, 1:cld(n * framesize, 8)))
            for m in 1:n
                f = frame + m # the frame of the file
                f < first && continue
                r = f + shift # the frame of the region
                batch = view(frames, :, :, (r - 1) ÷ batchsize + 1)
                if wholebytes
                    addbytes!(counts, buffer, ((m - 1) * framesize + x0 - 1) ÷ 8, width ÷ 8,
                        size(frames, 1) ÷ 8, y0)
                    incounts += 1
                    if incounts == 255 || r % batchsize == 0 || f == last
                        flushcounts!(batch, counts)
                        incounts = 0
                    end
                else
                    addbits!(batch, buffer, (m - 1) * framesize + x0 - 1, width, y0)
                end
            end
            frame += n
        end
    finally
        close(io)
    end
    return frames
end

# Add a frame's region to batch bit by bit; the region's first bit is bit `base` (0-based)
# of buffer.
function addbits!(batch, buffer, base, width, y0)
    @inbounds for j in axes(batch, 2)
        row = base + (y0 + j - 2) * width
        for i in axes(batch, 1)
            b = row + i - 1
            batch[i, j] += (buffer[b>>3+1] >> (7 - (b & 7))) & 0x01
        end
    end
    return batch
end

# Add a frame's region to counts byte by byte, 8 pixels at once; the region's first byte is
# byte `base` (0-based) of buffer, and its rows are `rowbytes` bytes long.
function addbytes!(counts, buffer, base, widthbytes, rowbytes, y0)
    @inbounds for j in 1:length(counts)÷rowbytes
        row = base + (y0 + j - 2) * widthbytes
        c = (j - 1) * rowbytes
        @simd for q in 1:rowbytes
            counts[c+q] += spreadbits(buffer[row+q])
        end
    end
    return counts
end

# A byte's 8 bits, most significant first, as the 8 bytes of a UInt64 (the first bit in the
# lowest byte), each 0 or 1.
function spreadbits(b::UInt8)
    x = (b * 0x0101010101010101) & 0x0102040810204080 # byte l keeps bit 7 - l
    return ((x + 0x7f7f7f7f7f7f7f7f) >> 7) & 0x0101010101010101 # nonzero bytes to 1
end

# Add the one-byte counters to batch, pixel by pixel, and zero them.
function flushcounts!(batch, counts)
    @inbounds for w in eachindex(counts)
        x = counts[w]
        for l in 0:7
            batch[8(w-1)+l+1] += (x >> (8l)) % UInt8
        end
        counts[w] = 0
    end
    return batch
end

# ============================================================================
# SPARSE FORMAT - Reading clicks into frames in one pass
# ============================================================================

"""
    readsparse(files, width, height, batchsize, roiorigin, roisize; chunk = 2^20)

Internal: the batched frames of the region of sparse `files` and its number of binary
frames, in one pass. The files are read in chunks of `chunk` values. A delimiter (a value
of at least `width·height`) ends a frame, so the frame of every click is known without
storing it; a table maps each pixel to its index in the region (0 outside), and each click
in the region adds 1 to its pixel in the current batch. Without `roisize[3]` the frames grow
until the recording ends; with it, reading stops after the region's last frame.
"""
function readsparse(files, width, height, batchsize, roiorigin, roisize; chunk=2^20)
    framesize = width * height
    sizes = filesize.(files)
    corrupted = sizes .% 4 .!= 0
    any(corrupted) && throw(DomainError(sizes[corrupted],
        "File size is not divisible by 4 in $(files[corrupted])"))

    x0, y0, f0 = roiorigin
    rw, rh, nframes = roisize
    npixels = rw * rh
    region = zeros(Int32, framesize) # 0-based pixel + 1 -> index in the region, 0 outside
    for j in 1:rh, i in 1:rw
        region[(y0+j-2)*width+x0+i-1] = (j - 1) * rw + i
    end
    lastframe = isnothing(nframes) ? typemax(Int) : f0 + nframes - 1
    frames = zeros(UInt16, npixels * (isnothing(nframes) ? 1 : cld(nframes, batchsize)))

    buffer = Vector{Int32}(undef, chunk)
    frame = 1 # the current frame of the recording
    nread = 0 # frames of the region read so far
    offset = 0 # where the current batch starts in frames
    delimited = true # whether the last value read was a delimiter
    for (file, nvalues) in zip(files, sizes .÷ 4)
        frame > lastframe && break
        io = open(file)
        try
            while nvalues > 0 && frame <= lastframe
                n = min(chunk, nvalues)
                read!(io, view(buffer, 1:n))
                nvalues -= n
                for k in 1:n
                    s = @inbounds buffer[k]
                    if s >= framesize # a delimiter ends the frame
                        frame >= f0 && (nread += 1)
                        frame += 1
                        frame > lastframe && break
                        offset = nread ÷ batchsize * npixels
                        if offset + npixels > length(frames) # only without nframes
                            oldlength = length(frames)
                            resize!(frames, max(2oldlength, offset + npixels))
                            fill!(view(frames, oldlength+1:length(frames)), 0)
                        end
                    elseif s >= 0
                        if frame >= f0
                            r = @inbounds region[s+1]
                            r != 0 && @inbounds frames[offset+r] += one(UInt16)
                        end
                    else
                        throw(ArgumentError("$file holds a negative value, $s."))
                    end
                end
                delimited = frame > lastframe || buffer[n] >= framesize
            end
        finally
            close(io)
        end
    end
    if !delimited # the recording stops mid-frame
        @warn "The last signal is not a frame delimiter: the clicks after the last delimiter are kept as a last frame, which may be incomplete."
        frame >= f0 && (nread += 1)
        frame += 1
    end

    if isnothing(nframes)
        nread >= 1 || throw(ArgumentError(
            "The region starts at frame $f0, after the recording's $(frame - 1) frames."))
    elseif nread < nframes
        throw(ArgumentError(
            "The region's frames ($f0 to $lastframe) are not all in the recording's $(frame - 1) frames."))
    end
    nbatches = cld(nread, batchsize)
    return reshape(resize!(frames, nbatches * npixels), rw, rh, nbatches), nread
end

# ============================================================================
# UNIFIED API - Main Entry Point
# ============================================================================

"""
    readbin(path; format=:sparse, width=512, height=512, batchsize=1000, roiorigin=(1, 1, 1), roisize=nothing, threshold=typemax(Int))

Read a SPAD recording from a `.bin` file, or from the `.bin` files of a directory (read in
name order, as one stream), and return the frames of a region: `(roiwidth, roiheight,
nbatches)` `UInt16` counts, each the sum of `batchsize` binary frames.

# Formats
- `:sparse` (default): `Int32` values. Each frame lists the 0-based indices of the pixels
  that clicked, then a delimiter (`width·height`). Frames may continue in the next file.
  Clicks after the last delimiter are kept as a last frame, with a warning.
- `:dense`: one bit per pixel, most significant bit first, frame after frame, in whole
  frames per file.

# Keywords
- `width`, `height`: the size of the recorded frames.
- `batchsize`: binary frames summed per frame, from 1 to 65535 (so the sums fit `UInt16`).
  A last batch of fewer binary frames is kept, with a warning.
- `roiorigin = (x, y, frame)`: the region's first pixel and first binary frame, 1-based.
- `roisize = (width, height, nframes)`: the region's size in pixels and binary frames.
  `(width, height)` reads to the end of the recording, and `nothing` (default) also to the
  edges of the frame. A region outside the frame or the recording is an `ArgumentError`.
- `threshold`: counts above it are set to 0.
"""
function readbin(
    path;
    format=:sparse,
    width::Integer=512,
    height::Integer=512,
    batchsize::Integer=1000,
    roiorigin::NTuple{3,Integer}=(1, 1, 1),
    roisize::Union{Nothing,NTuple{2,Integer},NTuple{3,Integer}}=nothing,
    threshold::Integer=typemax(Int),
)
    format in (:dense, :sparse) || throw(ArgumentError("format must be :dense or :sparse"))
    1 <= batchsize <= typemax(UInt16) || throw(ArgumentError(
        "batchsize must lie between 1 and 65535, so that the summed counts fit UInt16."))
    x0, y0, f0 = roiorigin
    rw, rh = isnothing(roisize) ? (width - x0 + 1, height - y0 + 1) : roisize[1:2]
    nframes = roisize isa NTuple{3,Integer} ? roisize[3] : nothing
    (min(x0, y0, f0, rw, rh, something(nframes, 1)) >= 1 && x0 + rw - 1 <= width &&
     y0 + rh - 1 <= height) || throw(ArgumentError(
        "The region (roiorigin = $roiorigin, roisize = $roisize) must lie within the $width × $height frames."))
    files = checkpath(path)
    @info "Reading $(length(files)) .bin file(s) from $path"

    reader = format == :dense ? readdense : readsparse
    frames, nread = reader(files, width, height, batchsize, (x0, y0, f0), (rw, rh, nframes))
    framesleft = rem(nread, batchsize)
    framesleft != 0 &&
        @warn "The last batch sums $framesleft frames, fewer than the batchsize ($batchsize)."
    threshold < batchsize && map!(c -> c > threshold ? zero(c) : c, frames, frames)
    return frames
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
