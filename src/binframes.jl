function sumbin!(binned::AbstractArray{T,3}, tobin::AbstractArray{S,3}) where {T<:Real,S<:Real}
    b = size(tobin, 3) ÷ size(binned, 3)
    @views for i in axes(binned, 3)
        sum!(
            binned[:, :, i],
            tobin[:, :, (i-1)*b+1:i*b],
        )
    end
    return binned
end

# Returns UInt16, the type of SP2T's readouts; the sums are taken in Int and checked.
function binframes(frames::AbstractArray{T,3}, batchsize::Integer) where {T<:Integer}
    width, height, nframes = size(frames)
    new_nframes = cld(nframes, batchsize)
    binned = similar(frames, Int, width, height, new_nframes)
    n = nframes ÷ batchsize * batchsize
    if n == nframes
        sumbin!(binned, frames)
    else
        # SP2T's detectors assume every frame sums batchsize raw frames
        @warn "The last batch sums $(nframes - n) frames, fewer than the batchsize ($batchsize)."
        @views sumbin!(binned[:, :, 1:end-1], frames[:, :, 1:n])
        @views sum!(binned[:, :, end], frames[:, :, n+1:end])
    end
    lo, hi = extrema(binned)
    0 <= lo && hi <= typemax(UInt16) || throw(ArgumentError(
        "Batch sums range from $lo to $hi, outside 0 to $(typemax(UInt16)), which SP2T's " *
        "UInt16 readouts can hold: use a smaller batchsize."))
    return convert(AbstractArray{UInt16}, binned)
end

# The mean position over each batch of frames, a partial last batch included, so the result
# matches binframes frame for frame. No warning: the mean of the frames left over is still
# the position for the last binned frame (binframes warns about that frame).
function bintracks(tracks::AbstractArray{<:AbstractFloat,3}, batchsize::Integer)
    nframes = size(tracks, 1)
    binned = similar(tracks, cld(nframes, batchsize), size(tracks, 2), size(tracks, 3))
    @views for i in axes(binned, 1)
        batch = (i-1)*batchsize+1:min(i * batchsize, nframes)
        mean!(binned[i:i, :, :], tracks[batch, :, :])
    end
    return binned
end

repeattracks(tracks::AbstractArray{<:Real,3}, batchsize::Integer) = repeat(tracks, inner=(batchsize, 1, 1))