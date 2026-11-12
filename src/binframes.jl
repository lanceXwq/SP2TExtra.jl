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

function binframes(frames::AbstractArray{T,3}, batchsize::Integer) where {T<:Integer}
    S = maximum(frames) * batchsize > typemax(T) ? UInt16 : T
    width, height, nframes = size(frames)
    new_nframes = cld(nframes, batchsize)
    binned = similar(frames, S, width, height, new_nframes)
    n = nframes ÷ batchsize * batchsize
    if n == nframes
        sumbin!(binned, frames)
    else
        @views sumbin!(binned[:, :, 1:end-1], frames[:, :, 1:n])
        @views sum!(binned[:, :, end], frames[:, :, n+1:end])
    end
    return binned
end

function meanbin!(binned::AbstractArray{T,3}, tobin::AbstractArray{T,3}, batchsize::Integer) where {T<:AbstractFloat}
    @views for i in axes(binned, 1)
        mean!(
            binned[i:i, :, :],
            tobin[(i-1)*batchsize+1:i*batchsize, :, :],
            weights(ones(T, batchsize)),
            dims=1,
        )
    end
    return binned
end

function bintracks(tracks::AbstractArray{<:AbstractFloat,3}, batchsize::Integer)
    binned = similar(
        tracks,
        size(tracks, 1) ÷ batchsize,
        size(tracks, 2),
        size(tracks, 3),
    )
    meanbin!(binned, tracks, batchsize)
    return binned
end

repeattracks(tracks::AbstractArray{<:Real,3}, batchsize::Integer) = repeat(tracks, inner=(batchsize, 1, 1))