"""
    tamsd(tracks; lags = 1:size(tracks, 1)-1)

Time-averaged mean squared displacement of each particle in `tracks`
(`(nframes, ndims, nparticles)`), as a `length(lags)` × `nparticles` matrix. The entry
for lag `L` averages `(x[t+L] - x[t])²` over every frame `t` and every coordinate, so it is
per coordinate like SP2T's MSD: for free diffusion it is about `L` times the one-frame MSD.
"""
function tamsd(tracks::AbstractArray{<:Real,3}; lags=1:size(tracks, 1)-1)
    nframes, ndims, nparticles = size(tracks)
    all(in(1:nframes-1), lags) ||
        throw(ArgumentError("lags must lie in 1:$(nframes - 1) for $nframes frames."))
    m = similar(tracks, float(eltype(tracks)), length(lags), nparticles)
    for (i, L) in enumerate(lags)
        Δ = @views tracks[1+L:end, :, :] .- tracks[1:end-L, :, :]
        m[i, :] .= vec(sum(abs2, Δ; dims=(1, 2))) ./ ((nframes - L) * ndims)
    end
    return m
end

"""
    msd(tracks)

Mean squared displacement between consecutive frames, averaged over frames, coordinates and
particles of `tracks` (`(nframes, ndims, nparticles)`): the per-coordinate MSD that SP2T
samples (`2·D·period`). For the sampled tracks of a chain, use `msd(sampledtracks(chain))`
(SP2T); the sampled MSDs themselves are `chain.msds`.
"""
msd(tracks::AbstractArray{<:Real,3}) = mean(abs2, diff(tracks; dims=1))
