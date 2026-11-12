function credible1D(S::AbstractVector{<:Sample}, i::Integer, xedges::AbstractRange, yedges::AbstractRange; factor::Real=1, xshift::Real=0, yshift::Real=0)
    N = size(S[1].tracks, 1)
    xcounts = zeros(Float64, N, length(xedges) - 1)
    ycounts = zeros(Float64, N, length(yedges) - 1)
    ntracks = length(S)
    for s in S
        @views for (n, x) in enumerate(eachslice(s.tracks, dims=1))
            histcounts!(xcounts[n, :], (x[1, i:i] .+ xshift) .* factor, xedges)
            histcounts!(ycounts[n, :], (x[2, i:i] .+ yshift) .* factor, yedges)
        end
    end
    xcounts ./= ntracks
    ycounts ./= ntracks
    xcounts, ycounts
end

function credible1D(S::AbstractVector{<:Sample}, xedges::AbstractRange, yedges::AbstractRange; factor::Real=1, xshift::Real=0, yshift::Real=0)
    N = size(S[1].tracks, 1)
    xcounts = zeros(Float64, N, length(xedges) - 1)
    ycounts = zeros(Float64, N, length(yedges) - 1)
    ntracks = 0
    for s in S
        ntracks += size(s.tracks, 3)
        @views for (n, x) in enumerate(eachslice(s.tracks, dims=1))
            histcounts!(xcounts[n, :], (x[1, :] .+ xshift) .* factor, xedges)
            histcounts!(ycounts[n, :], (x[2, :] .+ yshift) .* factor, yedges)
        end
    end
    xcounts ./= ntracks
    ycounts ./= ntracks
    xcounts, ycounts
end

function credible2D(S::AbstractVector{<:Sample}, i::Integer, xedges::AbstractRange, yedges::AbstractRange; factor::Real=1, xshift::Real=0, yshift::Real=0)
    N = size(S[1].tracks, 1)
    counts = zeros(Float64, length(xedges) - 1, length(yedges) - 1, N)
    ntracks = length(S)
    @views for s in S
        for (n, x) in enumerate(eachslice(s.tracks, dims=1))
            histcounts!(counts[:, :, n], (x[2, i:i] .+ yshift) .* factor, (x[1, i:i] .+ xshift) .* factor, yedges, xedges)
        end
    end
    counts ./= ntracks
end

function credible2D(S::AbstractVector{<:Sample}, xedges::AbstractRange, yedges::AbstractRange; factor::Real=1, xshift::Real=0, yshift::Real=0)
    N = size(S[1].tracks, 1)
    counts = zeros(Float64, length(xedges) - 1, length(yedges) - 1, N)
    ntracks = 0
    @views for s in S
        ntracks += size(s.tracks, 3)
        for (n, x) in enumerate(eachslice(s.tracks, dims=1))
            histcounts!(counts[:, :, n], (x[2, :] .+ yshift) .* factor, (x[1, :] .+ xshift) .* factor, yedges, xedges)
        end
    end
    counts ./= ntracks
end

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
samples (`2·D·period`). For the sampled tracks of a chain, use `msd(tracks(chain))`; the
sampled MSDs themselves are `chain.msds`.
"""
msd(tracks::AbstractArray{<:Real,3}) = mean(abs2, diff(tracks; dims=1))

function tracks(chain::Chain{T}; burn_in::Integer=0) where {T<:Real}
    N = sum(chain.emittercounts[burn_in+1:end])
    t = chain.samples[1].tracks
    x = Array{T}(undef, size(t, 1), size(t, 2), N)
    i = 1
    for s in chain.samples[burn_in+1:end]
        n = size(s.tracks, 3)
        copyto!(view(x, :, :, i:i+n-1), s.tracks)
        i += n
    end
    return x
end

"""
    localization_error(chain; burn_in = 0)

Localization uncertainty of a single emitter: the posterior standard deviation of its
position per coordinate, as the root mean square over frames and coordinates,
`sqrt(mean(var))`, the variance taken over the samples after the first `burn_in`. This is
the σ in the offset `2σ²` of the per-coordinate MSD curve (`tamsd`). It measures the spread
of the posterior, not the distance to a true position. Every sample after `burn_in` must
have exactly one emitter.
"""
function localization_error(chain::Chain; burn_in::Integer=0)
    counts = chain.emittercounts[burn_in+1:end]
    !isempty(counts) && all(==(1), counts) || throw(ArgumentError(
        "localization_error needs exactly one emitter in every sample after burn_in; " *
        "the emitter counts there are $(sort(unique(counts)))."))
    x = tracks(chain; burn_in=burn_in)
    return sqrt(mean(var(x; dims=3, corrected=false)))
end