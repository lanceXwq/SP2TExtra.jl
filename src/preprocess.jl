"""
    getdarkcounts(darkframes; batchsize = 1)

Dark counts of a SPAD array, as SP2T expects them: the expected number of dark photons per
binary frame in each pixel. `darkframes` (`(width, height, nframes)`) are frames recorded
without emitters, each summing `batchsize` binary frames. A pixel that clicks in `k` of
the `N = nframes·batchsize` binary frames clicks with probability `1 - exp(-λ)`, so its
maximum-likelihood dark count is `λ = -log(1 - k/N)`. A pixel that clicks in every binary
frame gets `Inf` (with a warning); SP2T filters such pixels out.
"""
function getdarkcounts(darkframes::AbstractArray{<:Integer,3}; batchsize::Integer=1)
    all(r -> 0 <= r <= batchsize, darkframes) || throw(ArgumentError(
        "Readouts must lie between 0 and the batchsize ($batchsize): is the batchsize right?"))
    darkcounts = Array{Float64}(undef, size(darkframes, 1), size(darkframes, 2))
    sum!(darkcounts, darkframes)
    N = size(darkframes, 3) * batchsize
    @. darkcounts .= -log1p(-darkcounts / N)
    nsaturated = count(isinf, darkcounts)
    nsaturated > 0 &&
        @warn "$nsaturated pixel(s) clicked in every binary frame: their dark counts are Inf, and SP2T filters them out."
    return darkcounts
end

# function getframes() end