"""
    viewframes(frames; batchsize = maximum(frames), tracks = nothing, pxsize = 1)

Show `frames`, `(width, height, nframes)` integer readouts, in an interactive viewer
with sliders for the frame, the contrast (up to `batchsize`) and a region of interest,
and overlay `tracks`, `(nframes, ndims, nparticles)` positions in the units of `pxsize`,
drawn up to the current frame. As in SP2T, pixel `(i, j)` covers
`[(i - 1)·pxsize, i·pxsize] × [(j - 1)·pxsize, j·pxsize]`; y points down, as in images.
The region is labeled with its first pixel `(x, y)` and its size, as `extractROI` takes
them.

Defined by an extension: load a Makie backend first (`using GLMakie` for an interactive
window).
"""
function viewframes end
