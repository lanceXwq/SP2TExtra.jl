# Dark count calibration: estimate each pixel's dark counts (expected dark photons per binary
# frame, as SP2T's detectors take them) from a sparse SPAD recording made without emitters,
# within the region that the tracking data will use, and save them for SP2T.
#
# Needs SP2TExtra, JLD2 and GLMakie in the active environment.
using SP2TExtra
using JLD2
using GLMakie

# Settings: edit these for your recording.
datapath = "path/to/dark/recording/" # a directory of sparse .bin files, or one file
fullframe = (512, 512) # width and height of the recorded frames
roiorigin = (198, 337, 201 * 255 + 1) # first pixel (x, y) and first frame of the region
roisize = (45, 45, 20 * 255) # its width, height and number of frames
batchsize = 255 # raw frames summed per frame; any value gives the same dark counts
output = "darkcounts.jld2"

indices = readbin(datapath; width = fullframe[1], height = fullframe[2])
roiindices = extractROI(indices, fullframe, roiorigin, roisize)
# nframes keeps the region's last frames even if they have no clicks: dropping them would
# count too few frames and overestimate the dark counts
darkframes = getframes(
    roiindices;
    width = roisize[1],
    height = roisize[2],
    batchsize = batchsize,
    nframes = roisize[3],
)
# Inf (with a warning) for a pixel that clicks in every frame; SP2T filters those out
darkcounts = getdarkcounts(darkframes; batchsize = batchsize)

fig = Figure()
ax = Axis(fig[1, 1]; aspect = DataAspect(), title = "dark counts per binary frame")
hm = heatmap!(ax, replace(darkcounts, Inf => NaN); colormap = :bone, nan_color = :red)
Colorbar(fig[1, 2], hm)
display(fig) # pixels that click in every frame are red

jldsave(output; darkcounts)
