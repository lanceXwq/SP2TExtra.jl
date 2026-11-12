module SP2TExtraMakieExt

using SP2TExtra
using Makie

theme_dataviewer = Theme(
    Axis = (
        aspect = DataAspect(),
        bottomspinevisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        xticklabelsvisible = false,
        xticksvisible = false,
        ygridvisible = false,
        yticklabelsvisible = false,
        yticksvisible = false,
    ),
    Poly = (; strokecolor = (:deepskyblue3, 0.5), strokewidth = 2),
    Heatmap = (; colormap = :bone),
)

function validaterange(s::String)
    parsed_vec = tryparse.(Int64, split(s, ','))
    return eltype(parsed_vec) == Int64 && length(parsed_vec) == 2
end

getrect(w_sl, h_sl) = @lift Rect(
    $(w_sl.interval)[1] - 0.5,
    $(h_sl.interval)[2] + 0.5,
    $(w_sl.interval)[2] - $(w_sl.interval)[1] + 1,
    $(h_sl.interval)[1] - $(h_sl.interval)[2],
)

getvertices(rect) = (@lift $(rect).origin), (@lift $(rect).origin .+ $(rect).widths)

function SP2TExtra.viewframes(
    frames::AbstractArray{<:Integer,3};
    batchsize::Union{Integer,Nothing} = nothing,
    tracks::Union{Nothing,AbstractArray} = nothing,
    pxsize::Real = 1,
)
    width, height, nframes = size(frames)
    isnothing(batchsize) && (batchsize = maximum(frames))
    isnothing(tracks) && (tracks = fill(NaN, (size(frames, 3), 3, 1)))

    set_theme!(theme_dataviewer)
    fig = Figure(; size = (1000, 300))
    ax = Axis(fig[1, 1])

    h_sl = IntervalSlider(
        fig[1, 2],
        range = 0:height,
        startvalues = (0, height),
        horizontal = false,
    )
    w_sl = IntervalSlider(fig[2, 1], range = 0:width, startvalues = (0, width))
    slgrid = SliderGrid(
        fig[3, :][1, 1],
        (
            label = "Frame",
            range = 1:nframes,
            startvalue = 1,
            snap = false,
            format = x -> "$x/$nframes",
        ),
        (
            label = "Contrast",
            range = 1:batchsize,
            startvalue = batchsize,
            snap = false,
            format = x -> "$x/$batchsize",
        ),
    )

    colorrange = @lift (0, $(slgrid.sliders[2].value))

    frameobs = @lift rotr90(view(frames, :, :, $(slgrid.sliders[1].value)))

    tracks2 = tracks ./ pxsize
    trajs2D = []
    for x in eachslice(view(tracks2, :, 1:2, :), dims = 3)
        push!(trajs2D, Point2f[eachslice(x, dims = 1)...])
    end
    trajs2Dobs = []
    foreach(x -> push!(trajs2Dobs, @lift view(x, 1:$(slgrid.sliders[1].value))), trajs2D)


    heatmap!(ax, 1:width, 1:height, frameobs, colorrange = colorrange)
    limits!(ax, 0.5, width + 0.5, 0.5, height + 0.5)

    foreach(x -> lines!(ax, x), trajs2Dobs)

    @show rect = getrect(w_sl, h_sl)
    poly!(ax, rect, color = (:white, 0))
    # @show w_sl.interval, rect[].origin
    # sl_tags = (@lift rect[].origin)
    # (@lift "($($(w_sl.interval)[2]-$(w_sl.interval)[1]+1), $($(h_sl.interval)[2]-$(h_sl.interval)[1]+1))")

    # @show rect[].orgin = 1

    # txts = @show @lift "($($(rect).origin[1] - 0.5), $($(rect).origin[2]+$(rect).widths[2]) )"
    topleft_txt =
        @lift "($(ceil(Int, $(rect).origin[1])), $(floor(Int, $(rect).origin[2])))"
    bottomright_txt = @lift "($($(rect).widths[1]), $($(rect).widths[2]))"

    @show verts = getvertices(rect)
    text!(ax, verts[1], align = (:left, :top), text = topleft_txt, color = :white)
    text!(ax, verts[2], align = (:right, :bottom), text = bottomright_txt, color = :white)

    Label(fig[0, 1:2], "Frame viewer")

    rowsize!(fig.layout, 1, Aspect(1, height / width))
    resize_to_layout!(fig)
    display(fig)
    set_theme!()
end

end
