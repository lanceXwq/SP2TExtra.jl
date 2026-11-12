# XMLDict gives an element that repeats as a vector, a single one as itself, and none at all
# as a missing key; list them the same way in each case.
function xmlelements(dict::XMLDict.XMLDictElement, key::String)
    haskey(dict, key) || return XMLDict.XMLDictElement[]
    elements = dict[key]
    return elements isa AbstractVector ? elements : [elements]
end

function parsedetection!(x::AbstractArray{<:Real,2}, particle::XMLDict.XMLDictElement, batchsize::Integer)
    for detection in xmlelements(particle, "detection")
        trange = tryparse(Int, detection[:t]) * batchsize .+ (1:batchsize)
        for (dim, coord) in enumerate([:x, :y, :z])
            x[trange, dim] .= tryparse(Float64, detection[coord])
        end
    end
    return x
end

function parseparticle!(x::AbstractArray{<:Real,3}, dict::XMLDict.XMLDictElement, batchsize::Integer)
    for (m, particle) in enumerate(xmlelements(dict, "particle"))
        parsedetection!(view(x, :, :, m), particle, batchsize)
    end
    return x
end

function xml2tracks(xml::String; batchsize::Integer=1, nframes::Integer)
    xmlcontent = parse_xml(String(read(xml)))
    ntracks = tryparse(Int, xmlcontent[:nTracks])
    tracks = fill!(Array{Float64}(undef, nframes, 3, ntracks), NaN)
    parseparticle!(tracks, xmlcontent, batchsize)
end

function xml2tracks(xmls::AbstractVector{String}; batchsizes::AbstractVector{<:Integer}, nframes::Integer)
    tracks = Vector{Array{Float64,3}}()
    for (xml, batchsize) in zip(xmls, batchsizes)
        push!(tracks, xml2tracks(xml; batchsize=batchsize, nframes=nframes))
    end
    tracks
end