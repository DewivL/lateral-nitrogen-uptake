# create the structure to vizualize the soil in Makie, based on VPL tutorial code (using the slicer)
module MyCube

using VirtualPlantLab
import ColorTypes: RGBA  

# grid setup 
struct GridParameters
    xcuts::Vector{Float64}  # X slice positions
    ycuts::Vector{Float64}  # Y slice positions
    zcuts::Vector{Float64}  # Z slice positions
end

# cube dimensions 
"A tile can be made into a VPL cube."
struct Tile <: Node
    length::Float64   # in m (x-dimension)
    width::Float64    # in m (y-dimension)
    height::Float64   # in m(z-dimension)
    colors::Array{RGBA{Float64},3}  
    grid_params::GridParameters     # contains slice positions in 3D
end    

function VirtualPlantLab.feed!(turtle::Turtle, t::Tile, data)
    # create a base solid cube with specified dimensions
    r = SolidCube(turtle, length=t.length, width=t.width, height=t.height)

    # apply slicing based on the cut planes
    slice!(r, X=t.grid_params.xcuts, Y=t.grid_params.ycuts, Z=t.grid_params.zcuts)

    # count how many voxel cells were created by the slicing
    num_slices = length(properties(r)[:slices])

    # prepare a color array
    all_colors = Vector{RGBA{Float64}}(undef, num_slices)

    # assign each voxel a color based on its position
    for (i, slice) in enumerate(properties(r)[:slices])
        x, y, z = slice[1], slice[2], slice[3]

        # avoid out-of-bounds errors by clamping
        x = clamp(x, 1, size(t.colors, 1))
        y = clamp(y, 1, size(t.colors, 2))
        z = clamp(z, 1, size(t.colors, 3))

        # assign the corresponding color
        all_colors[i] = t.colors[x, y, z]
    end

    # color the cube mesh based on the color array
    Mesh!(turtle, r, colors=all_colors)

    return nothing
end

export GridParameters, Tile

end  
