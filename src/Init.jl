# here, the roots and the soil is initialized
module Init
using Random, CSV, Statistics

export init_soil, save_root_distribution!

# roots, creates a new root distribution every run, every layer.
function generate_root_distribution(nx, ny, nz; coverage, layer_DW_targets, intercrop = false)
    DW_grid = zeros(Float64, nx, ny, nz)
    root_presence = zeros(Int, nx, ny, nz)  # 0 = no root (default)
    voxels_per_layer = nx * ny
    active_voxels_per_layer = round(Int, coverage * voxels_per_layer) 

    for z in 1:nz # for each layer
        layer_indices = [idx for idx in CartesianIndices(DW_grid) if idx[3] == z]
        selected = Random.shuffle(layer_indices)[1:active_voxels_per_layer] # choose random voxel for root
        DW_per_voxel = layer_DW_targets[z] / active_voxels_per_layer # distribute the layer dw over the voxels

        for idx in selected
            if intercrop
                # Left half = cereal (1), right half = legume (2)
                if idx[1] <= nx ÷ 2
                    root_presence[idx] = 1  # cereal
                    crop_dw_mult = 1.5 # increase cereal root dw
                else 
                    root_presence[idx] = 2  # legume
                    crop_dw_mult = 1.0
                end
            else 
                root_presence[idx] = 1  # all cereal (if not intercrop)
                crop_dw_mult = 1.5
            end
            noise = 1.1 + 0.1 * (rand() - 0.5) # add 5% variation on the dw
            DW_grid[idx] = DW_per_voxel * crop_dw_mult * noise
        end
    end
    return DW_grid, root_presence
end

# soil init with the roots
function init_soil(nx, ny, nz; Nconc, intercrop = false, coverage, layer_DW_targets = nothing)
    # if no DW targets given, assign uniform target of 1.0 per layer
    layer_DW_targets = layer_DW_targets === nothing ? fill(1.0, nz) : layer_DW_targets
    
    # call the root distribution function
    DW_grid, root_presence = generate_root_distribution(nx, ny, nz; coverage, layer_DW_targets, intercrop)

    # fill the N concentration as well. now soil is ready for simulation.
    N_initial = fill(Nconc, nx, ny, nz)

    return DW_grid, N_initial, root_presence
end

end
