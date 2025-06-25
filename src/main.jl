module SimulationMain  


include("MyCube.jl")
include("Nmovement.jl")
include("Nuptake.jl")
include("Init.jl")
import VirtualPlantLab as VPL

using CSV, DataFrames, Random, DelimitedFiles, GLMakie, Statistics, ColorTypes, UUIDs
using .MyCube, .Nmovement, .Nuptake, .Init



########################
### File Management  ###
########################

function cleanup_files()
    patterns = [
        r"^output(_fullN|_Nlimited)?\.csv$",
        r"^zone_timeseries(_fullN|_Nlimited)?\.csv$",
        r"^zone_cumuptake(_fullN|_Nlimited)?\.csv$",
        r"^makie_viz_treatment\d+_rep\d+(_t\d+)?\.png$" 
    ]

    for file in readdir()
        if any(p -> occursin(p, file), patterns)
            @info "Removing: $file"
            rm(file; force = true)
        end
    end
end


############################################
### Simulation Configuration & Constants ###
############################################

const custom_DW = [0.8, 1.0, 1.0, 0.9, 0.8, 0.6, 0.4, 0.3, 0.1, 0.1] # this is the total dw in grams per layer in the soil, top- bottom

const treatment = [
    # High soil N, 100% coverage
    (id = 1,  intercrop = false, soilN = 140.0, lateral = true,  influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = 1.0),
    (id = 2,  intercrop = false, soilN = 140.0, lateral = false, influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = 1.0),
    (id = 3,  intercrop = true,  soilN = 140.0, lateral = true,  influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = 1.0),
    (id = 4,  intercrop = true,  soilN = 140.0, lateral = false, influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = 1.0),

    # High soil N, ~70% coverage
    (id = 5,  intercrop = false, soilN = 140.0, lateral = true,  influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = rand([0.68, 0.69, 0.70, 0.71, 0.72])),
    (id = 6,  intercrop = false, soilN = 140.0, lateral = false, influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = rand([0.68, 0.69, 0.70, 0.71, 0.72])),
    (id = 7,  intercrop = true,  soilN = 140.0, lateral = true,  influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = rand([0.68, 0.69, 0.70, 0.71, 0.72])),
    (id = 8,  intercrop = true,  soilN = 140.0, lateral = false, influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = rand([0.68, 0.69, 0.70, 0.71, 0.72])),

    # Low soil N, 100% coverage
    (id = 9,  intercrop = false, soilN = 20.0, lateral = true,  influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = 1.0),
    (id = 10, intercrop = false, soilN = 20.0, lateral = false, influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = 1.0),
    (id = 11, intercrop = true,  soilN = 20.0, lateral = true,  influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = 1.0),
    (id = 12, intercrop = true,  soilN = 20.0, lateral = false, influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = 1.0),

    # Low soil N, ~70% coverage
    (id = 13, intercrop = false, soilN = 20.0, lateral = true,  influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = rand([0.68, 0.69, 0.70, 0.71, 0.72])),
    (id = 14, intercrop = false, soilN = 20.0, lateral = false, influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = rand([0.68, 0.69, 0.70, 0.71, 0.72])),
    (id = 15, intercrop = true,  soilN = 20.0, lateral = true,  influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = rand([0.68, 0.69, 0.70, 0.71, 0.72])),
    (id = 16, intercrop = true,  soilN = 20.0, lateral = false, influx_max = 250.0, NCmin = 1.0, KM = 50, K2 = 1.46, diffusion = 0.00018, coverage = rand([0.68, 0.69, 0.70, 0.71, 0.72]))
]
# here the size of the soil is set in meters
xsize, ysize, zsize = 1.0, 1.0, 1.0

function get_grid_parameters()
    xsize, ysize, zsize = 1.0, 1.0, 1.0
    xmin, xmax = -xsize / 2, xsize / 2
    ymin, ymax = -ysize / 2, ysize / 2

    # create the soil cell boundaries
    xcuts = collect(range(start = (xmin + 0.1), stop = xmax, step = 0.1))
    ycuts = collect(range(start = (ymin + 0.1), stop = ymax, step = 0.1))
    zcuts = collect(0.1:0.1:zsize)

    return length(xcuts), length(ycuts), length(zcuts), MyCube.GridParameters(xcuts, ycuts, zcuts)
end

# this fucntion saves the simulation state as a makie image.
function save_makie_image(N, grid_params, treatment_id, rep; timestep = nothing)
    N_min, N_max = minimum(N), maximum(N)
    # relative colors
    norm_N = (N .- N_min) ./ (N_max - N_min + eps())
    # black grey white colors
    colors = RGBA.(1 .- norm_N, 1 .- norm_N, 1 .- norm_N, 1.0) 
    colors = reverse(colors, dims = 3)
    # render the scene
    graph = VPL.Graph(axiom = MyCube.Tile(zsize, xsize, ysize, colors, grid_params))
    #display(VPL.render(VPL.Mesh(graph)))
    # only needed if want to save the pictures themselves.
    fig = VPL.render(VPL.Mesh(graph))

    filename = isnothing(timestep) ?
        "makie_viz_treatment$(treatment_id)_rep$(rep).png" :
        "makie_viz_treatment$(treatment_id)_rep$(rep)_t$(timestep).png"
    save(filename, fig)
end

########################
###     Output List  ###
########################
# Write all simulation outputs (summary + time series + voxel map)
function write_outputs(treatment, rep, DW_grid, N, total_uptake, root_presence,
    avg_N_legume, avg_N_cereal, avg_N_none,
    zone_uptake_legume, zone_uptake_cereal, zone_uptake_none,
    start_soil_N_layer_string;
    output_suffix = "", coverage= treatment.coverage)

    nx, ny, nz = size(N)
    UMOL_TO_MG = 14.0 / 1000
    total_biomass = sum(DW_grid)
    uptake_legume = sum(total_uptake[root_presence .== 2])
    uptake_cereal  = sum(total_uptake[root_presence .== 1])
    uptake_per_gram_legume_mg = uptake_legume / sum(DW_grid[root_presence .== 2]) * UMOL_TO_MG
    uptake_per_gram_cereal_mg  = uptake_cereal  / sum(DW_grid[root_presence .== 1]) * UMOL_TO_MG
    uptake_per_gram_root_mg = sum(total_uptake) / total_biomass * UMOL_TO_MG
    total_uptake_mg = sum(total_uptake) * UMOL_TO_MG
    final_soil_N = sum(N)

    df_map = DataFrame(x=Int[], y=Int[], z=Int[], N=Float64[], DW=Float64[], uptake=Float64[])
    for x in 1:nx, y in 1:ny, z in 1:nz
    push!(df_map, (x, y, z, N[x, y, z], DW_grid[x, y, z], total_uptake[x, y, z]))
    end

    df = DataFrame([(simulation_id = "Treatment$(treatment.id)_Rep$(rep)",
    treatment_id = treatment.id,
    intercrop = treatment.intercrop,
    soilN = treatment.soilN,
    diffusion = treatment.diffusion,
    lateral_enabled = treatment.lateral,
    influx_max = treatment.influx_max,
    KM = treatment.KM,
    NC_MIN = treatment.NCmin,
    K2 = treatment.K2,
    coverage = treatment.coverage,
    total_root_biomass = total_biomass,
    total_uptake_μmol = sum(total_uptake),
    total_uptake_mg = total_uptake_mg,
    m3_uptake_cereal = uptake_cereal,
    uptake_per_gram_mg = uptake_per_gram_root_mg,
    uptake_per_gram_legume_mg = uptake_per_gram_legume_mg,
    uptake_per_gram_cereal_mg = uptake_per_gram_cereal_mg,
    total_legume_uptake_mg = uptake_legume * UMOL_TO_MG,
    total_cereal_uptake_mg = uptake_cereal * UMOL_TO_MG,
    final_soil_N_μmol = final_soil_N,
    final_soil_N_mg = final_soil_N * UMOL_TO_MG,
    start_soil_N_per_layer = start_soil_N_layer_string)])

    # === Time series data ===
    zone_df = DataFrame(
        treatment_id = fill(treatment.id, length(avg_N_legume)),
        soilN = fill(treatment.soilN, length(avg_N_legume)),   # ADD THIS
        replicate = fill(rep, length(avg_N_legume)),
        timestep = 1:length(avg_N_legume),
        legume_root = avg_N_legume,
        cereal_root = avg_N_cereal,
        no_root = avg_N_none
    )
    
    CSV.write("zone_timeseries$(output_suffix).csv", zone_df;
    append = isfile("zone_timeseries$(output_suffix).csv"),
    writeheader = !isfile("zone_timeseries$(output_suffix).csv"))
    

    # === Cumulative uptake data ===
    cumu_df = DataFrame(
        treatment_id = fill(treatment.id, length(zone_uptake_legume)),
        soilN = fill(treatment.soilN, length(zone_uptake_legume)),  
        replicate = fill(rep, length(zone_uptake_legume)),
        timestep = 1:length(zone_uptake_legume),
        legume_root = cumsum(zone_uptake_legume),
        cereal_root = cumsum(zone_uptake_cereal),
        no_root = cumsum(zone_uptake_none)
    )
    
    CSV.write("zone_cumuptake$(output_suffix).csv", cumu_df;
    append = isfile("zone_cumuptake$(output_suffix).csv"),
    writeheader = !isfile("zone_cumuptake$(output_suffix).csv"))

    # === Summary output ===
    CSV.write("output$(output_suffix).csv", df;
    append = isfile("output$(output_suffix).csv"),
    writeheader = !isfile("output$(output_suffix).csv"))
end

# === SIMULATION RUNNER ===
function run_simulation(treatment, rep, custom_DW,)
    # grid setup
    nx, ny, nz, grid_params = get_grid_parameters()
    DW_grid, N_initial, root_presence = Init.init_soil(
        nx, ny, nz; Nconc = treatment.soilN, intercrop = treatment.intercrop,
        coverage = treatment.coverage, layer_DW_targets = custom_DW)

    N_initial .= treatment.soilN

     # Track initial N per layer for later output
    start_soil_N_layers = [mean(N_initial[:, :, z]) for z in 1:nz]

    # Prepare simulation variables
    N = copy(N_initial)
    N_new = similar(N)
    total_uptake = zeros(size(N))

    avg_N_legume, avg_N_cereal, avg_N_none = Float64[], Float64[], Float64[]
    zone_uptake_legume, zone_uptake_cereal, zone_uptake_none = Float64[], Float64[], Float64[]
    # Used to distinguish full-N vs N-limited treatments in filenames
    output_suffix = treatment.soilN < 100 ? "_Nlimited" : "_fullN"
    # === Optional: initial visualization at t = 0 ===
    # === Before loop, t = 0 ===
    if rep ==1
        #save_makie_image(N, grid_params, treatment.id, rep, timestep = 0)
    end

    for t in 1:40
        uptake = Nuptake.n_uptake!(N, treatment.influx_max, treatment.KM, treatment.NCmin, treatment.K2, DW_grid)
        if treatment.lateral
            Nmovement.n_movement!(N, N_new, treatment.diffusion)
        else
            N_new .= N  
        end
        #if rep == 2 && t == 20
            #ave_makie_image(N, grid_params, treatment.id, rep, timestep = 20)
        #end
        
        # track the N by root type
        push!(avg_N_legume, mean(N[root_presence .== 2]))
        push!(avg_N_cereal, mean(N[root_presence .== 1]))
        push!(avg_N_none, mean(N[root_presence .== 0]))
        
        push!(zone_uptake_legume, sum(uptake[root_presence .== 2]))
        push!(zone_uptake_cereal, sum(uptake[root_presence .== 1]))
        push!(zone_uptake_none, sum(uptake[root_presence .== 0]))
        
        # update for nex timestep
        N, N_new = N_new, N
        total_uptake .+= uptake
    end
    # === optional: final or mid-simulation visualization ===
    #=if rep == 1
        save_makie_image(N, grid_params, treatment.id, rep, timestep = 40)
    end =#
    
    
    #save_makie_image(N, grid_params, treatment.id, rep)
    write_outputs(treatment, rep, DW_grid, N, total_uptake, root_presence, avg_N_legume, avg_N_cereal, avg_N_none, zone_uptake_legume, zone_uptake_cereal, zone_uptake_none, string(start_soil_N_layers), output_suffix = output_suffix, coverage= treatment.coverage)
end

# run all simulations 
function run_all_simulations(num_reps)
    for treatment in treatment
        for rep in 1:num_reps
            run_simulation(treatment, rep, custom_DW)
        end
    end
end


end # module
SimulationMain.cleanup_files()
SimulationMain.run_all_simulations(20)
