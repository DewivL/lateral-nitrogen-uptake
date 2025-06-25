# this is the main function for the nitrogen uptake
module Nuptake
export n_uptake!


function n_uptake!(N::AbstractArray, INFLUX_MAX, KM, NC_MIN, K2, root_DW::AbstractArray)
    uptake_per_cell = similar(N) # array to store N uptake per cell

    @inbounds for i in eachindex(N)
        # only nitrogen higher than nc_min is available for uptake
        available = max(N[i] - NC_MIN, 0.0)
        hats = (INFLUX_MAX * available) / (KM + available) #μmol/(m2·day)
        lats = K2 * available # μmol/(g day)

        # high affinity depends on area, low affinity depends on root dw. Maximum uptake is what is in the cell.
        uptake = (hats+lats) * root_DW[i]
        uptake = min(uptake, N[i])
        N[i] -= uptake
        uptake_per_cell[i] = uptake # store uptake value
    end

    return uptake_per_cell
end
end
