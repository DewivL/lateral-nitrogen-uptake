# this is the main function for the nitrogen movement 

module Nmovement
export n_movement!

using Printf


function n_movement!(N::AbstractArray, N_new::AbstractArray, diffusion_coefficient)
    nx, ny, nz = size(N) # get the size of the soil 
    @assert size(N) == size(N_new)
    net_flux = zeros(size(N))
    N_new .= N  # Copy current state

    @inbounds for k in 1:nz # loop over the coordinates
        for j in 1:ny
            for i in 1:nx
                for (dx, dy) in ((1,0), (-1,0), (0,1), (0,-1)) # check vertical neighbors
                    ni, nj = i + dx, j + dy

                    # skip out-of-bounds neighbors
                    if 1 ≤ ni ≤ nx && 1 ≤ nj ≤ ny

                        if N[ni, nj, k] < N[i, j, k] # only diffuse if neighbor has lower conc
                            flux = diffusion_coefficient * ((N[i, j, k] - N[ni, nj, k]) / (0.1^2))  # Fick's law, distance between two cells 0.1 meter
                            if flux > 0
                                flux = min(flux, N[i, j, k])  # can't send more than avaialable
                            else
                                flux = max(flux, -N[ni, nj, k]) # can't remove more than neighbor has
                            end   
                            # update the flux                         
                            net_flux[ni, nj, k] += flux
                            net_flux[i, j, k] -= flux
                        end
                    end
                end
            end
        end
    end

    # apply the flux
    N_new .= N .+ net_flux
    N_new .= max.(N_new, 0.0)  # prevent negative concentrations

    return nothing
end

end

