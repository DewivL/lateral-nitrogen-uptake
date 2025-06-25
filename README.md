# Lateral Nitrogen Movement Uptake Effect

This repository contains a voxel-based simulation of nitrogen uptake and lateral movement in soil, designed to compare cropping systems under varying nitrogen conditions.
The model is implemented in Julia and outputs zone-specific and voxel-level nitrogen  for post-analysis in R.
This model was created by Dewi van Lier for an MSc thesis at Wageningen University.

## Repository Structure

- src/
  - Init.jl (Soil and root biomass initialization)
  - Nuptake.jl (Nitrogen uptake calculation (HATS + LATS))
  - Nmovement.jl (Lateral nitrogen diffusion)
  - MyCube.jl (Grid and visualization in Makie)
  - main.jl (Simulation loop)
- simulation_output/
- data_analysis/
- README.md

## How It Works

- Initializes a 3D soil domain with root dry weight distribution.
- Simulates nitrogen uptake using both high-affinity and low-affinity transport.
- Optionally includes lateral diffusion based on concentration gradients.
- Runs multiple treatments (e.g., monocrop vs intercrop, high vs low N) with replicates.
- Outputs:
  - Summary CSVs
  - Time-series N concentration by zone
  - Cumulative uptake by root type
  - Optional voxel-level maps for visualization
