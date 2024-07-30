#!/bin/bash
#SBATCH -p serc
#SBATCH -c 15
#SBATCH -G 6
#SBATCH --output=<input_fnm".out">
#SBATCH --error=<input_fnm".err">
#SBATCH --mail-type=ALL
#SBATCH --mail-type=ALL
#SBATCH --time=42:00:00
#SBATCH --mem=30G

ml load devel julia
julia HydroMech2D_Main.jl
