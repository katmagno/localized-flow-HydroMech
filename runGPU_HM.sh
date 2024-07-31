#!/bin/bash
#SBATCH -p serc
#SBATCH -c 15
#SBATCH -G 14
#SBATCH --output="np21.out"
#SBATCH --error="np21.err"
#SBATCH --mail-type=ALL
#SBATCH --mail-type=ALL
#SBATCH --time=42:00:00
#SBATCH --mem=40G

ml load devel julia
julia HydroMech2D_Main.jl
