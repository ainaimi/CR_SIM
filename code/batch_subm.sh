#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=30-00:00:00
#SBATCH --job-name=cr_simulation
#SBATCH --mem=120g
#SBATCH --partition=naimi

module purge
module load R

# Simulation Run Script
Rscript --no-save --no-restore --verbose ~/CR_SIM/code/simulation_run.R 200 1200  > ./misc/sim_run.Rout 2>&1