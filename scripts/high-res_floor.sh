#!/bin/bash
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J high-res_floor
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH -o high-res_floor_%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=robin.dymer@it.uu.se

# Misfit of the true a against the high-resolution data at m = 76, 151 and 251
# (notes/high_res_plan.md, step 0.3). One forward solve per grid.

module load MATLAB  # check the version with `module avail MATLAB`

REPO_DIR="$SLURM_SUBMIT_DIR/.."
RUN_DIR="$REPO_DIR/adjoint_rate_state_inversion"

cd "$RUN_DIR"
srun matlab -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); high_res_misfit_floor"
