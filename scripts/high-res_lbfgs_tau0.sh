#!/bin/bash
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J high-res_lbfgs_tau0
#SBATCH -t 48:00:00
#SBATCH -c 8
#SBATCH -o high-res_lbfgs_tau0_%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=robin.dymer@it.uu.se

module load MATLAB  # check the version with `module avail MATLAB`

REPO_DIR="$SLURM_SUBMIT_DIR/.."
RUN_DIR="$REPO_DIR/adjoint_rate_state_inversion"

cd "$RUN_DIR"
srun matlab -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); inv_par='tau0'; high_res_data_inversion"
