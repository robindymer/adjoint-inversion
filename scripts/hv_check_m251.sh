#!/bin/bash
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J hv_check_m251
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH -o hv_check_m251_%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=robin.dymer@it.uu.se

# Finite-difference check of the Hessian-vector product at m = 251 (about 20-30 min).
# Submit from the scripts/ folder. Fails (non-zero exit) if the FD error is too large.

module load MATLAB  # check the version with `module avail MATLAB`

REPO_DIR="$SLURM_SUBMIT_DIR/.."
BENCH_DIR="$REPO_DIR/adjoint_rate_state_inversion/optimizer_benchmark"

srun matlab -batch "addpath('$BENCH_DIR'); hv_check(251, $SLURM_CPUS_PER_TASK)"
