#!/bin/bash
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J newton_tr_m251
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH -o newton_tr_m251_%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=robin.dymer@it.uu.se

# newton-tr inversion for a at m = 251 (estimated 5-7 h with 3 threads).
# Submit from the scripts/ folder. History is saved to
# optimizer_benchmark/results/res_newton_tr_repo_m251.mat when the run ends.

module load MATLAB  # check the version with `module avail MATLAB`

REPO_DIR="$SLURM_SUBMIT_DIR/.."
BENCH_DIR="$REPO_DIR/adjoint_rate_state_inversion/optimizer_benchmark"
# Stop 4 h before the Slurm limit so the results are saved (one iteration can take a while at m = 251)
MAX_TIME=$((20*3600))

srun matlab -batch "addpath('$BENCH_DIR'); run_method('newton_tr_repo', $MAX_TIME, $SLURM_CPUS_PER_TASK, 251)"
