#!/bin/bash
# High-res Bayesian run, step 2b: four H-MALA chains (seeds 1-4) with the Laplace approximation of
# highres_map.sh as preconditioner (bayes/run_highres_mala.m), m = 151, m_p = 26, 500 iterations.
# Chains are saved every 10 iterations. Override e.g. sbatch --export=ALL,NITER=300 highres_mala.sh
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J highres_mala
#SBATCH -t 48:00:00
#SBATCH -c 4
#SBATCH --array=1-4
#SBATCH -o highres_mala_%A_%a.out

module load MATLAB

REPO_DIR="$SLURM_SUBMIT_DIR/.."
BAYES_DIR="$REPO_DIR/bayes"
NITER="${NITER:-500}"

echo "Task $SLURM_ARRAY_TASK_ID: H-MALA high-res, seed $SLURM_ARRAY_TASK_ID"
srun matlab -batch "addpath('$BAYES_DIR'); maxNumCompThreads($SLURM_CPUS_PER_TASK); m=151; m_p=26; nIter=$NITER; seed=$SLURM_ARRAY_TASK_ID; run_highres_mala"
