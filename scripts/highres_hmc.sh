#!/bin/bash
# High-res Bayesian run, step 2c: four H-HMC chains (seeds 1-4) with mass matrix H from
# highres_map.sh (bayes/run_highres_hmc.m), m = 151, m_p = 26, 100 iterations of 5 leapfrog steps.
# Chains are saved every 5 iterations. Override e.g. sbatch --export=ALL,NITER=60,NLEAP=5 highres_hmc.sh
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J highres_hmc
#SBATCH -t 48:00:00
#SBATCH -c 4
#SBATCH --array=1-4
#SBATCH -o highres_hmc_%A_%a.out

module load MATLAB

REPO_DIR="$SLURM_SUBMIT_DIR/.."
BAYES_DIR="$REPO_DIR/bayes"
NITER="${NITER:-100}"
NLEAP="${NLEAP:-5}"

echo "Task $SLURM_ARRAY_TASK_ID: H-HMC high-res, seed $SLURM_ARRAY_TASK_ID, nLeap = $NLEAP"
srun matlab -batch "addpath('$BAYES_DIR'); maxNumCompThreads($SLURM_CPUS_PER_TASK); m=151; m_p=26; nIter=$NITER; nLeap=$NLEAP; seed=$SLURM_ARRAY_TASK_ID; run_highres_hmc"
