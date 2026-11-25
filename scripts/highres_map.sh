#!/bin/bash
# High-res Bayesian run, step 1: MAP by Newton-TR, noise inflation (discrepancy principle) and the
# Hessian at the MAP, m = 151, m_p = 26 (bayes/run_highres_map.m).
# Result: bayes/results/highres_map_m151_mp26_eta0.02.mat, needed by highres_hpcn.sh.
# If it hits the time limit or ends with "Hessian not positive definite", resubmit with
#   sbatch --export=ALL,RESUME=true highres_map.sh
# Submit from scripts/ (or use highres_submit.sh, which also queues the chains).
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J highres_map
#SBATCH -t 48:00:00
#SBATCH -c 8
#SBATCH -o highres_map_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robin.dymer@it.uu.se

module load MATLAB  # check the version with `module avail MATLAB`

REPO_DIR="$SLURM_SUBMIT_DIR/.."
BAYES_DIR="$REPO_DIR/bayes"
RESUME="${RESUME:-false}"
MAXITER="${MAXITER:-150}"

srun matlab -batch "addpath('$BAYES_DIR'); maxNumCompThreads($SLURM_CPUS_PER_TASK); m=151; m_p=26; resume=$RESUME; maxIter=$MAXITER; run_highres_map"
