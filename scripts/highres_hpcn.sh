#!/bin/bash
# High-res Bayesian run, step 2: four H-pCN chains (seeds 1-4) with the Laplace approximation of
# highres_map.sh as reference (bayes/run_highres_hpcn.m), m = 151, m_p = 26, s = 1, 1000 samples.
# Chains are saved every 10 iterations, so a job that hits the time limit keeps its samples.
# Check the acceptance rate and std(phi) in the log: if std(phi) >> 1, rerun with a smaller s,
#   sbatch --export=ALL,S=0.5 highres_hpcn.sh
# Diagnostics: cd bayes/python; python diagnostic.py "../results/hpcn_highres_m151_mp26_eta0.02_s1_seed*.mat"
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J highres_hpcn
#SBATCH -t 48:00:00
#SBATCH -c 4
#SBATCH --array=1-4
#SBATCH -o highres_hpcn_%A_%a.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=robin.dymer@it.uu.se

module load MATLAB  # check the version with `module avail MATLAB`

REPO_DIR="$SLURM_SUBMIT_DIR/.."
BAYES_DIR="$REPO_DIR/bayes"
S="${S:-1}"
NITER="${NITER:-1000}"

echo "Task $SLURM_ARRAY_TASK_ID: H-pCN high-res, seed $SLURM_ARRAY_TASK_ID, s = $S"
srun matlab -batch "addpath('$BAYES_DIR'); maxNumCompThreads($SLURM_CPUS_PER_TASK); m=151; m_p=26; s=$S; nIter=$NITER; seed=$SLURM_ARRAY_TASK_ID; run_highres_hpcn"
