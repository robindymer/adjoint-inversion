#!/bin/bash
# Inverse crime (m = 36): H-MALA and H-HMC, 4 chains each (seeds 1-4), Laplace-preconditioned and
# started from an overdispersed Laplace draw, to compare with the existing H-pCN chains at about the
# same cost (~1000 solves): MALA 500 iterations, HMC 100 x 5 leapfrog steps.
# Array task i: sampler = SAMPLERS[(i-1)/4], seed = (i-1)%4 + 1. Submit from scripts/.
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J laplace_cmp
#SBATCH -t 06:00:00
#SBATCH -c 1
#SBATCH --array=1-8
#SBATCH -o laplace_cmp_%A_%a.out

module load MATLAB

REPO_DIR="$SLURM_SUBMIT_DIR/.."
BAYES_DIR="$REPO_DIR/bayes"
SAMPLERS=(mala hmc)
I=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLER=${SAMPLERS[$((I / 4))]}
SEED=$((I % 4 + 1))
if [ "$SAMPLER" = "mala" ]; then
    OPTS="nIter=500; nAdapt=50;"
else
    OPTS="nIter=100; nLeap=5; nAdapt=20; eps0=0.3;"
fi

echo "Task $SLURM_ARRAY_TASK_ID: Laplace $SAMPLER, seed $SEED"
srun matlab -batch "addpath('$BAYES_DIR'); maxNumCompThreads(1); seed=$SEED; start='laplace'; $OPTS run_laplace_$SAMPLER"
