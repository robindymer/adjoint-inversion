#!/bin/bash
# Prior-preconditioned pCN, MALA and HMC from prior draws, 4 chains each (seeds 1-4), one core per
# chain, about equal cost (~24 h). Array task i: sampler = SAMPLERS[(i-1)/4], seed = (i-1)%4 + 1.
# Submit from scripts/: sbatch sampler_comparison.sh
# Chains are saved to bayes/results/ as they run, so a job that hits the time limit keeps its samples.
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J sampler_cmp
#SBATCH -t 36:00:00
#SBATCH -c 1
#SBATCH --array=1-12
#SBATCH -o sampler_cmp_%A_%a.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=robin.dymer@it.uu.se

module load MATLAB  # check the version with `module avail MATLAB`

REPO_DIR="$SLURM_SUBMIT_DIR/.."
BAYES_DIR="$REPO_DIR/bayes"

SAMPLERS=(pcn mala hmc)
I=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLER=${SAMPLERS[$((I / 4))]}
SEED=$((I % 4 + 1))

EXTRA=""
if [ "$SAMPLER" = "pcn" ]; then
    EXTRA="nIter=20000;"   # 20000 forward solves, ~4.3 s each
fi                         # mala: 15000 iterations x 1 gradient, hmc: 1500 x 10 gradients, ~5.8 s each

echo "Task $SLURM_ARRAY_TASK_ID: $SAMPLER, seed $SEED"
srun matlab -batch "addpath('$BAYES_DIR'); maxNumCompThreads(1); seed=$SEED; start='prior'; $EXTRA run_$SAMPLER"
