#!/bin/bash
# Run several chains of run_pcn.m or run_hpcn.m in parallel on this machine, one MATLAB process
# (one core) per chain, with seeds 1..n. Logs go to results/logs/.
# Usage: ./run_chains.sh <pcn|hpcn|mala|hmc> <n_chains> [start] [extra MATLAB assignments]
# e.g.   ./run_chains.sh pcn 4 prior
#        ./run_chains.sh hpcn 4 laplace "s=0.3; nIter=2000;"
# start: pcn, mala, hmc: true or prior, hpcn: map or laplace

SAMPLER="$1"
N="$2"
START="${3:-}"
EXTRA="${4:-}"
if [ -n "$START" ]; then
    START="start='$START';"
fi

cd "$(dirname "$0")"
mkdir -p results/logs

# The chains of H-pCN share the cached Hessian: compute it once before starting them
if [ "$SAMPLER" = "hpcn" ]; then
    matlab -batch "hessianOnly=true; run_hpcn" > results/logs/hessian.log 2>&1 || exit 1
fi

for SEED in $(seq 1 "$N"); do
    LOG="results/logs/${SAMPLER}_${3:-default}_seed${SEED}.log"
    matlab -batch "maxNumCompThreads(1); seed=$SEED; $START $EXTRA run_$SAMPLER" > "$LOG" 2>&1 &
    echo "Started $SAMPLER chain seed $SEED, log $LOG"
done
wait
echo "All $N $SAMPLER chains finished"
