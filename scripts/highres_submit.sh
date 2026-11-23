#!/bin/bash
# Queue the whole high-res Bayesian run: MAP + Hessian, then H-pCN, H-MALA and H-HMC chains (4 each)
# once it has finished successfully. Run from scripts/: ./highres_submit.sh
# With an existing MAP (e.g. after a RESUME run): ./highres_submit.sh --skip-map
cd "$(dirname "$0")"
DEP=""
if [ "$1" != "--skip-map" ]; then
    MAP_JOB=$(sbatch --parsable highres_map.sh)
    echo "MAP + Hessian: job $MAP_JOB"
    DEP="--dependency=afterok:$MAP_JOB"
fi
for s in hpcn mala hmc; do
    JOB=$(sbatch --parsable $DEP highres_$s.sh)
    echo "$s chains: job $JOB"
done
