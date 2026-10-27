#!/bin/bash
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J high-res_reg_a
#SBATCH -t 48:00:00
#SBATCH -c 8
#SBATCH -o high-res_reg_a_%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=robin.dymer@it.uu.se

# Regularized inversion for a with high-resolution data, see adjopt.regularization.
# Usage: sbatch high-res_reg_a.sh <optimizer> <reg_type> <reg_beta> [reg_eps]
# e.g.   sbatch high-res_reg_a.sh CG tv 1e-1
#        sbatch high-res_reg_a.sh lbfgs lasso 1e-1 1e-4
# reg_type: tikhonov, tikhonov-grad, tv or lasso

module load MATLAB  # check the version with `module avail MATLAB`

OPT="$1"
REG_TYPE="$2"
REG_BETA="$3"
if [ -n "$4" ]; then
    REG_EPS="reg_eps=$4;"
fi

REPO_DIR="$SLURM_SUBMIT_DIR/.."
RUN_DIR="$REPO_DIR/adjoint_rate_state_inversion"

cd "$RUN_DIR"
srun matlab -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); inv_par='a'; opt_routine_str='$OPT'; reg_type='$REG_TYPE'; reg_beta=$REG_BETA; $REG_EPS high_res_data_inversion"
