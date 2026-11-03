#!/bin/bash
#SBATCH -A uppmax2026-1-45
#SBATCH -p pelle
#SBATCH -J high-res_sweep_a
#SBATCH -t 24:00:00
#SBATCH -c 4
#SBATCH -o high-res_sweep_a_%A_%a.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=robin.dymer@it.uu.se

# Array of inversions for a with high-resolution data at m = 76 (notes/high_res_plan.md, steps 2-4).
# Each array task runs one line of the run file: <optimizer> <reg_type> <reg_beta> [reg_eps].
# Lines starting with # are skipped. All tasks use the same number of cores, so that the wall
# times can be compared.
# The data are the m = 1001 high-res data, or the data generated with data_m (second argument).
# The optional third argument sets the parameter grid size m_p (default 26); pass "" as data_m to
# keep the m = 1001 data.
# Usage: sbatch --array=1-$(grep -cv '^#' high-res_runs_step23.txt) high-res_sweep_a.sh high-res_runs_step23.txt [data_m] [m_p]

module load MATLAB  # check the version with `module avail MATLAB`

RUN_FILE="$SLURM_SUBMIT_DIR/$1"
LINE=$(grep -v '^#' "$RUN_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")
read -r OPT REG_TYPE REG_BETA REG_EPS <<< "$LINE"
if [ -n "$REG_EPS" ]; then
    REG_EPS="reg_eps=$REG_EPS;"
fi
if [ -n "$2" ]; then
    DATA_M="data_m=$2;"
fi
if [ -n "$3" ]; then
    M_P="m_p_in=$3;"
fi
echo "Task $SLURM_ARRAY_TASK_ID: $LINE $DATA_M $M_P"

REPO_DIR="$SLURM_SUBMIT_DIR/.."
RUN_DIR="$REPO_DIR/adjoint_rate_state_inversion"

cd "$RUN_DIR"
srun matlab -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); m=76; $DATA_M inv_par='a'; $M_P opt_routine_str='$OPT'; reg_type='$REG_TYPE'; reg_beta=$REG_BETA; $REG_EPS high_res_data_inversion"
