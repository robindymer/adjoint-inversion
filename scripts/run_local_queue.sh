#!/bin/bash
# Runs the lines of a run file locally (the local alternative to high-res_sweep_a.sh). A new run is
# started whenever fewer than <max_parallel> MATLAB batch jobs are running on the machine, counting
# jobs started by other queues, so that several queues can share the memory.
# Each line of the run file is <optimizer> <reg_type> <reg_beta> [reg_eps] [x0_iter], as in
# high-res_sweep_a.sh; lines starting with # are skipped. Logs are written to adjoint_rate_state_inversion/logs/.
# Usage: run_local_queue.sh <driver> <run_file> <max_parallel> <threads> [matlab_assignments]
# e.g.   nohup scripts/run_local_queue.sh high_res_data_inversion scripts/high-res_runs_step23.txt 3 3 "m=76; data_m=251; inv_par='a';" &

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUN_DIR="$SCRIPT_DIR/../adjoint_rate_state_inversion"
DRIVER="$1"; RUN_FILE="$(realpath "$2")"; MAX_PARALLEL="$3"; THREADS="$4"; EXTRA="$5"

running() {
    pgrep -fc "[M]ATLAB -batch maxNumCompThreads"
}

cd "$RUN_DIR"
mkdir -p logs
grep -v '^#' "$RUN_FILE" | grep -v '^\s*$' | while read -r OPT REG_TYPE REG_BETA REG_EPS X0_ITER; do
    NAME="${DRIVER}_${OPT}_${REG_TYPE}_${REG_BETA}"
    EPS_STR=""
    if [ -n "$REG_EPS" ] && [ "$REG_EPS" != "-" ]; then
        EPS_STR="reg_eps=$REG_EPS;"
        NAME="${NAME}_eps${REG_EPS}"
    fi
    X0_STR=""
    if [ -n "$X0_ITER" ]; then
        X0_STR="x0_iter=$X0_ITER;"
        NAME="${NAME}_x0iter${X0_ITER}"
    fi
    while [ "$(running)" -ge "$MAX_PARALLEL" ]; do
        sleep 60
    done
    echo "$(date '+%F %T') starting: $NAME"
    matlab -batch "maxNumCompThreads($THREADS); $EXTRA opt_routine_str='$OPT'; reg_type='$REG_TYPE'; reg_beta=$REG_BETA; $EPS_STR $X0_STR $DRIVER" \
        > "logs/$NAME.txt" 2>&1 < /dev/null &
    sleep 120 # Let MATLAB start before counting the running jobs again
done
wait
echo "$(date '+%F %T') all runs finished: $RUN_FILE"
