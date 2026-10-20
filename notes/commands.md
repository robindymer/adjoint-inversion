cd adjoint_rate_state_inversion
## Simple solo run
matlab -batch "opt_routine_str='newton-tr'; direct_effect_inverse_crime" > log_newton-tr.txt 2>&1

## Solo run with thread cap
matlab -batch "opt_routine_str='newton-tr'; maxNumCompThreads(3); direct_effect_inverse_crime" > log_newton-tr.txt 2>&1

## Sequential
for s in lbfgs CG newton-tr newton-tr-noprec fminunc-tr newton-cg-ls; do
  matlab -batch "opt_routine_str='$s'; direct_effect_inverse_crime" > log_$s.txt 2>&1
done

## Parallel (with thread cap)
for s in lbfgs CG newton-tr newton-tr-noprec fminunc-tr newton-cg-ls; do
  matlab -batch "opt_routine_str='$s'; maxNumCompThreads(3); direct_effect_inverse_crime" > log_$s.txt 2>&1 &
done
wait