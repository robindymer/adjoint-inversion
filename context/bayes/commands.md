## Start 4 pCN chains, different initial guesses (drawn from prior)
nohup ./bayes/run_chains.sh pcn 4 prior > /dev/null 2>&1 &

## Start 4 H-pCN chains, different initial guesses (drawn from laplace with doubled spread)
nohup ./bayes/run_chains.sh hpcn 4 laplace "s=1;" > /dev/null 2>&1 &