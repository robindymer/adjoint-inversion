""" Diagnostic and verification, based on 4.6 in Glaubitz course litterature.
Want to answer
1) Has the chain forgotten its start (diagnostic)
2) Value of ESS, i.e. how many independent samples is it worth (diagnostic)
3) Does it sample the right possterior (verification)
"A sampler never reports its own failure, and every result needs several chains,
R hat, effective sample sizes, and a verification against a known answer" - chapter 4 summary, Glaubitz course litterature.
"""

# First, plot QOI (a(x)) against iterations for each chain
def load_chains(pattern): # e.g. pattern is hpcn_m36_eta0.02_s1_laplace_seed*.mat
    files = sorted(glob.glob(pattern))