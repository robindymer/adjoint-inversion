% Hamiltonian Monte Carlo (HMC) for the posterior post (see posterior.m), MATLAB port of hmc() in the
% course code (course_code/samplers.py). Works in whitened coordinates a = mu + R\w of a Gaussian
% reference N(mu, (R'R)^-1) with identity mass matrix in w, i.e. mass matrix R'R in a. By default
% the reference is the prior, so the mass matrix is the prior precision C^-1.
%
% Each iteration draws a momentum r ~ N(0, I), takes nLeap leapfrog steps of size eps (jittered
% uniformly in [(1 - jitter) eps, (1 + jitter) eps], Neal 2011) and accepts with probability
% min(1, exp(H0 - H1)), H = U + |r|^2/2. nLeap forward and adjoint solves per iteration. A trajectory
% that reaches a <= 0 or a failed solve is stopped and rejected.
%
% a0     starting point (column vector at the m_p parameter nodes), must have U(a0) < inf
% nIter  number of iterations
% eps    leapfrog step size, nLeap number of leapfrog steps
% opts.ref         reference Gaussian (mass matrix), as in mala.m. Default: the prior.
% opts.nAdapt      adapt eps during the first nAdapt iterations towards the acceptance rate
%                  opts.target (Robbins-Monro on log eps). These iterations are burn-in.
% opts.target      target acceptance rate, default 0.65
% opts.jitter      relative jitter of eps, default 0.2
% opts.rs, opts.saveFile, opts.saveEvery, opts.printEvery   as in pcn.m
%
% chain.a          samples as columns (m_p x nIter), chain.U their potential
% chain.alpha      acceptance probabilities, chain.accepted whether the proposal was accepted
% chain.dH         energy error H1 - H0 (inf for a stopped trajectory)
% chain.eps        (unjittered) step size used in each iteration, chain.nLeap
% chain.nFwd, chain.nAdj   number of forward and adjoint solves, chain.time wall-clock time [s]
function chain = hmc(post, a0, nIter, eps, nLeap, opts)
    default_arg('opts', struct());
    default_field(opts, 'ref', []);
    default_field(opts, 'nAdapt', 0);
    default_field(opts, 'target', 0.65);
    default_field(opts, 'jitter', 0.2);
    default_field(opts, 'rs', RandStream('mt19937ar', 'Seed', 0));
    default_field(opts, 'saveFile', '');
    default_field(opts, 'saveEvery', 10);
    default_field(opts, 'printEvery', 1);

    ref = opts.ref;
    if isempty(ref)
        ref = post.prior;
    end
    R = ref.R;
    rs = opts.rs;
    m_p = numel(a0);
    w = R*(a0 - ref.mu);
    [u, g] = post.Ug(a0);
    if ~isfinite(u)
        error('U(a0) is not finite');
    end
    g = R'\g;   % gradient with respect to w
    a = a0;

    chain.a = zeros(m_p, nIter);
    chain.U = zeros(1, nIter);
    chain.alpha = zeros(1, nIter);
    chain.accepted = false(1, nIter);
    chain.dH = zeros(1, nIter);
    chain.eps = zeros(1, nIter);
    chain.nLeap = nLeap;
    chain.opts = rmfield(opts, 'rs');
    [nFwd0, nAdj0] = solveCounts(post);
    t0 = tic;

    logEps = log(eps);
    for k = 1:nIter
        r = randn(rs, m_p, 1);
        H0 = u + 0.5*(r'*r);
        hk = eps*(1 + opts.jitter*(2*rand(rs) - 1));

        % Leapfrog
        wn = w; gn = g;
        rn = r - 0.5*hk*gn;
        ok = true;
        for j = 1:nLeap
            wn = wn + hk*rn;
            an = ref.mu + R\wn;
            [un, gan] = post.Ug(an);
            if ~isfinite(un)
                ok = false;
                break
            end
            gn = R'\gan;
            if j < nLeap
                rn = rn - hk*gn;
            end
        end

        if ok
            rn = rn - 0.5*hk*gn;
            dH = un + 0.5*(rn'*rn) - H0;
            alpha = min(1, exp(-dH));
        else
            dH = inf;
            alpha = 0;
        end
        if rand(rs) < alpha
            w = wn; a = an; u = un; g = gn;
            chain.accepted(k) = true;
        end
        chain.a(:,k) = a;
        chain.U(k) = u;
        chain.alpha(k) = alpha;
        chain.dH(k) = dH;
        chain.eps(k) = eps;

        if k <= opts.nAdapt
            logEps = logEps + k^(-0.6)*(alpha - opts.target);
            eps = exp(logEps);
        end

        if mod(k, opts.printEvery) == 0
            fprintf('[%6.0fs] it %5d  U %10.4e  eps %8.2e  dH %9.2e  acc. rate (last %d) %.2f\n', ...
                toc(t0), k, chain.U(k), eps, dH, opts.printEvery, ...
                mean(chain.accepted(k-opts.printEvery+1:k)));
        end
        if ~isempty(opts.saveFile) && (mod(k, opts.saveEvery) == 0 || k == nIter)
            saveChain(k);
        end
    end
    [chain.nFwd, chain.nAdj] = solveCounts(post, nFwd0, nAdj0);
    chain.time = toc(t0);

    function saveChain(k)
        c = chain;
        c.a = c.a(:,1:k); c.U = c.U(1:k); c.alpha = c.alpha(1:k);
        c.accepted = c.accepted(1:k); c.dH = c.dH(1:k); c.eps = c.eps(1:k);
        [c.nFwd, c.nAdj] = solveCounts(post, nFwd0, nAdj0);
        c.time = toc(t0);
        S.chain = c;
        save(opts.saveFile, '-struct', 'S');
    end
end
