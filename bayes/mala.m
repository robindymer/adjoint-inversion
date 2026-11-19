% Preconditioned Metropolis-adjusted Langevin algorithm (MALA) for the posterior post (see
% posterior.m). Works in whitened coordinates a = mu + R\w of a Gaussian reference N(mu, (R'R)^-1),
% so the preconditioner is the reference covariance (R'R)^-1, by default the prior covariance C:
%
%   w' = w - h/2 grad_w U(w) + sqrt(h) xi,   xi ~ N(0, I),   grad_w U = R'\grad_a U
%
% accepted with the Metropolis-Hastings probability, which includes the asymmetric proposal
% density. One forward and one adjoint solve per iteration. Proposals with a <= 0 or a failed
% solve have U = inf and are rejected.
%
% a0     starting point (column vector at the m_p parameter nodes), must have U(a0) < inf
% nIter  number of iterations
% h      step size (> 0)
% opts.ref         reference Gaussian, struct with mean ref.mu and upper triangular ref.R, where
%                  ref.R'*ref.R is its precision. Default: the prior.
% opts.nAdapt      adapt h during the first nAdapt iterations towards the acceptance rate
%                  opts.target (Robbins-Monro on log h). These iterations are burn-in.
% opts.target      target acceptance rate, default 0.574 (optimal for MALA)
% opts.rs, opts.saveFile, opts.saveEvery, opts.printEvery   as in pcn.m
%
% chain.a          samples as columns (m_p x nIter), chain.U their potential
% chain.alpha      acceptance probabilities, chain.accepted whether the proposal was accepted
% chain.h          step size used in each iteration
% chain.nFwd, chain.nAdj   number of forward and adjoint solves, chain.time wall-clock time [s]
function chain = mala(post, a0, nIter, h, opts)
    default_arg('opts', struct());
    default_field(opts, 'ref', []);
    default_field(opts, 'nAdapt', 0);
    default_field(opts, 'target', 0.574);
    default_field(opts, 'rs', RandStream('mt19937ar', 'Seed', 0));
    default_field(opts, 'saveFile', '');
    default_field(opts, 'saveEvery', 100);
    default_field(opts, 'printEvery', 10);

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
    chain.h = zeros(1, nIter);
    chain.opts = rmfield(opts, 'rs');
    [nFwd0, nAdj0] = solveCounts(post);
    t0 = tic;

    logH = log(h);
    for k = 1:nIter
        xi = randn(rs, m_p, 1);
        wp = w - 0.5*h*g + sqrt(h)*xi;
        ap = ref.mu + R\wp;
        [up, gp] = post.Ug(ap);
        if isfinite(up)
            gp = R'\gp;
            logqFwd = -0.5*(xi'*xi);                             % log q(w' | w)
            logqBwd = -sum((w - wp + 0.5*h*gp).^2)/(2*h);        % log q(w | w')
            alpha = min(1, exp(u - up + logqBwd - logqFwd));
        else
            alpha = 0;
        end
        if rand(rs) < alpha
            w = wp; a = ap; u = up; g = gp;
            chain.accepted(k) = true;
        end
        chain.a(:,k) = a;
        chain.U(k) = u;
        chain.alpha(k) = alpha;
        chain.h(k) = h;

        if k <= opts.nAdapt
            logH = logH + k^(-0.6)*(alpha - opts.target);
            h = exp(logH);
        end

        if mod(k, opts.printEvery) == 0
            fprintf('[%6.0fs] it %5d  U %10.4e  h %8.2e  acc. rate (last %d) %.2f\n', toc(t0), k, ...
                chain.U(k), h, opts.printEvery, mean(chain.accepted(k-opts.printEvery+1:k)));
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
        c.accepted = c.accepted(1:k); c.h = c.h(1:k);
        [c.nFwd, c.nAdj] = solveCounts(post, nFwd0, nAdj0);
        c.time = toc(t0);
        S.chain = c;
        save(opts.saveFile, '-struct', 'S');
    end
end
