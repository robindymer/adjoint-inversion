% Preconditioned Crank-Nicolson (pCN) sampler for the posterior post (see posterior.m), whose
% prior is Gaussian. Works in whitened coordinates a = mu + R\w of a Gaussian reference N(mu, (R'R)^-1),
% which is N(0, I) in w:
%
%   w' = sqrt(1 - s^2) w + s xi,   xi ~ N(0, I),   accepted with prob. min(1, exp(phi(a) - phi(a')))
%
% The proposal keeps the reference, so only the rest of the posterior, phi(a) = U(a) - |w|^2/2,
% enters the acceptance (script, Sec. 4.2.2 and the end of Sec. 4.5 on elliptical slice sampling).
%   - Reference = prior (default): plain pCN, phi = Phi is the likelihood.
%   - Reference = Laplace approximation N(a_MAP, H^-1), with H the Hessian of U (opts.ref):
%     H-pCN. The steps then have the shape of the posterior, and phi is nearly constant if the
%     posterior is nearly Gaussian, so that large s are accepted.
% One forward solve per iteration. Proposals with a <= 0 or a failed solve have phi = inf and are
% rejected.
%
% a0     starting point (column vector at the m_p parameter nodes), must have U(a0) < inf
% nIter  number of iterations
% s      step size in (0, 1]
% opts.ref         reference Gaussian, struct with mean ref.mu and upper triangular ref.R, where
%                  ref.R'*ref.R is its precision, e.g. ref.R = chol(H). Default: the prior.
% opts.nAdapt      adapt s during the first nAdapt iterations towards the acceptance rate
%                  opts.target (Robbins-Monro on log s). These iterations are not a valid chain
%                  and must be discarded as burn-in. Default 0, no adaptation.
% opts.target      target acceptance rate, default 0.25
% opts.rs          RandStream, default mt19937ar with seed 0
% opts.saveFile    if non-empty, the chain so far is saved there every opts.saveEvery iterations
%                  and at the end, so a killed job keeps its samples
% opts.printEvery  print progress every printEvery iterations, default 10
%
% chain.a          samples as columns (m_p x nIter), chain.U their potential and chain.phi their
%                  phi (= Phi for the prior reference)
% chain.alpha      acceptance probabilities, chain.accepted whether the proposal was accepted
% chain.s          step size used in each iteration
% chain.nFwd       number of forward solves, chain.time wall-clock time [s]
function chain = pcn(post, a0, nIter, s, opts)
    default_arg('opts', struct());
    default_field(opts, 'ref', []);
    default_field(opts, 'nAdapt', 0);
    default_field(opts, 'target', 0.25);
    default_field(opts, 'rs', RandStream('mt19937ar', 'Seed', 0));
    default_field(opts, 'saveFile', '');
    default_field(opts, 'saveEvery', 100);
    default_field(opts, 'printEvery', 10);

    ref = opts.ref;
    if isempty(ref)
        ref = post.prior;
        phiFun = @(a, w) post.Phi(a);   % U(a) - |w|^2/2 = Phi(a) for the prior
    else
        phiFun = @(a, w) post.U(a) - 0.5*(w'*w);
    end
    rs = opts.rs;
    m_p = numel(a0);
    w = ref.R*(a0 - ref.mu);
    a = a0;
    phi = phiFun(a, w);
    if ~isfinite(phi)
        error('U(a0) is not finite');
    end

    chain.a = zeros(m_p, nIter);
    chain.U = zeros(1, nIter);
    chain.phi = zeros(1, nIter);
    chain.alpha = zeros(1, nIter);
    chain.accepted = false(1, nIter);
    chain.s = zeros(1, nIter);
    chain.opts = rmfield(opts, 'rs');
    nFwd0 = post.prob.nFwd;
    t0 = tic;

    logS = log(s);
    for k = 1:nIter
        wp = sqrt(1 - s^2)*w + s*randn(rs, m_p, 1);
        ap = ref.mu + ref.R\wp;
        phip = phiFun(ap, wp);
        alpha = min(1, exp(phi - phip));   % 0 if phip = inf
        if rand(rs) < alpha
            w = wp; a = ap; phi = phip;
            chain.accepted(k) = true;
        end
        chain.a(:,k) = a;
        chain.U(k) = phi + 0.5*(w'*w);
        chain.phi(k) = phi;
        chain.alpha(k) = alpha;
        chain.s(k) = s;

        % Idea: if a above 25 percent target, increase s, otherwise decrease s (The Robbins-Monro update)
        if k <= opts.nAdapt
            logS = logS + k^(-0.6)*(alpha - opts.target);
            s = min(exp(logS), 1);
        end

        if mod(k, opts.printEvery) == 0
            fprintf('[%6.0fs] it %5d  U %10.4e  s %8.2e  acc. rate (last %d) %.2f\n', toc(t0), k, ...
                chain.U(k), s, opts.printEvery, mean(chain.accepted(k-opts.printEvery+1:k)));
        end
        if ~isempty(opts.saveFile) && (mod(k, opts.saveEvery) == 0 || k == nIter)
            saveChain(k);
        end
    end
    chain.nFwd = post.prob.nFwd - nFwd0;
    chain.time = toc(t0);

    % Saves the first k iterations under the variable name chain
    function saveChain(k)
        c = chain;
        c.a = c.a(:,1:k); c.U = c.U(1:k); c.phi = c.phi(1:k); c.alpha = c.alpha(1:k);
        c.accepted = c.accepted(1:k); c.s = c.s(1:k);
        c.nFwd = post.prob.nFwd - nFwd0;
        c.time = toc(t0);
        S.chain = c;
        save(opts.saveFile, '-struct', 'S');
    end
end
