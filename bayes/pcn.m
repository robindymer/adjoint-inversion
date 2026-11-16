% Preconditioned Crank-Nicolson (pCN) sampler for the posterior post (see posterior.m), whose
% prior is Gaussian. Works in whitened coordinates a = mu + R\w, where the prior is N(0, I):
%
%   w' = sqrt(1 - s^2) w + s xi,   xi ~ N(0, I),   accepted with prob. min(1, exp(Phi(a) - Phi(a')))
%
% The proposal keeps the prior, so only the likelihood Phi enters the acceptance (script, Sec. 4.2.2).
% One forward solve per iteration. Proposals with a <= 0 or a failed solve have Phi = inf and are
% rejected.
%
% a0     starting point (column vector at the m_p parameter nodes), must have Phi(a0) < inf
% nIter  number of iterations
% s      step size in (0, 1]
% opts.nAdapt      adapt s during the first nAdapt iterations towards the acceptance rate
%                  opts.target (Robbins-Monro on log s). These iterations are not a valid chain
%                  and must be discarded as burn-in. Default 0, no adaptation.
% opts.target      target acceptance rate, default 0.25
% opts.rs          RandStream, default mt19937ar with seed 0
% opts.saveFile    if non-empty, the chain so far is saved there every opts.saveEvery iterations
%                  and at the end, so a killed job keeps its samples
% opts.printEvery  print progress every printEvery iterations, default 10
%
% chain.a          samples as columns (m_p x nIter), chain.phi their Phi
% chain.alpha      acceptance probabilities, chain.accepted whether the proposal was accepted
% chain.s          step size used in each iteration
% chain.nFwd       number of forward solves, chain.time wall-clock time [s]
function chain = pcn(post, a0, nIter, s, opts)
    default_arg('opts', struct());
    default_field(opts, 'nAdapt', 0);
    default_field(opts, 'target', 0.25);
    default_field(opts, 'rs', RandStream('mt19937ar', 'Seed', 0));
    default_field(opts, 'saveFile', '');
    default_field(opts, 'saveEvery', 100);
    default_field(opts, 'printEvery', 10);

    prior = post.prior;
    rs = opts.rs;
    m_p = numel(a0);
    w = prior.R*(a0 - prior.mu);
    a = a0;
    phi = post.Phi(a);
    if ~isfinite(phi)
        error('Phi(a0) is not finite');
    end

    chain.a = zeros(m_p, nIter);
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
        ap = prior.fromWhitened(wp);
        phip = post.Phi(ap);
        alpha = min(1, exp(phi - phip));   % 0 if phip = inf
        if rand(rs) < alpha
            w = wp; a = ap; phi = phip;
            chain.accepted(k) = true;
        end
        chain.a(:,k) = a;
        chain.phi(k) = phi;
        chain.alpha(k) = alpha;
        chain.s(k) = s;

        % Idea: if a above 25 percent target, increase s, otherwise decrease s (The Robbins-Monro update)
        if k <= opts.nAdapt
            logS = logS + k^(-0.6)*(alpha - opts.target);
            s = min(exp(logS), 1);
        end

        if mod(k, opts.printEvery) == 0
            fprintf('[%6.0fs] it %5d  Phi %10.4e  s %8.2e  acc. rate (last %d) %.2f\n', toc(t0), k, ...
                phi, s, opts.printEvery, mean(chain.accepted(k-opts.printEvery+1:k)));
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
        c.a = c.a(:,1:k); c.phi = c.phi(1:k); c.alpha = c.alpha(1:k);
        c.accepted = c.accepted(1:k); c.s = c.s(1:k);
        c.nFwd = post.prob.nFwd - nFwd0;
        c.time = toc(t0);
        S.chain = c;
        save(opts.saveFile, '-struct', 'S');
    end
end
