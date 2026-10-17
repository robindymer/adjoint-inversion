% Trust-region Newton-CG (Steihaug-Toint) using exact Hessian-vector products.
%
% Each outer iteration solves the trust-region subproblem
%   min_s g'*s + 0.5*s'*H*s,  ||s||_M <= Delta
% with preconditioned CG, stopping at negative curvature or at the trust-region boundary.
% The preconditioner M^{-1} is an L-BFGS matrix built from the CG pairs (d, H*d) of the
% previous accepted iteration (Morales & Nocedal, 2000).
%
% fg      -- handle, [f, g] = fg(x). Called once per trial point.
% Hv      -- handle, HV = Hv(x, V), Hessian at x applied to the columns of V.
%            Called only at accepted points, after fg has been called there.
% x       -- initial guess (column vector)
% lb      -- lower bound, x > lb is kept by a fraction-to-boundary rule. [] for no bound.
% opts    -- struct with optional fields
%            maxIter    maximum number of outer iterations (default 100)
%            maxCG      maximum CG iterations per outer iteration (default 20)
%            etaMax     largest CG forcing term (default 0.1)
%            gtolRel    stop when ||g|| <= gtolRel*max_k ||g_k|| (default 1e-8)
%            Delta0     initial trust-region radius (default 0.1*norm(x))
%            DeltaMax   largest trust-region radius (default 10*norm(x))
%            maxStep    largest Euclidean step, guards against huge steps along
%                       near-null directions of the Hessian (default 0.25*norm(x))
%            precond    use the L-BFGS preconditioner (default true)
%            memL       number of stored L-BFGS pairs (default 20)
%            ftb        fraction to boundary (default 0.9)
%            OutputFcn  handle, stop = OutputFcn(x, f, g, info), called at accepted points
%
% flag    -- 0: gradient tolerance met, 1: maxIter reached, 2: trust region collapsed,
%            3: stopped by OutputFcn
function [x, f, flag, output] = trustRegionNewtonCG(fg, Hv, x, lb, opts)
    default_arg('lb', []);
    default_arg('opts', struct);
    opts = defaultField(opts, 'maxIter', 100);
    opts = defaultField(opts, 'maxCG', 20);
    opts = defaultField(opts, 'etaMax', 0.1);
    opts = defaultField(opts, 'gtolRel', 1e-8);
    opts = defaultField(opts, 'Delta0', 0.1*norm(x));
    opts = defaultField(opts, 'DeltaMax', 10*norm(x));
    opts = defaultField(opts, 'maxStep', 0.25*norm(x));
    opts = defaultField(opts, 'precond', true);
    opts = defaultField(opts, 'memL', 20);
    opts = defaultField(opts, 'ftb', 0.9);
    opts = defaultField(opts, 'OutputFcn', []);

    Delta = opts.Delta0;
    [f, g] = fg(x);
    gmax = norm(g);
    S = zeros(numel(x), 0); % L-BFGS pairs for the preconditioner
    Y = zeros(numel(x), 0);
    memo = emptyMemo(x); % Hessian-vector products at x, reused when a step is rejected
    output = struct('iterations', 0, 'funcCount', 1, 'hessVecCount', 0, 'cgIterations', 0);
    flag = 1;
    if callOutput(opts, x, f, g, 'init')
        flag = 3; return
    end

    for it = 1:opts.maxIter
        gmax = max(gmax, norm(g));
        if norm(g) <= opts.gtolRel*gmax
            flag = 0; break
        end
        if Delta < 1e-12*norm(x)
            flag = 2; break
        end

        % Solve the trust-region subproblem
        if opts.precond && ~isempty(S)
            Minv = @(r) lbfgsTwoLoop(r, S, Y);
        else
            Minv = @(r) r;
        end
        eta = min(opts.etaMax, sqrt(norm(g)/gmax));
        nHv = size(memo.D, 2);
        [s, Hs, cg, memo] = steihaugCG(@(v) Hv(x, v), memo, g, Minv, Delta, eta*norm(g), opts.maxCG);
        output.hessVecCount = output.hessVecCount + size(memo.D, 2) - nHv;
        output.cgIterations = output.cgIterations + cg.iter;

        % Shorten the step to keep x > lb and ||s|| <= maxStep
        th = min(1, opts.maxStep/norm(s));
        if ~isempty(lb)
            dec = s < 0;
            th = min([th; opts.ftb*(x(dec) - lb(dec))./(-s(dec))]);
        end
        pred = -(th*(g'*s) + 0.5*th^2*(s'*Hs));

        % Evaluate the trial point and update the trust region
        xt = x + th*s;
        [ft, gt] = fg(xt);
        output.funcCount = output.funcCount + 1;
        rho = (f - ft)/pred;
        if ~isfinite(rho)
            rho = -inf;
        end
        if rho < 0.25
            Delta = 0.25*min(Delta, th*cg.sMnorm);
        elseif rho > 0.75 && cg.boundary && th == 1
            Delta = min(4*Delta, opts.DeltaMax);
        end

        info = sprintf('CG %2d %-8s rho %7.3f th %4.2f Delta %8.2e', cg.iter, cg.flag, rho, th, Delta);
        if rho > 1e-4 % Accept
            if opts.precond
                S = [S, memo.D(:, memo.positive)]; %#ok<AGROW>
                Y = [Y, memo.HD(:, memo.positive)]; %#ok<AGROW>
                S = S(:, max(1, end-opts.memL+1):end);
                Y = Y(:, max(1, end-opts.memL+1):end);
            end
            x = xt; f = ft; g = gt;
            memo = emptyMemo(x);
            output.iterations = it;
            if callOutput(opts, x, f, g, info)
                flag = 3; break
            end
        else
            fprintf('rejected: %s, f_trial = %e\n', info, ft);
        end
    end
end

% Preconditioned Steihaug CG for min g'*s + 0.5*s'*H*s subject to ||s||_M <= Delta.
% Hessian-vector products are looked up in memo first, so that a re-solve with a smaller
% Delta after a rejected step replays the same CG iterations without new products.
function [s, Hs, info, memo] = steihaugCG(Hfun, memo, g, Minv, Delta, tol, maxit)
    n = numel(g);
    s = zeros(n, 1);
    Hs = zeros(n, 1);
    r = g;
    z = Minv(r);
    d = -z;
    rz = r'*z;
    % ||s||_M^2, <s, M d> and ||d||_M^2 updated by recurrences, M is never formed
    sMs = 0;
    sMd = 0;
    dMd = rz;
    info.flag = 'maxit';
    info.boundary = false;
    info.iter = maxit;
    for j = 1:maxit
        [Hd, memo] = memoHv(Hfun, memo, d);
        kappa = d'*Hd;
        if kappa <= 0 % Negative curvature, follow d to the boundary
            tau = toBoundary(sMs, sMd, dMd, Delta);
            s = s + tau*d;
            Hs = Hs + tau*Hd;
            info.flag = 'negcurv'; info.boundary = true; info.iter = j; info.sMnorm = Delta;
            return
        end
        alpha = rz/kappa;
        sMsNew = sMs + 2*alpha*sMd + alpha^2*dMd;
        if sMsNew >= Delta^2 % Step leaves the trust region
            tau = toBoundary(sMs, sMd, dMd, Delta);
            s = s + tau*d;
            Hs = Hs + tau*Hd;
            info.flag = 'boundary'; info.boundary = true; info.iter = j; info.sMnorm = Delta;
            return
        end
        s = s + alpha*d;
        Hs = Hs + alpha*Hd;
        sMs = sMsNew;
        r = r + alpha*Hd;
        if norm(r) <= tol
            info.flag = 'conv'; info.iter = j;
            break
        end
        z = Minv(r);
        rzNew = r'*z;
        beta = rzNew/rz;
        sMd = beta*(sMd + alpha*dMd);
        dMd = rzNew + beta^2*dMd;
        d = -z + beta*d;
        rz = rzNew;
    end
    info.sMnorm = sqrt(sMs);
end

function tau = toBoundary(sMs, sMd, dMd, Delta)
    tau = (-sMd + sqrt(sMd^2 + dMd*(Delta^2 - sMs)))/dMd;
end

function memo = emptyMemo(x)
    memo.D = zeros(numel(x), 0);
    memo.HD = zeros(numel(x), 0);
    memo.positive = false(1, 0);
end

function [Hd, memo] = memoHv(Hfun, memo, d)
    k = find(all(memo.D == d, 1), 1);
    if isempty(k)
        Hd = Hfun(d);
        memo.D(:, end+1) = d;
        memo.HD(:, end+1) = Hd;
        memo.positive(end+1) = d'*Hd > 0;
    else
        Hd = memo.HD(:, k);
    end
end

% L-BFGS two-loop recursion, returns H*q where H approximates the inverse Hessian.
function r = lbfgsTwoLoop(q, S, Y)
    L = size(S, 2);
    rho = 1./sum(S.*Y, 1);
    a = zeros(L, 1);
    for i = L:-1:1
        a(i) = rho(i)*(S(:, i)'*q);
        q = q - a(i)*Y(:, i);
    end
    r = (S(:, L)'*Y(:, L))/(Y(:, L)'*Y(:, L))*q;
    for i = 1:L
        b = rho(i)*(Y(:, i)'*r);
        r = r + (a(i) - b)*S(:, i);
    end
end

function stop = callOutput(opts, x, f, g, info)
    stop = false;
    fprintf('f = %12.6e  |g| = %10.3e  %s\n', f, norm(g), info);
    if ~isempty(opts.OutputFcn)
        stop = opts.OutputFcn(x, f, g, info);
    end
end
