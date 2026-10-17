% Trust-region Newton-CG (Steihaug-Toint) with optional L-BFGS preconditioner and
% optional switch to L-BFGS once the Hessian looks positive definite.
%
% P     -- InvProblem (f, fg, Hv, log, timeUp)
% opts  -- maxIter, Delta0, DeltaMax, maxCG, gtol, precond (true/false), memL,
%          switchToLBFGS (true/false), nConvexToSwitch, ftb (fraction to boundary for x > 0)
function [x, P] = tr_ncg(P, x, opts)
    opts = defaultField(opts, 'maxIter', 100);
    opts = defaultField(opts, 'Delta0', 0.1*norm(x));
    opts = defaultField(opts, 'DeltaMax', 10*norm(x));
    opts = defaultField(opts, 'maxCG', 20);
    opts = defaultField(opts, 'gtol', 1e-10);
    opts = defaultField(opts, 'precond', false);
    opts = defaultField(opts, 'memL', 20);
    opts = defaultField(opts, 'switchToLBFGS', false);
    opts = defaultField(opts, 'nConvexToSwitch', 3);
    opts = defaultField(opts, 'ftb', 0.9);
    opts = defaultField(opts, 'etaMax', 0.1);

    Delta = opts.Delta0;
    [f, g] = P.fg(x);
    gmax = norm(g);
    P.log(x, f, g, 'init');
    S = zeros(numel(x), 0); Y = zeros(numel(x), 0); % curvature pairs for the preconditioner
    nConvex = 0;
    memoX = []; memoD = []; memoHD = [];

    function Hd = Hmemo(d)
        k = find(all(memoD == d, 1), 1);
        if ~isempty(k)
            Hd = memoHD(:,k);
        else
            Hd = P.Hv(x, d);
            memoD(:,end+1) = d; memoHD(:,end+1) = Hd;
        end
    end

    for it = 1:opts.maxIter
        if P.timeUp() || norm(g) < opts.gtol
            break
        end
        if opts.precond && ~isempty(S)
            Minv = @(r) lbfgs_twoloop(r, S, Y);
        else
            Minv = @(r) r;
        end
        gmax = max(gmax, norm(g));
        eta = min(opts.etaMax, sqrt(norm(g)/gmax)); % forcing term, tightens as |g| drops below its peak
        if ~isequal(x, memoX) % Hv memo, so a re-solve after a rejected step replays CG for free
            memoX = x; memoD = zeros(numel(x),0); memoHD = zeros(numel(x),0);
        end
        [s, Hs, cgInfo, Sp, Yp] = steihaug(@Hmemo, g, Minv, Delta, eta*norm(g), opts.maxCG);

        % Fraction to boundary for positivity, x + th*s >= (1-ftb)*x
        neg = s < 0;
        th = min([1; opts.ftb*x(neg)./(-s(neg))]);
        gs = g'*s; sHs = s'*Hs;
        pred = -(th*gs + 0.5*th^2*sHs);
        xt = x + th*s;
        ft = P.f(xt);
        rho = (f - ft)/pred;

        sMnorm = cgInfo.sMnorm*th;
        if rho < 0.25
            Delta = 0.25*sMnorm;
        elseif rho > 0.9 && cgInfo.boundary && th == 1
            Delta = min(4*Delta, opts.DeltaMax);
        elseif rho > 0.75 && cgInfo.boundary && th == 1
            Delta = min(2*Delta, opts.DeltaMax);
        end

        info = sprintf('CG %2d %-8s rho %6.3f th %4.2f Delta %8.2e', cgInfo.iter, cgInfo.flag, rho, th, Delta);
        if rho > 1e-4
            if opts.precond % Morales-Nocedal: precondition next solve with the CG pairs from this x
                S = [S, Sp]; Y = [Y, Yp];
                S = S(:, max(1,end-opts.memL+1):end); Y = Y(:, max(1,end-opts.memL+1):end);
            end
            x = xt;
            [f, g] = P.fg(x);
            P.log(x, f, g, info);
        else
            fprintf('    rejected: %s  f_trial %10.4e\n', info, ft);
        end

        % Interior Newton step with no negative curvature and a good model -> locally convex
        if (strcmp(cgInfo.flag, 'conv') || strcmp(cgInfo.flag, 'maxit')) && cgInfo.iter >= 3
            if rho > 0.75, nConvex = nConvex + 1; else, nConvex = 0; end
        else
            nConvex = 0;
        end
        if opts.switchToLBFGS && nConvex >= opts.nConvexToSwitch
            fprintf('Switching to L-BFGS at iteration %d, seeded with %d CG pairs\n', it, size(S,2));
            lopts = opts; lopts.maxIter = opts.maxIter - it;
            [x, P] = lbfgs_ls(P, x, lopts, S, Y);
            return
        end
    end
end

% Preconditioned Steihaug CG for min g's + 0.5 s'Hs, ||s||_M <= Delta, with M^{-1} = Minv.
function [s, Hs, info, Sp, Yp] = steihaug(Hfun, g, Minv, Delta, tol, maxit)
    n = numel(g);
    s = zeros(n,1); Hs = zeros(n,1);
    r = g; z = Minv(r); d = -z; rz = r'*z;
    sMs = 0; sMd = 0; dMd = rz;
    Sp = zeros(n,0); Yp = zeros(n,0);
    info.flag = 'maxit'; info.boundary = false;
    for j = 1:maxit
        Hd = Hfun(d);
        kappa = d'*Hd;
        if kappa <= 0
            tau = toBoundary(sMs, sMd, dMd, Delta);
            s = s + tau*d; Hs = Hs + tau*Hd;
            info.flag = 'negcurv'; info.boundary = true; info.iter = j;
            info.sMnorm = Delta;
            return
        end
        Sp(:,end+1) = d; Yp(:,end+1) = Hd; %#ok<AGROW>
        alpha = rz/kappa;
        sMs_new = sMs + 2*alpha*sMd + alpha^2*dMd;
        if sMs_new >= Delta^2
            tau = toBoundary(sMs, sMd, dMd, Delta);
            s = s + tau*d; Hs = Hs + tau*Hd;
            info.flag = 'boundary'; info.boundary = true; info.iter = j;
            info.sMnorm = Delta;
            return
        end
        s = s + alpha*d; Hs = Hs + alpha*Hd; sMs = sMs_new;
        r = r + alpha*Hd;
        if norm(r) <= tol
            info.flag = 'conv'; info.iter = j; info.sMnorm = sqrt(sMs);
            return
        end
        z = Minv(r); rz_new = r'*z; beta = rz_new/rz;
        sMd = beta*(sMd + alpha*dMd);
        dMd = rz_new + beta^2*dMd;
        d = -z + beta*d; rz = rz_new;
    end
    info.iter = maxit; info.sMnorm = sqrt(sMs);
end

function tau = toBoundary(sMs, sMd, dMd, Delta)
    tau = (-sMd + sqrt(sMd^2 + dMd*(Delta^2 - sMs)))/dMd;
end
