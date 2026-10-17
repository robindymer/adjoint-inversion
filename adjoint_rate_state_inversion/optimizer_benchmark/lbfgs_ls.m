% L-BFGS with backtracking (quadratic/cubic, Dennis-Schnabel) line search and fraction-to-boundary for x > 0.
% S0, Y0: optional seed pairs (e.g. from Newton-CG iterations).
function [x, P] = lbfgs_ls(P, x, opts, S, Y)
    default_field(opts, 'maxIter', 100);
    default_field(opts, 'memL', 20);
    default_field(opts, 'gtol', 1e-10);
    default_field(opts, 'ftb', 0.9);
    default_field(opts, 'step0', 0.1*norm(x)); % length of the first (steepest descent) step
    default_arg('S', zeros(numel(x),0));
    default_arg('Y', zeros(numel(x),0));
    [f, g] = P.fg(x);
    if isempty(P.hist.f), P.log(x, f, g, 'init'); end
    for it = 1:opts.maxIter
        if P.timeUp() || norm(g) < opts.gtol, break, end
        if isempty(S)
            d = -g*(opts.step0/norm(g));
        else
            d = -lbfgs_twoloop(g, S, Y);
        end
        gd = g'*d;
        if gd >= 0
            d = -g*(opts.step0/norm(g)); gd = g'*d;
        end
        neg = d < 0;
        stp = min([1; opts.ftb*x(neg)./(-d(neg))]);
        [stp, ft, nls] = backtrack(P, x, d, f, gd, stp);
        if isempty(ft)
            fprintf('line search failed\n'); break
        end
        xn = x + stp*d;
        [fn, gn] = P.fg(xn);
        s = xn - x; y = gn - g;
        if s'*y > 1e-12*norm(s)*norm(y)
            S = [S, s]; Y = [Y, y];
            S = S(:, max(1,end-opts.memL+1):end); Y = Y(:, max(1,end-opts.memL+1):end);
        end
        x = xn; f = fn; g = gn;
        P.log(x, f, g, sprintf('stp %8.2e ls %d', stp, nls));
    end
end

function [stp, f, iter] = backtrack(P, x, d, f0, dg, stp)
    ftol = 1e-4;
    iter = 0; fp = []; stp_p = [];
    while stp > 1e-10
        f = P.f(x + stp*d); iter = iter + 1;
        if f < f0 + ftol*stp*dg
            return
        end
        if isempty(fp) || ~isfinite(f) || ~isfinite(fp)
            stn = -dg*stp^2/(2*(f - f0 - dg*stp));
            if ~isfinite(f), stn = 0.1*stp; end
        else
            t1 = f - f0 - stp*dg; t2 = fp - f0 - stp_p*dg; t3 = 1/(stp - stp_p);
            a = t3*(t1/stp^2 - t2/stp_p^2); b = t3*(t2*stp/stp_p^2 - t1*stp_p/stp^2);
            if a ~= 0, stn = (-b + sqrt(max(b^2 - 3*a*dg, 0)))/(3*a); else, stn = -dg/(2*b); end
        end
        stp_p = stp; fp = f;
        stp = min(max(stn, 0.1*stp), 0.5*stp);
    end
    f = [];
end
