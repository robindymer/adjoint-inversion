function run_method(method, maxTime, nThreads, m)
    if nargin < 4, m = 36; end
    benchDir = fileparts(mfilename('fullpath'));
    cd(fileparts(benchDir));
    addSubpaths;
    maxNumCompThreads(nThreads);
    [optObj, x0, ~, xtrue] = setup_problem(m, 26, 0.0135);
    P = InvProblem(optObj, xtrue);
    P.maxTime = maxTime;
    P.maxCost = 250;
    lb = zeros(size(x0));
    fprintf('==== %s ====\n', method);
    switch method
        case 'fmincon_ip_lbfgs' % current 'lbfgs' option
            o = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Algorithm', 'interior-point', ...
                'HessianApproximation', 'lbfgs', 'MaxIterations', 1000, 'OptimalityTolerance', 1e-10, ...
                'Display', 'iter', 'OutputFcn', @(x,ov,st) outfun(P,x,ov,st));
            fmincon(@(x) P.fg(x), x0, [],[],[],[], lb, [], [], o);
        case 'fmincon_ip_cg' % current 'CG' option
            o = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Algorithm', 'interior-point', ...
                'SubproblemAlgorithm', 'cg', 'HessianMultiplyFcn', @(x,lambda,v) P.Hv(x,v), ...
                'MaxIterations', 1000, 'OptimalityTolerance', 1e-10, ...
                'Display', 'iter', 'OutputFcn', @(x,ov,st) outfun(P,x,ov,st));
            fmincon(@(x) P.fg(x), x0, [],[],[],[], lb, [], [], o);
        case 'fmincon_trr' % built-in trust-region Newton-CG with Hv, bounds via reflection
            o = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Algorithm', 'trust-region-reflective', ...
                'HessianMultiplyFcn', @(Hinfo,Y) P.Hv(Hinfo,Y), 'SubproblemAlgorithm', 'cg', ...
                'MaxIterations', 1000, 'OptimalityTolerance', 1e-10, 'FunctionTolerance', 1e-14, ...
                'Display', 'iter', 'OutputFcn', @(x,ov,st) outfun(P,x,ov,st));
            fmincon(@(x) fgx(P,x), x0, [],[],[],[], lb, [], [], o);
        case 'fminunc_qn' % BFGS, no bounds (f = inf for a <= 0)
            o = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'Algorithm', 'quasi-newton', ...
                'MaxIterations', 1000, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-14, ...
                'Display', 'iter', 'OutputFcn', @(x,ov,st) outfun(P,x,ov,st));
            fminunc(@(x) fg_nonan(P,x), x0, o);
        case 'newton_tr_repo'
            o.maxIter = 1000;
            o.OutputFcn = @(x, f, g, info) logStop(P, x, f, g, info);
            adjopt.trustRegionNewtonCG(@(x) P.fg(x), @(x,V) P.Hv(x,V), x0, lb, o);
        case 'lbfgs'
            lbfgs_ls(P, x0, struct('maxIter', 1000));
        case 'trncg'
            tr_ncg(P, x0, struct('maxIter', 1000));
        case 'trncg_pc'
            tr_ncg(P, x0, struct('maxIter', 1000, 'precond', true));
        case 'trncg_pc_switch'
            tr_ncg(P, x0, struct('maxIter', 1000, 'precond', true, 'switchToLBFGS', true));
        otherwise
            error('unknown method %s', method);
    end
    hist = P.hist;
    save(fullfile(benchDir, 'results', sprintf('res_%s_m%d.mat', method, m)), 'hist', 'method');
    fprintf('==== done %s: f %g err %g cost %g time %g\n', method, hist.f(end), hist.err(end), P.cost(), hist.t(end));
end

function stop = outfun(P, x, ov, state)
    stop = false;
    if strcmp(state, 'iter') || strcmp(state, 'init')
        g = ov.gradient; if isempty(g), g = nan; end
        P.log(x, ov.fval, g, state);
        stop = P.timeUp();
    end
end

function [f, g, Hinfo] = fgx(P, x)
    [f, g] = P.fg(x); Hinfo = x;
end

function [f, g] = fg_nonan(P, x)
    [f, g] = P.fg(x); g(isnan(g)) = 0;
end

function stop = logStop(P, x, f, g, info)
    P.log(x, f, g, info);
    stop = P.timeUp();
end
