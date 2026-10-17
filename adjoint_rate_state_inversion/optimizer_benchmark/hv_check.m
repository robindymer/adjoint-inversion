function hv_check(m, nThreads, tol)
    % Finite-difference check of Hv at grid size m, at the initial guess a = 0.0135.
    % Compares Hv with central differences of the gradient for several step sizes.
    % Errors if the best relative error exceeds tol, so a dependent Slurm job (afterok) does not start.
    if nargin < 2, nThreads = maxNumCompThreads; end
    if nargin < 3, tol = 1e-3; end
    benchDir = fileparts(mfilename('fullpath'));
    cd(fileparts(benchDir));
    addSubpaths;
    maxNumCompThreads(nThreads);
    fprintf('==== Hv FD check, m = %d, %d threads ====\n', m, maxNumCompThreads);
    tic; [optObj, x0] = setup_problem(m, 26, 0.0135); fprintf('setup: %.1f s\n', toc);
    tic; [f, g] = adjopt.vectorGradient(x0, optObj); tg = toc;
    fprintf('grad: %.1f s, f = %g, |g| = %g\n', tg, f, norm(g));
    rng(0);
    v = randn(size(x0)); v = v/norm(v)*1e-3;
    tic; Hv = adjopt.hessianVector(x0, v, optObj); th = toc;
    fprintf('Hv (cached 1st order): %.1f s, |Hv| = %g\n', th, norm(Hv));
    hs = [1, 1e-1, 1e-2];
    relErr = zeros(size(hs));
    for i = 1:numel(hs)
        h = hs(i);
        [~, gp] = adjopt.vectorGradient(x0 + h*v, optObj);
        [~, gm] = adjopt.vectorGradient(x0 - h*v, optObj);
        relErr(i) = norm(Hv - (gp - gm)/(2*h))/norm(Hv);
        fprintf('h = %g (|step| = %g): rel err %g\n', h, h*norm(v), relErr(i));
    end
    save(fullfile(benchDir, 'results', sprintf('hv_check_m%d.mat', m)), 'Hv', 'v', 'hs', 'relErr', 'g', 'f', 'x0');
    if min(relErr) > tol
        error('Hv FD check failed at m = %d: best rel err %g > %g', m, min(relErr), tol);
    end
    fprintf('==== Hv FD check passed: best rel err %g ====\n', min(relErr));
end
