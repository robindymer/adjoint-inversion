benchDir = fileparts(mfilename('fullpath'));
cd(fileparts(benchDir));
addSubpaths;
maxNumCompThreads
[optObj, x0, truePars, xtrue] = setup_problem(36, 26, 0.0135);
tic; [f, g] = adjopt.vectorGradient(x0, optObj); tg = toc;
fprintf('grad: %.2f s, f = %g, |g| = %g\n', tg, f, norm(g));
v = randn(size(x0)); v = v/norm(v)*1e-3;
tic; Hv = adjopt.hessianVector(x0, v, optObj); th = toc;
fprintf('Hv (cached 1st order): %.2f s\n', th);
tic; optObj.setParameters(optObj.parametersVectorToNative(x0)); optObj.updateForwardDiscr(); optObj.runForward(false); f2 = optObj.computeMisfit(); tf = toc;
fprintf('fval only: %.2f s (f=%g)\n', tf, f2);
tic; Hv2 = adjopt.hessianVector(x0, v, optObj); th2 = toc;
fprintf('Hv (uncached): %.2f s, rel diff %g\n', th2, norm(Hv-Hv2)/norm(Hv));
% FD check of Hv against gradients
h = 1e-1;
[~, gp] = adjopt.vectorGradient(x0 + h*v, optObj);
[~, gm] = adjopt.vectorGradient(x0 - h*v, optObj);
fprintf('Hv vs central FD of grad: rel err %g\n', norm(Hv - (gp-gm)/(2*h))/norm(Hv));
% Full Hessian at x0 (26 Hv)
tic; Hfull = adjopt.hessianVector(x0, eye(numel(x0)), optObj); tH = toc;
fprintf('Full Hessian: %.1f s, asym %g\n', tH, norm(Hfull-Hfull','fro')/norm(Hfull,'fro'));
e = eig((Hfull+Hfull')/2); fprintf('eig: min %g max %g, #neg %d, cond %g\n', min(e), max(e), sum(e<0), max(abs(e))/min(abs(e)));
disp(e')
save(fullfile(benchDir,'results','H0.mat'),'Hfull','g','x0','xtrue','f');
