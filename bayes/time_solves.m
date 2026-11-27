% Wall-clock time of one forward, adjoint, second-order forward and second-order adjoint solve at a_MAP,
% for the inverse crime (m = 36) and the high-res setup (m = 151), median over nRep repetitions.
% Run alone on the machine: matlab -batch "time_solves"
bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('nRep', 'var'), nRep = 3; end
resDir = fullfile(bayesDir, 'results');
cases = {36, 'hessian_m36_eta0.02.mat'; 151, 'highres_map_m151_mp26_eta0.02.mat'};
stages = {'forward', 'adjoint', '2nd-order forward', '2nd-order adjoint'};
T = zeros(size(cases, 1), numel(stages));

for c = 1:size(cases, 1)
    m = cases{c, 1};
    L = load(fullfile(resDir, cases{c, 2}), 'aMap');
    if m == 36
        prob = setup_bayes_problem(m, 26, 0.02);
    else
        prob = setup_highres_problem(m, 26, 0.02);
    end
    o = prob.obj;
    v = randn(26, 1);
    t = zeros(nRep, numel(stages));
    for r = 1:nRep
        % Each time includes the discretization update preceding the solve
        tic; prob.setX(L.aMap); o.runForward(false); t(r,1) = toc;
        tic; o.updateAdjointDiscr(); o.runAdjoint(false); t(r,2) = toc;
        o.firstOrderParsVec = o.parametersNativeToVector(o.forwardDiscr.friction.rsParams);
        o.setPerturbationDirection(o.parametersVectorToNative(v));
        tic; o.updateSecondOrderForwardDiscr(); o.runSecondOrderForward(false, o.T, [], false); t(r,3) = toc;
        tic; o.updateSecondOrderAdjointDiscr(); o.runSecondOrderAdjoint(false, o.T, [], false); t(r,4) = toc;
        fprintf('m = %d, rep %d: %s\n', m, r, sprintf('%.1f s  ', t(r,:)));
    end
    T(c,:) = median(t, 1);
end

fprintf('\nMedian wall-clock time per solve [s]\n%6s', 'm');
fprintf('%20s', stages{:}); fprintf('\n');
for c = 1:size(cases, 1)
    fprintf('%6d', cases{c, 1}); fprintf('%20.1f', T(c,:)); fprintf('\n');
end
save(fullfile(resDir, 'solve_times.mat'), 'T', 'stages', 'cases', 'nRep');
