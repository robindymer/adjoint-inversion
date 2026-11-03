% Summary table of saved inversions (metrics in notes/high_res_plan.md).
% pattern -- glob for the result files in mat/inversion, e.g.
%            matlab -batch "pattern='rough_*_m76_mp26_velocity_synth_*'; summarize_inversions"
% For each run: final misfit F and regularization R, the exact objective F + R with eps = 0 for
% 'tv' and 'lasso', the relative error in a over all nodes and over the resolved nodes 6-25,
% iterations, wall time, and the wall time to reach (J - J*)/(J0 - J*) = 1e-1, 1e-2, 1e-3, where J
% is the exact objective and J* the smallest J over all listed runs that use the same regularization.
% High-res runs have no true parameters saved, so a_true is taken from an inverse-crime run.
addSubpaths;
if ~exist('pattern', 'var')
    pattern = 'rough_*_m76_mp26_velocity_*';
end
resolved = 6:25;
levels = [1e-1, 1e-2, 1e-3];

files = dir(fullfile('mat/inversion', [pattern, '.mat']));
files = files(~contains({files.name}, 'checkpoint'));
n = numel(files);
runs = [];
aTrueSynth = [];
for i = 1:n
    path = fullfile(files(i).folder, files(i).name);
    vars = intersect({'history', 'truePars', 'regOpts', 'initialPars'}, who('-file', path));
    S = load(path, vars{:});
    h = S.history;
    if ~isfield(h, 'reg') % Runs from before the regularization was added
        h.reg = zeros(size(h.fval));
    end
    if ~isfield(S, 'regOpts')
        S.regOpts = struct('type', 'none', 'beta', 0);
    end
    regOpts = S.regOpts;
    default_field(regOpts, 'eps', 1e-4);

    % Exact objective along the history, with eps = 0 for 'tv' and 'lasso'
    p0 = S.initialPars.a(:);
    regExact = adjopt.regularization(setfield(regOpts, 'eps', 0), p0, numel(p0));
    misfit = h.fval - h.reg;
    Jhist = misfit + arrayfun(@(k) regExact.R(h.x(:,k)), (1:size(h.x,2))');

    if ~isempty(S.truePars)
        aTrueSynth = S.truePars.a(:);
    end
    r.name = erase(files(i).name, '.mat');
    r.regKey = sprintf('%s_%g', regOpts.type, regOpts.beta);
    r.F = misfit(end);
    r.R = h.reg(end);
    r.Jexact = Jhist(end);
    r.a = h.x(:,end);
    r.iter = numel(h.fval) - 1;
    r.T = h.t(end)/3600;
    r.Jhist = Jhist;
    r.t = h.t;
    if i == 1
        runs = r;
    else
        runs(i) = r;
    end
end
if isempty(aTrueSynth)
    S = load('mat/inversion/rough_newton-cg-ls_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat', 'truePars');
    aTrueSynth = S.truePars.a(:);
end

fprintf('%-72s %9s %9s %9s %7s %7s %5s %6s | t (h) to rel. level %s\n', 'run', 'F', 'R', ...
        'F+R exact', 'err', 'err6-25', 'iter', 'T (h)', mat2str(levels));
for i = 1:n
    r = runs(i);
    sameReg = strcmp({runs.regKey}, r.regKey);
    Jstar = min(cellfun(@min, {runs(sameReg).Jhist}));
    relJ = (r.Jhist - Jstar)/(r.Jhist(1) - Jstar);
    tLevel = '';
    for L = levels
        k = find(relJ <= L, 1);
        if isempty(k)
            tLevel = [tLevel, sprintf(' %6s', '-')];
        else
            tLevel = [tLevel, sprintf(' %6.2f', r.t(k)/3600)];
        end
    end
    err = norm(r.a - aTrueSynth)/norm(aTrueSynth);
    errRes = norm(r.a(resolved) - aTrueSynth(resolved))/norm(aTrueSynth(resolved));
    fprintf('%-72s %9.3e %9.3e %9.3e %7.4f %7.4f %5d %6.2f |%s\n', r.name, r.F, r.R, r.Jexact, ...
            err, errRes, r.iter, r.T, tLevel);
end
