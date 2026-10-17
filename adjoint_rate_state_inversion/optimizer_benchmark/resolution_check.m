benchDir = fileparts(mfilename('fullpath'));
cd(fileparts(benchDir));
addSubpaths; maxNumCompThreads(3);
SP = fullfile(benchDir, 'results');
ms = [36 76 151 251];
R = cell(size(ms)); Tt = cell(size(ms));
for i = 1:numel(ms)
    o.initialGuessScalings = []; o.initialGuessValues = 0.0135; o.misfitType = 'velocity';
    o.inversionParameters = {'a'}; o.m = ms(i); o.m_p = 26;
    [~, tp] = pars.rsFriction2DFractalFaultInversion(o);
    ob = AntiplaneShear2DRSFrictionOptHessian(tp, {'a'}, 8);
    tic; ob.runForward(false); tf = toc;
    R{i} = ob.forwardReceiverRecordings; Tt{i} = ob.forwardTimeIntegrationData.T; Ht = ob.Ht;
    fprintf('m=%d forward %.1f s, k=%g, nT=%d, max|v_rec|=%g\n', ms(i), tf, ob.k, numel(Tt{i}), max(abs(R{i}(:))));
    save(fullfile(SP,'resolution_data.mat'), 'R', 'Tt', 'ms');
    clear ob
end
% All m use the same time levels (fixed k), including repeated RK stage times,
% so recordings are compared directly. Integrate over unique times only.
ref = R{end}; Tref = Tt{end};
[~, idx] = unique(Tref);
w = zeros(size(Tref)); w(idx) = gradient(Tref(idx));
for i = 1:numel(ms)-1
    assert(isequal(Tt{i}, Tref), 'time levels differ between grids');
    d = R{i} - ref;
    fprintf('m=%d vs m=%d: rel L2 diff of recordings %.3e, misfit-equivalent %.3e\n', ...
        ms(i), ms(end), norm(d(:))/norm(ref(:)), 0.5*sum((d(:,idx).^2)*w(idx)'));
end
fprintf('data 0.5*int v^2: %.3e\n', 0.5*sum((ref(:,idx).^2)*w(idx)'));
