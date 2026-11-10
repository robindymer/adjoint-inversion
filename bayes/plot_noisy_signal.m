% Plot noisy vs clean high-res receiver data (output of add_noise.m)

dataDir = '../adjoint_rate_state_inversion/mat/fractal_fault_m1001_mp1001';
eta     = 0.02;
recs    = [1 100 200 300];   % receivers to plot

warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
S = load(fullfile(dataDir, sprintf('receiverData_noisy_eta%g.mat', eta)));
noise = S.receiverData.noise;

figure;
tiledlayout(numel(recs), 1, 'TileSpacing', 'compact');
for k = 1:numel(recs)
    r = recs(k);
    nexttile;
    plot(noise.t, noise.d_noisy(r,:), '.', 'MarkerSize', 4); hold on
    plot(noise.t, noise.d_clean(r,:), 'LineWidth', 1.5); hold off
    ylabel(sprintf('Rec %d', r));
    xlim([noise.t(1), noise.t(end)]);
end
xlabel('t [s]');
legend({'Noisy', 'Clean'}, 'Location', 'best');
sgtitle(sprintf('\\sigma = %.2e (\\eta = %g)', noise.sigma, noise.eta));
