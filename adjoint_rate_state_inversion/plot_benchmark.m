S = load('optimizer_benchmark/results/res_newton_tr_repo_m251.mat');
h = S.hist;
m_p = size(h.x,1);
x = linspace(-10,10,m_p);   % same parameter grid as direct_effect_inverse_crime.m
vw = (x >= -5) & (x <= 6);
a_true = 0.009*vw + 0.013*~vw;
b = 0.011*ones(size(x));
cost = (4.3*h.nFwd + 1.3*h.nAdj + 4.1*h.nHv)/5.6;   % gradient-equivalents, as in summarize.m

% Parameter profile
figure();
plot(x,a_true, x,h.x(:,1),':x', x,h.x(:,end),'k--x', x,b,'k:','linewidth',2);
xlabel('$x$ (km)','interpreter','latex');
legend('$\mathbf{a}$, true','$\mathbf{a}$, initial','$\mathbf{a}$, final','$\mathbf{b}$','interpreter','latex','Location','southwest');
title(strrep(S.method,'_','\_')); grid on; set(gca,'fontsize',16);

% Convergence
figure();
tiledlayout(1,3);
nexttile; semilogy(cost, h.f, '-o'); xlabel('cost (grad-equiv.)'); ylabel('misfit f'); grid on;
nexttile; semilogy(cost, h.gnorm, '-o'); xlabel('cost (grad-equiv.)'); ylabel('|g|'); grid on;
nexttile; plot(cost, h.err, '-o'); xlabel('cost (grad-equiv.)'); ylabel('rel. error in a'); grid on;
