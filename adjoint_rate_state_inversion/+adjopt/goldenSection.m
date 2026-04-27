function adjopt = goldenSection(adjopt, gradient, L, tol)
% Minimize function in 1D using golden section search
% L:    Length of first search interval
% tol:  Relative tolerance. Search stops when interval is less than tol*L;
default_arg('tol', 1e-3);

% Store initial interval length for printing.
L0 = L;

% Left and right end points of search interval
a = 0;
b = L;

% Golden section coefficient, tau is approx 0.6
tau = (-1+sqrt(5))/2;

% Location of first probe point (about 40% of interval)
alpha_l = (1-tau)*L;

% Location of second probe point (about 60% of interval)
alpha_r = tau*L;

pars = adjopt.pars;

% Evaluate misfit at left end point and set dummy value for right end points
f_a = adjopt.computeMisfit();
f_b = 1.1*f_a;

% Evaluate misfits at left probe point
adjopt.updateParameters(pars, gradient, -alpha_l);
adjopt.updateForwardDiscr()
adjopt.runForward();
f_l = adjopt.computeMisfit();

% Step to right probe point and evaluate misfit
adjopt.updateParameters(pars, gradient, -alpha_r);
adjopt.updateForwardDiscr()
adjopt.runForward();
f_r = adjopt.computeMisfit();

k = 0;
while L/L0>tol

  k = k+1;

  if f_l<f_r
    % Minimum bracketed in (a, alpha_r)
    % Move right endpoint
    b = alpha_r;
    f_b = f_r;
    L = b-a;

    % Use the left probe point as right probe point in new interval
    alpha_r = alpha_l;
    f_r = f_l;

    % New left probe point
    alpha_l = a + (1-tau)*L;
    adjopt.updateParameters(pars, gradient, -alpha_l);
    adjopt.updateForwardDiscr()
    adjopt.runForward();
    f_l = adjopt.computeMisfit();

  else
    % Minimum bracketed in (alpha_l, b)
    % Move left endpoint
    a = alpha_l;
    f_a = f_l;
    L = b-a;

    % Use the right probe point as left probe point in new interval
    alpha_l = alpha_r;
    f_l = f_r;

    % New right probe point
    alpha_r = a + tau*L;
    adjopt.updateParameters(pars, gradient, -alpha_r);
    adjopt.updateForwardDiscr()
    adjopt.runForward();
    f_r = adjopt.computeMisfit();
  end
end
% stopping criterion satisfied

% Choose point with smallest function value of those evaluated
f_values = [f_a, f_l, f_r, f_b];
alpha_values = [a, alpha_l, alpha_r, b];
[~, i] = min(f_values);
alpha = alpha_values(i);

adjopt.updateParameters(pars, gradient, -alpha);
adjopt.updateForwardDiscr()

fprintf('GS iter: %d, interval length: %f, rel step length: %f \n', k, L0, alpha/L0);

end
