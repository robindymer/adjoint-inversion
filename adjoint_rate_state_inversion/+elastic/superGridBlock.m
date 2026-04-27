
%   xPars.side              -- 'l' or 'r' or for stretching left/right sides.
%   xPars.x0                -- Coordinate where layer attaches to non-physical domain.
%   xPars.measuresAtEnds    -- Cell array of "blow-up" of reference domain at different ends of layer
function ti = superGridBlock(xPars, yPars, lines)
    default_arg('x0', 0);
    default_arg('y0', 0);

    defaultPars = struct;
    defaultPars.x0 = 0;
    defaultPars.side = 'r';

    % No stretching by default
    defaultPars.measuresAtEnds = {1,1};

    default_struct('xPars', defaultPars);
    default_struct('yPars', defaultPars);

    lines_default = struct;
    lines_default.L = [];
    lines_default.B = [];
    lines_default.R = [];
    lines_default.T = [];
    default_struct('lines', lines_default);

    x0 = xPars.x0;
    y0 = yPars.x0;
    mAEX = xPars.measuresAtEnds;
    mAEY = yPars.measuresAtEnds;

    x = polynomialStretch(mAEX{1}, mAEX{2});
    y = polynomialStretch(mAEY{1}, mAEY{2});

    switch xPars.side
    case 'r'
        x = @(xi) x0 + x(xi);
    case 'l'
        x = @(xi) x0 - x(1-xi);
    end

    switch yPars.side
    case 'r'
        y = @(xi) y0 + y(xi);
    case 'l'
        y = @(xi) y0 - y(1-xi);
    end

    if isempty(lines.L)
        line_l = @(xi) [x(0)+0*xi; y(xi)];
        line_l = parametrization.Curve(line_l).reverse;
    else
        line_l = lines.L;
    end

    if isempty(lines.B)
        line_b = @(xi) [x(1-xi); y(0)+0*xi];
        line_b = parametrization.Curve(line_b).reverse;
    else
        line_b = lines.B;
    end

    if isempty(lines.T)
        line_t = @(xi) [x(xi); y(1)+0*xi];
        line_t = parametrization.Curve(line_t).reverse;
    else
        line_t = lines.T;
    end

    if isempty(lines.R)
        line_r = @(xi) [x(1)+0*xi; y(1-xi)];
        line_r = parametrization.Curve(line_r).reverse;
    else
        line_r = lines.R;
    end

    ti = parametrization.Ti(line_b, line_r, line_t, line_l);

end


% Creates stretching x = x(xi), where x is a polynomial of degree
% deg + 1 in xi.
%
% The stretching function has as many zero derivatives (at both ends)
% as are allowed by deg.
function x = polynomialStretch(x_xi_min, x_xi_max, deg)

    default_arg('deg', []);

    % Stretching function from 0 to 1
    x_unit = stretchingPolynomial(deg);

    x = @(xi) x_xi_min*xi + (x_xi_max - x_xi_min)*x_unit(xi);

end

function P = stretchingPolynomial(deg)
    default_arg('deg', 3);

    % Load MAT-file with polynomial if it exists,
    % otherwise create it.
    filename = ['stretchingPol_deg' num2str(deg) '.mat'];
    if (exist(filename,'file') == 2 )
        load(filename);
        return;
    end

    syms x
    a = sym('a', [1, deg + 2]);

    P = symfun(poly2sym(a,x), x);
    p = diff(P,x);
    px = diff(p,x);

    % Arbitrary base line
    eq_P = P(0) == 0;

    % Condition on grid spacing
    eq_p = [p(0) == 0; p(1) == 1];

    % Continuity of derivatives
    p_der = p;
    eq_px = [];
    for i = 1:(deg-1)/2
        p_der = diff(p_der,x);
        eq_px = [eq_px; p_der(0) == 0; p_der(1) == 0];
    end

    eq = [eq_P; eq_p; eq_px];
    S = solve(eq, a, 'ReturnConditions', true);
    DOF = remainingPar(eq, a);
    % fprintf('Remaining free parameters in a: %d \n',DOF);
    [a, var] = mySubs(a, S, a);

    % Construct polynomials
    P = symfun(poly2sym(a,x), x);
    P = matlabFunctionSizePreserving(P);

    % Save to .mat-file
    save(filename, 'P');

end

function DOF = remainingPar(eq, pars)
    [A, ~] = equationsToMatrix(eq, pars);
    nullSpace = null(A);
    DOF = size(nullSpace,2);
end

function [A,var] = mySubs(A,S,var)
    %Substitutes variables in var for variables in S
    %in symbolic object A.
    % Typically, S is the result of S=solve(...);

    if(isempty(S.conditions))
        disp('OBS! Matlab didnt find a solution...');
    end

    par = S.parameters;
    for i = 1:numel(var)
        var_str = char(var(i));
        callstr = ['A = subs(A,var(i),S.' var_str ');'];
        eval(callstr);
    end
    var = par;
end