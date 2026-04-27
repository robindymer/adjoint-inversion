function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = d2_variable_periodic_8(m,h)
    % m = number of unique grid points, i.e. h = L/m;

    if(m<9)
        error(['Operator requires at least ' num2str(8) ' grid points']);
    end

    % Norm
    Hv = ones(m,1);
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Dummy boundary operators
    e_l = sparse(m,1);
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_r = -rot90(d1_l, 2);


    % D1 operator
    diags   = -4:4;
    stencil = [1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280];
    D1 = stripeMatrixPeriodic(stencil, diags, m);
    D1 = D1/h;

    % Undivided differences
    diags   = -3:2;
    stencil = [-1 5 -10 10 -5 1];
    DD_5 = stripeMatrixPeriodic(stencil, diags, m);

    diags   = -3:3;
    stencil = [1 -6 15 -20 15 -6 1];
    DD_6 = stripeMatrixPeriodic(stencil, diags, m);

    diags   = -4:3;
    stencil = [-1 7 -21 35 -35 21 -7 1];
    DD_7 = stripeMatrixPeriodic(stencil, diags, m);

    diags   = -4:4;
    stencil = [1 -8 28 -56 70 -56 28 -8 1];
    DD_8 = stripeMatrixPeriodic(stencil, diags, m);

    % D2 operator
    function D2 = D2_fun(c)

        C1 = sparse(diag(c));
        
        diags   = -1:0;
        stencil = [1/2, 1/2];
        C2 = stripeMatrixPeriodic(stencil, diags, m);

        diags   = -1:1;
        stencil = [3/10, 2/5, 3/10];
        C3 = stripeMatrixPeriodic(stencil, diags, m);

        diags   = -2:1;
        stencil = [1/4, 1/4, 1/4, 1/4];
        C4 = stripeMatrixPeriodic(stencil, diags, m);

        C2 = sparse(diag(C2 * c));
        C3 = sparse(diag(C3 * c));
        C4 = sparse(diag(C4 * c));

        % Remainder term added to wide second derivative operator
        R = (1/78400 / h) * transpose(DD_8) * C1 * DD_8 + (1/14700 / h) * transpose(DD_7) * C2 * DD_7 + (1/2520 / h) * transpose(DD_6) * C3 * DD_6 + (1/350 / h) * transpose(DD_5) * C4 * DD_5;

        D2 = D1 * C1 * D1 - H \ R;
    end
    D2 = @D2_fun;


end
