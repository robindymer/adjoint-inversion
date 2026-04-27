function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = d2_variable_hollow_4(m,h)

    BP = 6;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end




    % Norm
    Hv = ones(m,1);
    Hv(1:4) = [17/48 59/48 43/48 49/48];
    Hv(m-3:m) = rot90(Hv(1:4),2);
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Boundary operators
    e_l = sparse(m,1);
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_l(1:4) = 1/h*[-11/6 3 -3/2 1/3];
    d1_r = -rot90(d1_l, 2);




    S = d1_l*d1_l' + d1_r*d1_r';

    stencil = [1/12 -2/3 0 2/3 -1/12];
    diags = -2:2;

    Q_U = [
        0 0.59e2/0.96e2 -0.1e1/0.12e2 -0.1e1/0.32e2;
        -0.59e2/0.96e2 0 0.59e2/0.96e2 0;
        0.1e1/0.12e2 -0.59e2/0.96e2 0 0.59e2/0.96e2;
        0.1e1/0.32e2 0 -0.59e2/0.96e2 0;
    ];

    Q = stripeMatrix(stencil, diags, m);
    Q(1:4,1:4) = Q_U;
    Q(m-3:m,m-3:m) = -rot90(Q_U, 2);

    D1 = HI*(Q - 1/2*e_l*e_l' + 1/2*e_r*e_r');


    % Second derivative
    nBP = 6;
    M = sparse(m,m);
    coeffs = load('sbplib/+sbp/+implementations/coeffs_d2_variable_4.mat');

    function D2 = D2_fun(c)
        M_l = zeros(nBP, coeffs.nBPC);
        M_r = zeros(nBP, coeffs.nBPC);

        for i=1:coeffs.nBPC
            M_l = M_l + coeffs.C_l{i}*c(i);
            M_r = M_r + coeffs.C_r{i}*c(m-coeffs.nBPC+i);
        end

        M(1:nBP, 1:coeffs.nBPC) = M_l;
        M(m-nBP+1:m, m-coeffs.nBPC+1:m) = M_r;

        D2 = M/h^2;
    end
    D2 = @D2_fun;

end