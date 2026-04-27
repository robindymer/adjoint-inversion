function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = d2_variable_hollow_2(m,h)

    BP = 1;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    % Norm
    Hv = ones(m,1);
    Hv(1) = 1/2;
    Hv(m:m) = 1/2;
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Boundary operators
    e_l = sparse(m,1);
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_l(1:3) = 1/h*[-3/2 2 -1/2];
    d1_r = -rot90(d1_l, 2);

    % D1 operator
    diags   = -1:1;
    stencil = [-1/2 0 1/2];
    D1 = stripeMatrix(stencil, diags, m);

    D1(1,1)=-1;D1(1,2)=1;D1(m,m-1)=-1;D1(m,m)=1;
    D1(m,m-1)=-1;D1(m,m)=1;
    D1=D1/h;
    %Q=H*D1 + 1/2*(e_1*e_1') - 1/2*(e_m*e_m');


    nBP = 2;
    M = sparse(m,m);
    coeffs = load('sbplib/+sbp/+implementations/coeffs_d2_variable_2.mat');

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