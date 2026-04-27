function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_periodic_3(m,h)

    if(m<6)
        error('Operator requires at least 6 grid points');
    end

    Hv = ones(m,1);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-1, 0, 1, 2];
    q_stencil = [-1/3 -1/2 1 -1/6];
    Qp = stripeMatrixPeriodic(q_stencil, q_diags,m);

    Qm=-Qp';

    e_1=sparse(m,1);
    e_1(1)=1;
    e_m=sparse(m,1);
    e_m(m)=1;

    Dp=HI*(Qp) ;

    Dm=HI*(Qm) ;
end

