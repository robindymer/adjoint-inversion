function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_periodic_2(m,h)
    
    if(m<4)
        error('Operator requires at least 4 grid points');
    end

    Hv=ones(m,1);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [0 1 2];
    q_stencil = [-3/2 +2 -1/2];
    Qp = stripeMatrixPeriodic(q_stencil, q_diags,m);

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp) ;

    Dm=HI*(Qm) ;
end

