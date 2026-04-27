function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_periodic_6(m,h)
    
    if(m<12)
        error('Operator requires at least 12 grid points');
    end

    Hv = ones(m,1);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-2 -1 0 1 2 3 4];
    q_stencil = [1/30 -2/5 -7/12 4/3 -1/2 2/15 -1/60];
    Qp = stripeMatrixPeriodic(q_stencil, q_diags,m);

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp) ;

    Dm=HI*(Qm) ;

end
