function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_periodic_5(m,h)
    
    if(m<8)
        error('Operator requires at least 8 grid points');
    end

    Hv=ones(m,1);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-2 -1 0 1 2 3];
    q_stencil = [1/20 -1/2 -1/3 +1 -1/4 +1/30];
    Qp = stripeMatrixPeriodic(q_stencil, q_diags,m);

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp) ;

    Dm=HI*(Qm) ;
end
