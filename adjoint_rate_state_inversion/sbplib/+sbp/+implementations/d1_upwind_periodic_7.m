function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_periodic_7(m,h)

    BP = 6;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    Hv=ones(m,1);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags = [-3 -2 -1 0 1 2 3 4];
    q_stencil = [-1/105 +1/10 -3/5 -1/4 +1 -3/10 +1/15 -1/140];
    Qp = stripeMatrixPeriodic(q_stencil, q_diags,m);

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp) ;

    Dm=HI*(Qm) ;
end
