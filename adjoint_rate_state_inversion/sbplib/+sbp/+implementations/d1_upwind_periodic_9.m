function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_periodic_9(m,h)

    BP = 8;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    Hv=ones(m,1);
    Hv=Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-4 -3 -2 -1 0 1 2 3 4 5];
    q_stencil = [1/504 -1/42 +1/7 -2/3 -1/5 +1 -1/3 +2/21 -1/56 +1/630];
    Qp = stripeMatrixPeriodic(q_stencil, q_diags,m);

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp) ;

    Dm=HI*(Qm) ;
end
