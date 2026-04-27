function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_periodic_8(m,h)
    
    BP = 8;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    Hv=ones(m,1);
    Hv=Hv*h;

    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-3 -2 -1 0 1 2 3 4 5];
    q_stencil = [-1/168 +1/14 -1/2 -9/20 +5/4 -1/2 +1/6 -1/28 +1/280];
    Qp = stripeMatrixPeriodic(q_stencil, q_diags,m);


    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp) ;

    Dm=HI*(Qm) ;
end
