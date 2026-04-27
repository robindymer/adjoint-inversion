function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_optimal_2(m,h)
    
    if(m<4)
        error('Operator requires at least 4 grid points');
    end

    Hv=ones(m,1);
    Hv(1:3)=[0.32580924224615181860e0; 0.11249698536693687180e1; 0.10650275929771488741e1];
    Hv(m-2:m)=rot90(Hv(1:3),2);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [0 1 2];
    q_stencil = [-3/2 +2 -1/2];
    Qp = stripeMatrix(q_stencil, q_diags,m);


    Q_U =[-0.4487186143192082430e-1 0.77254710061771387915e0 -0.22767523918579305475e0; 
          -0.59164110378265189425e0 -0.43233575924045701580e0 0.15239768630231089100e1;
           0.1365129652145727186e0 -0.34021134137725686335e0 -0.12963016238373158552e1;
    ];

    Qp(1:3,1:3)=Q_U;
    Qp(m-2:m,m-2:m)=rot90(Q_U,2)'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    Dm=HI*(Qm-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;
end

