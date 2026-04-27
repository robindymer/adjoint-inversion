function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_optimal_3(m,h)
    
    if(m<4)
        error('Operator requires at least 4 grid points');
    end

    Hv=ones(m,1);
    Hv(1:3)=[0.27352620305030932792e0; 0.10237033707827460770e1; 0.98920442283432623495e0;];
    Hv(m-2:m)=rot90(Hv(1:3),2);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-1 0 1 2];
    q_stencil = [-0.1e1 / 0.3e1 -0.1e1 / 0.2e1 1 -0.1e1 / 0.6e1;];
    Qp = stripeMatrix(q_stencil, q_diags,m);

    Q_U =[-0.1745959018893113752e-1 0.65088120076476834894e0 -0.13342161057583721138e0;
        -0.58850038980199463513e0 -0.13905293405512346370e0 0.89421999052378476544e0;
        0.10595997999092577264e0 -0.51182826670964488524e0 -0.42746504661461422073e0;
    ];

    Qp(1:3,1:3)=Q_U;
    Qp(m-2:m,m-2:m)=rot90(Q_U,2)'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    Dm=HI*(Qm-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;
end

