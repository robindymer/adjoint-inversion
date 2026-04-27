function [H, HI, D1, D2, DI] = d2_noneq_variable_4(N, h, options)
    % N: Number of grid points
    % h: grid spacing
    % options: struct containing options for constructing the operator
    %          current options are: 
    %               options.stencil_type ('minimal','nonminimal','wide')
    %               options.AD ('upwind', 'op')

    % BP: Number of boundary points
    % order: Accuracy of interior stencil
    BP = 4;
    order = 4;
    if(N<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    %%%% Norm matrix %%%%%%%%
    P = zeros(BP, 1);
    P0 = 2.1259737557798e-01;
    P1 = 1.0260290400758e+00;
    P2 = 1.0775123588954e+00;
    P3 = 9.8607273802835e-01;

    for i = 0:BP - 1
        P(i + 1) = eval(['P' num2str(i)]);
    end

    Hv = ones(N, 1);
    Hv(1:BP) = P;
    Hv(end - BP + 1:end) = flip(P);
    Hv = h * Hv;
    H = spdiags(Hv, 0, N, N);
    HI = spdiags(1 ./ Hv, 0, N, N);
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Q matrix %%%%%%%%%%%
    d = [1/12, -2/3, 0, 2/3, -1/12];
    d = repmat(d, N, 1);
    Q = spdiags(d, -order / 2:order / 2, N, N);

    % Boundaries
    Q0_0 = -5.0000000000000e-01;
    Q0_1 = 6.5605279837843e-01;
    Q0_2 = -1.9875859409017e-01;
    Q0_3 = 4.2705795711740e-02;
    Q0_4 = 0.0000000000000e+00;
    Q0_5 = 0.0000000000000e+00;
    Q1_0 = -6.5605279837843e-01;
    Q1_1 = 0.0000000000000e+00;
    Q1_2 = 8.1236966439895e-01;
    Q1_3 = -1.5631686602052e-01;
    Q1_4 = 0.0000000000000e+00;
    Q1_5 = 0.0000000000000e+00;
    Q2_0 = 1.9875859409017e-01;
    Q2_1 = -8.1236966439895e-01;
    Q2_2 = 0.0000000000000e+00;
    Q2_3 = 6.9694440364211e-01;
    Q2_4 = -8.3333333333333e-02;
    Q2_5 = 0.0000000000000e+00;
    Q3_0 = -4.2705795711740e-02;
    Q3_1 = 1.5631686602052e-01;
    Q3_2 = -6.9694440364211e-01;
    Q3_3 = 0.0000000000000e+00;
    Q3_4 = 6.6666666666667e-01;
    Q3_5 = -8.3333333333333e-02;

    for i = 1:BP

        for j = 1:BP
            Q(i, j) = eval(['Q' num2str(i - 1) '_' num2str(j - 1)]);
            Q(N + 1 - i, N + 1 - j) = -eval(['Q' num2str(i - 1) '_' num2str(j - 1)]);
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Undivided difference operators %%%%
    % Closed with zeros at the first boundary nodes.
    m = N;

    DD_2 = (diag(ones(m - 1, 1), -1) - 2 * diag(ones(m, 1), 0) + diag(ones(m - 1, 1), 1));
    DD_2(1:3, 1:4) = [0 0 0 0; 0.16138369498429727170e1 -0.26095138364100825853e1 0.99567688656710986834e0 0; 0 0.84859980956172494512e0 -0.17944203477786665350e1 0.94582053821694158989e0; ];
    DD_2(m - 2:m, m - 3:m) = [0.94582053821694158989e0 -0.17944203477786665350e1 0.84859980956172494512e0 0; 0 0.99567688656710986834e0 -0.26095138364100825853e1 0.16138369498429727170e1; 0 0 0 0; ];
    DD_2 = sparse(DD_2);

    DD_3 = (-diag(ones(m - 2, 1), -2) + 3 * diag(ones(m - 1, 1), -1) - 3 * diag(ones(m, 1), 0) + diag(ones(m - 1, 1), 1));
    DD_3(1:4, 1:5) = [0 0 0 0 0; 0 0 0 0 0; -0.17277463987989539852e1 0.37021976718569105700e1 -0.29870306597013296050e1 0.10125793866433730203e1 0; 0 -0.81738495424057284493e0 0.26916305216679998025e1 -0.28374616146508247697e1 0.96321604722339781208e0; ];
    DD_3(m - 2:m, m - 4:m) = [-0.96321604722339781208e0 0.28374616146508247697e1 -0.26916305216679998025e1 0.81738495424057284493e0 0; 0 -0.10125793866433730203e1 0.29870306597013296050e1 -0.37021976718569105700e1 0.17277463987989539852e1; 0 0 0 0 0; ];
    DD_3 = sparse(DD_3);

    DD_4 = (diag(ones(m - 2, 1), 2) - 4 * diag(ones(m - 1, 1), 1) + 6 * diag(ones(m, 1), 0) - 4 * diag(ones(m - 1, 1), -1) + diag(ones(m - 2, 1), -2));
    DD_4(1:4, 1:6) = [0 0 0 0 0 0; 0 0 0 0 0 0; 0.18176226052481525189e1 -0.47546882767009058782e1 0.59740613194026592100e1 -0.40503175465734920811e1 0.10133218986235862303e1 0; 0 0.79462567299107735362e0 -0.35888406955573330700e1 0.56749232293016495393e1 -0.38528641888935912483e1 0.97215598215819742539e0; ];
    DD_4(m - 3:m, m - 5:m) = [0.97215598215819742539e0 -0.38528641888935912483e1 0.56749232293016495393e1 -0.35888406955573330700e1 0.79462567299107735362e0 0; 0 0.10133218986235862303e1 -0.40503175465734920811e1 0.59740613194026592100e1 -0.47546882767009058782e1 0.18176226052481525189e1; 0 0 0 0 0 0; 0 0 0 0 0 0; ];
    DD_4 = sparse(DD_4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Difference operators %%%
    D1 = H \ Q;

    % Helper functions for constructing D2(c)
    % TODO: Consider changing sparse(diag(...)) to spdiags(....)
    min_inds = sbp.implementations.d2_sparsity_pattern_inds(m, order, BP, 0, 2);
    nonmin_inds = sbp.implementations.d2_sparsity_pattern_inds(m, order, BP, 1, 2);

    % Minimal 5 point stencil width
    function D2 = D2_fun_minimal(c)
        % Here we add variable diffusion
        C1 = sparse(diag(c));
        C2 = 1/2 * diag(ones(m - 1, 1), -1) + 1/2 * diag(ones(m, 1), 0); C2(1, 2) = 1/2;

        C2 = sparse(diag(C2 * c));

        % Remainder term added to wide second drivative opereator, to obtain a 5
        % point narrow stencil.
        R = (1/144 / h) * transpose(DD_4) * C1 * DD_4 + (1/18 / h) * transpose(DD_3) * C2 * DD_3;
        D2 = D1 * C1 * D1 - H \ R;

        % Remove potential round off zeros
        D2tmp = sparse(m,m);
        D2tmp(min_inds) = D2(min_inds);
        D2 = D2tmp;
    end

    %  Non-minimal 7 point stencil width
    function D2 = D2_fun_nonminimal(c)
        % Here we add variable diffusion
        C1 = sparse(diag(c));

        % Remainder term added to wide second derivative operator
        R = (1/144 / h) * transpose(DD_4) * C1 * DD_4;
        D2 = D1 * C1 * D1 - H \ R;
        
        % Remove potential round off zeros
        D2tmp = sparse(m,m);
        D2tmp(nonmin_inds) = D2(nonmin_inds);
        D2 = D2tmp;
    end

    % Wide stencil
    function D2 = D2_fun_wide(c)
        % Here we add variable diffusion
        C1 = sparse(diag(c));
        D2 = D1 * C1 * D1;
    end

    switch options.stencil_width
        case 'minimal'
            D2 = @D2_fun_minimal;
        case 'nonminimal'
            D2 = @D2_fun_nonminimal;
        case 'wide'
            D2 = @D2_fun_wide;
        otherwise
            error('No option %s for stencil width', options.stencil_width)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Artificial dissipation operator %%%
    switch options.AD
        case 'upwind'
            % This is the choice that yield 3rd order Upwind
            DI = H \ (transpose(DD_2) * DD_2) * (-1/12);
        case 'op'
            % This choice will preserve the order of the underlying
            % Non-dissipative D1 SBP operator
            DI = H \ (transpose(DD_3) * DD_3) * (-1 / (5 * 12));
        otherwise
            error("Artificial dissipation options '%s' not implemented.", option.AD)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end