function testDiscrMatrix(order, m, discr_fun, discr, curvilinear)

	default_arg('order',2);
	default_arg('m',[21,31]);
	default_arg('discr_fun',@elasticDiscr);
    default_arg('discr',[]);
    default_arg('curvilinear',false);

    if isempty(discr)
        discr = discr_fun(m, order);
    end

	D = discr.D;
	H = discr.H;
    try
        % Works for single block discr
        RHO = discr.diffOp.RHO;
        RHO_kron = kron(RHO,speye(2));

        if curvilinear
            J = discr.diffOp.J;
            RHO_kron = kron(RHO*J,speye(2));
        end
    catch
        % Multiblock discr
        RHO_kron = discr.RHO_kron;
    end

    M = RHO_kron*H*D;

    SR = abs(eigs(D, 1, 'largestabs'));
    computed = true;
    if isnan(SR)
        disp('Will have to compute spectral radius using eig, which may be very slow.')
        answer = input('Continue(y) or abort(n)?', 's');
        if strcmp(answer, 'y')
            disp('Computing spectral radius using eig');
            SR = max(abs(eig(full(D))));
        else
            computed = false;
        end
    end

    if computed
        fprintf("Spectral radius: %.4e \n", SR);
    end
    helpers.checkSymmetryAndEigenvalues(M);

end

