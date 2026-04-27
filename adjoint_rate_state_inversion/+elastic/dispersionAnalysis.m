function [vs, vp, symbols] = dispersionAnalysis(kx, ky, lambda, mu, order, type)

	default_arg('lambda', 1);
	default_arg('mu', 1);
	default_arg('type', 'narrow');
	default_arg('order', 2);

	switch order
	case 2
		P1 = @(k) 1i*sin(k);
		P1_staggered = @(k) 2*P1(k/2);
		P2_narrow = @(k) -4*(sin(k/2)).^2;
		P2_staggered = P2_narrow;

		P1_plus = @(k) cos(k) +1i*sin(k) - 1;
		P1_minus = @(k) -cos(k) +1i*sin(k) + 1;

		P1_remez = P1;
	case 4
		P1 = @(k) 1i*sin(k).*(1 + 2/3*sin(k/2).^2 ) ;
		P2_narrow = @(k) -2/12*cos(2*k) + 2*4/3*cos(k) - 5/2;

		% Staggered D1 = [1/24 -9/8 9/8 -1/24]
		P1_staggered = @(k) 2*1/24*1i*sin(3*k/2) -2*9/8*1i*sin(k/2);
		P2_staggered = @(k) P1_staggered(k).^2;

		P1_remez = P1;

		d4 = [1 -4 6 -4 1];
		P4 = @(k) 2*(d4(1)*cos(2*k) + d4(2)*cos(k)) + d4(3);

		P1_plus = @(k) P1(k) - 1/12*P4(k);
		P1_minus = @(k) P1(k) + 1/12*P4(k);
	case 6
		P1 = @(k) 1i*sin(k).*(1 + 2/3*sin(k/2).^2 + 8/15*sin(k/2).^4) ;

		% D2 stencil: D2 = 1/h^2* [1/90 −3/20 3/2 −49/18 3/2 −3/20 1/90]
		P2_narrow = @(k) 2*(1/90*cos(3*k) -3/20*cos(2*k) + 3/2*cos(k)) -49/18;

		% Staggered D1 = [-0.3e1 / 0.640e3 0.25e2 / 0.384e3 -0.75e2 / 0.64e2 0.75e2 / 0.64e2 -0.25e2 / 0.384e3 0.3e1 / 0.640e3]
		d1 = [-0.3e1 / 0.640e3 0.25e2 / 0.384e3 -0.75e2 / 0.64e2 0.75e2 / 0.64e2 -0.25e2 / 0.384e3 0.3e1 / 0.640e3];
		P1_staggered = @(k) 2*1i*(d1(1)*sin(5/2*k) + d1(2)*sin(3/2*k) + d1(3)*sin(k/2));
		P2_staggered = @(k) P1_staggered(k).^2;

		% Remez 27 DRP
		d1 = [0.78028389, -0.17585010, 0.02380544];
		P1_remez = @(k) 2*1i*(d1(1)*sin(k) + d1(2)*sin(2*k) + d1(3)*sin(3*k));

		d6 = [1 -6 15 -20 15 -6 1];
		P6 = @(k) 2*(d6(1)*cos(3*k) + d6(2)*cos(2*k) + d6(3)*cos(k)) + d6(4);

		P1_plus = @(k) P1(k) + 1/60*P6(k);
		P1_minus = @(k) P1(k) - 1/60*P6(k);
	end

	symbols.D1 = P1;
	symbols.D2 = P2_narrow;
	symbols.D1_staggered = P1_staggered;
	symbols.D1_remez = P1_remez;
	symbols.D1_plus = P1_plus;
	symbols.D1_minus = P1_minus;

	switch type
		case 'narrow'
			P1_l = P1;
			P1_mu = P1;
			P2_l = P2_narrow;
			P2_mu = P2_narrow;
		case 'mixed'
			P1_l = P1;
			P1_mu = P1;
			P2_l = @(k) P1(k).^2;
			P2_mu = P2_narrow;
		case 'D1D1'
			P1_l = P1;
			P1_mu = P1;
			P2_l = @(k) P1(k).^2;
			P2_mu = @(k) P1(k).^2;
		case 'staggered'
			P1_l = P1_staggered;
			P1_mu = P1_staggered;
			P2_l = P2_staggered;
			P2_mu = P2_staggered;
		case 'exact'
			P1_l = @(k) 1i*k;
			P1_mu = @(k) 1i*k;
			P2_l = @(k) -k.^2;
			P2_mu = @(k) -k.^2;
		case 'remez'
			P1_l = P1_remez;
			P1_mu = P1_remez;
			P2_l = @(k) P1_remez(k).^2;
			P2_mu = @(k) P1_remez(k).^2;
		case 'upwind'
			P1p_l = P1_plus;
			P1m_l = P1_minus;
			P1p_mu = P1_plus;
			P1m_mu = P1_minus;
			P2_l = @(k) P1_plus(k).*P1_minus(k);
			P2_mu = @(k) P1_plus(k).*P1_minus(k);
		case 'mixed-upwind'
			P1p_l = P1_plus;
			P1m_l = P1_minus;
			P1p_mu = P1;
			P1m_mu = P1;
			P2_l = @(k) P1_plus(k).*P1_minus(k);
			P2_mu = P2_narrow;
	end

	switch type
	case {'narrow', 'mixed', 'D1D1', 'staggered', 'exact', 'remez'}
		P1p_l = P1_l;
		P1m_l = P1_l;
		P1p_mu = P1_mu;
		P1m_mu = P1_mu;
	end

	P1px_l = P1p_l(kx);
	P1mx_l = P1m_l(kx);
	P1px_mu = P1p_mu(kx);
	P1mx_mu = P1m_mu(kx);

	P1py_l = P1p_l(ky);
	P1my_l = P1m_l(ky);
	P1py_mu = P1p_mu(ky);
	P1my_mu = P1m_mu(ky);

	P2x_l = P2_l(kx);
	P2x_mu = P2_mu(kx);

	P2y_l = P2_l(ky);
	P2y_mu = P2_mu(ky);


	% Q: 2-by-2 spatial operator after Fourier transform
	Q11 = lambda*P2x_l + mu*(2*P2x_mu + P2y_mu);
	Q12 = lambda*P1px_l.*P1my_l + mu*P1px_mu.*P1my_mu;
	Q21 = lambda*P1py_l.*P1mx_l + mu*P1py_mu.*P1mx_mu;
	Q22 = lambda*P2y_l + mu*(2*P2y_mu + P2x_mu);

	% detP = P2(kx).*P2(ky) - (P1(kx).*P1(ky)).^2;

	% Eigenvalues of Q
	% p = -(P2(kx) + P2(ky));
	% q = P2(kx).*P2(ky) - (P1(kx).*P1(ky)).^2;

	p = -(Q11 + Q22);
	q = Q11.*Q22 - Q12.*Q21;

	lm = -p/2 - ((p/2).^2 -  q).^(1/2);
	lp = -p/2 + ((p/2).^2 -  q).^(1/2);

	k2 = kx.^2 + ky.^2 + 1e-10;

	vp = abs( sqrt(-lm./k2) );
	vs = abs( sqrt(-lp./k2) );

end