function residual = nonlin_solve_fun_asinh(V, Psi, eta, tau_l, a, sigma0, V0)
	F = sigma0.*a.*asinh(V.*exp(Psi./a)./(2*V0));
	residual = eta.*V + F + tau_l;