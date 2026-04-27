function V = asinh_inv(tau, Psi, a, sigma0, V0)
	V = 2*V0./exp(Psi./a).*sinh(tau./(sigma0.*a));