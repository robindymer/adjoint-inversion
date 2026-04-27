function F = asinh(V, Psi, a, sigma0, V0)
	F = sigma0.*a.*asinh(V.*exp(Psi./a)./(2*V0));