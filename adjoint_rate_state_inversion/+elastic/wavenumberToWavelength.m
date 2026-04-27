function wvl = wavenumberToWavelength(k)

	% Norm of wavenumber vector
	k = norm(k);

	wvl = 2*pi/k;

end