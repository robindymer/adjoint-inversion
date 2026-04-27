function etaInv = basic(sigma)

	% Compute an invariant
	tauBar = sqrt(sigma{1,2}.^2 + sigma{2,1}.^2);

	n = 3;
	% etaInv = tauBar.^(n-1);
	etaInv = 0.01*tauBar.^(n-1);

end