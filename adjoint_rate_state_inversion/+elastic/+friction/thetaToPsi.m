% Converts between state variables theta (units of time) and psi (dimensionless)
function psi = thetaToPsi(theta, b, V0, f0, L)

	psi = f0 + b.*log(V0.*theta./L);

end