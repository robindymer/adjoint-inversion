% Converts between state variables theta (units of time) and psi (dimensionless)
function theta = psiToTheta(psi, b, V0, f0, L)

	theta = (L./V0).*exp((psi-f0)./b);

end