% Copmutes angular frequencey from wave speed and wavelength
function w = angularFrequency(c, lambda)

	f = c/lambda;
	w = 2*pi*f;

end