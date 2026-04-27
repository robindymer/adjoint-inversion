function k = wavelengthToWavenumber(wvl, dir)

	% Normalize direction
	dir = dir/norm(dir);

	% Magnitude of wavenumber
	k = 2*pi/wvl;

	% Multiply magnitude by direction to get the wavenumber vector
	k = k*dir;

end