% Computes the max wave speed over all the directions in thetaVec
% in a general anistropic elastic medium
% rho   	- density (vector)
% C     	- stiffness tensor, C{i,j,k,l} (cell tensor of vectors)
% thetaVec  - angle of propagation (vector)
% thetaVec may also be set to 'optimize', in which case an optimization routine is applied
function v_max = maxAnisotropicWaveSpeed(rho, C, thetaVec)

	if isa(thetaVec, 'double')
	    max_vs = zeros(size(thetaVec));
	    for th = 1:length(thetaVec)

	        theta = thetaVec(th);
	        v_max_vec = elastic.anisotropicWaveSpeed(rho, C, theta);
	        max_vs(th) = max(v_max_vec);

	    end
	    v_max = max(max_vs);

	elseif strcmp(thetaVec, 'optimize')

		f = @(theta) -max( elastic.anisotropicWaveSpeed(rho, C, theta) );
		options = optimoptions('fmincon', 'Display','off');

		% Search for theta in [0, pi]
		[~, v_max_neq] = fmincon(f, pi/2, [], [], [], [], 0, pi, [], options);
		v_max = -v_max_neq;
	else
		error('Incorrect thetaVec type');
	end

end