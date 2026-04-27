% Computes the wave speed in the direction given by theta
% in a general anistropic elastic medium
% rho   - density (vector)
% C     - stiffness tensor, C{i,j,k,l} (cell tensor of vectors)
% theta - angle of propagation (scalar)
function [v_max, v_min] = anisotropicWaveSpeed(rho, C, theta)

    dim = 2;
    y = [cos(theta), sin(theta)];

    R = cell(dim,dim);
    for j = 1:dim
        for l = 1:dim
            R{j,l} = 0*C{1,1,1,1};
        end
    end

    % R_jl = y_i C_ijkl y_k
    for i = 1:dim
        for j = 1:dim
            for k = 1:dim
                for l = 1:dim
                    R{j,l} = R{j,l} + y(i)*C{i,j,k,l}*y(k);
                end
            end
        end
    end

    % Solve quadratic for v^2.
    p = R{1,1} + R{2,2};
    q = R{1,1}.*R{2,2} - R{1,2}.*R{2,1};

    v_max_2 = 1./rho .* (p/2 + sqrt( 1/4*p.^2 - q ) );
    v_min_2 = 1./rho .* (p/2 - sqrt( 1/4*p.^2 - q ) );

    v_max = sqrt(v_max_2);
    v_min = sqrt(abs(v_min_2));

end