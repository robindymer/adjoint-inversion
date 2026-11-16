% Hessian of the potential U at a, from m_p Hessian-vector products (one per column of the identity):
%
%   H = c*H_M + P,   H_M the Hessian of the misfit (Gauss-Newton + second-order terms), P the prior precision
%
% Also returns the eigenvalues lambda of the prior-preconditioned data Hessian R'^-1 (c*H_M) R^-1,
% sorted descending: a direction with lambda >> 1 is informed by the data, with lambda << 1 by the
% prior. The posterior width in whitened prior coordinates is about 1/sqrt(1 + lambda).
function [H, lambda] = posterior_hessian(post, a)
    m_p = numel(a);
    H = post.Uhv(a, eye(m_p));
    H = (H + H')/2;   % remove round-off asymmetry
    if min(eig(H)) <= 0
        error('The Hessian at a is not positive definite (min eigenvalue %.3e)', min(eig(H)));
    end
    R = post.prior.R;
    lambda = sort(eig((R'\(H - post.prior.P))/R), 'descend');
end
