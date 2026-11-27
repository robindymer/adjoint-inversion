% Negative log posterior (potential) for a, from an InvProblem, its noise and a prior:
%
%   U(a) = M(a)/(sigma^2*dt_d) - log p(a)
%
% post.Phi(a)      negative log likelihood M(a)/(sigma^2*dt_d); inf if a <= 0 somewhere or the
%                  forward solve fails (NaN). pCN needs only this, since its proposal keeps the prior.
% post.U(a)        potential Phi(a) - log p(a)
% post.Ug(a)       potential and gradient (forward + adjoint solve)
% post.Uhv(a, V)   Hessian-vector products (Gauss-Newton + second-order terms of the misfit)
% post.c           likelihood scale 1/(sigma^2*dt_d)
function post = posterior(prob, noise, prior)
    c = 1/(noise.sigma^2*noise.dt_d);
    post.c = c;
    post.prob = prob;
    post.noise = noise;
    post.prior = prior;
    post.Phi = @Phi;
    post.U = @(a) Phi(a) - prior.logp(a);
    post.Ug = @Ug;
    if isfield(prior, 'hv')   % non-Gaussian prior, e.g. tv_prior.m
        post.Uhv = @(a, V) c*prob.Hv(a, V) + prior.hv(a, V);
    else
        post.Uhv = @(a, V) c*prob.Hv(a, V) + prior.P*V;
    end

    function phi = Phi(a)
        phi = c*prob.f(a);
        if isnan(phi)
            phi = inf;
        end
    end

    function [u, g] = Ug(a)
        [f, gf] = prob.fg(a);
        u = c*f - prior.logp(a);
        g = c*gf - prior.grad(a);
        if isnan(u) || any(isnan(g))
            u = inf;
            g = nan(size(a));
        end
    end
end
