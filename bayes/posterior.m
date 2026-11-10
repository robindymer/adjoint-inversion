% Negative log posterior (potential) for a, from an InvProblem, its noise and a prior:
%
%   U(a) = M(a)/(sigma^2*dt_d) - log p(a)
%
% post.U(a)        potential; inf if a <= 0 somewhere or the forward solve fails (NaN)
% post.Ug(a)       potential and gradient (forward + adjoint solve)
% post.Uhv(a, V)   Hessian-vector products (Gauss-Newton + second-order terms of the misfit)
% post.c           likelihood scale 1/(sigma^2*dt_d)
function post = posterior(prob, noise, prior)
    c = 1/(noise.sigma^2*noise.dt_d);
    post.c = c;
    post.prob = prob;
    post.noise = noise;
    post.prior = prior;
    post.U = @U;
    post.Ug = @Ug;
    post.Uhv = @(a, V) c*prob.Hv(a, V) + prior.P*V;

    function u = U(a)
        u = c*prob.f(a) - prior.logp(a);
        if isnan(u)
            u = inf;
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
