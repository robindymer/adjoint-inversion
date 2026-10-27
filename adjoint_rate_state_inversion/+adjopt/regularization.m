% Regularization term R(p) added to the misfit, F(p) + R(p), for the parameter vector p.
% Returns a struct with handles
%   reg.R(p)          -- value of R
%   reg.grad(p)       -- gradient of R, column vector
%   reg.hessVec(p, V) -- Hessian of R at p applied to each column of V
%
% opts.type  -- 'none'           R = 0
%               'tikhonov'       R = beta/2*sum((p - p0).^2)
%               'tikhonov-grad'  R = beta/2*sum((D*p).^2), penalizes non-smooth p
%               'tv'             R = beta*sum(sqrt((D*p).^2 + eps^2) - eps), total variation
%               'lasso'          R = beta*sum(sqrt((p - p0).^2 + eps^2) - eps), L1 norm of p - p0
% opts.beta  -- regularization weight
% opts.eps   -- smoothing of the absolute value in 'tv' and 'lasso', which makes R twice
%               differentiable. Should be small compared to the differences/deviations that are
%               expected, but a smaller eps gives a stiffer Hessian (beta/eps at zero).
% p0         -- prior parameter vector for 'tikhonov' and 'lasso', e.g. the initial guess
% blockSizes -- number of entries in p of each inversion parameter. D is the forward difference
%               within each block, so that no differences are taken across parameters.
%
% The sums are over parameter grid points, without grid spacing, so beta depends on m_p.
function reg = regularization(opts, p0, blockSizes)
    default_field(opts, 'type', 'none');
    default_field(opts, 'beta', 0);
    default_field(opts, 'eps', 1e-4);
    beta = opts.beta;
    eps = opts.eps;

    % Forward differences within each block
    Dblocks = cell(numel(blockSizes), 1);
    for i = 1:numel(blockSizes)
        n = blockSizes(i);
        Dblocks{i} = spdiags([-ones(n-1,1), ones(n-1,1)], [0, 1], n-1, n);
    end
    D = blkdiag(Dblocks{:});

    switch opts.type
        case 'none'
            reg.R = @(p) 0;
            reg.grad = @(p) zeros(size(p));
            reg.hessVec = @(p, V) zeros(size(V));
        case 'tikhonov'
            reg.R = @(p) beta/2*sum((p - p0).^2);
            reg.grad = @(p) beta*(p - p0);
            reg.hessVec = @(p, V) beta*V;
        case 'tikhonov-grad'
            reg.R = @(p) beta/2*sum((D*p).^2);
            reg.grad = @(p) beta*(D'*(D*p));
            reg.hessVec = @(p, V) beta*(D'*(D*V));
        case 'tv'
            reg.R = @(p) beta*sum(sqrt((D*p).^2 + eps^2) - eps);
            reg.grad = @(p) beta*(D'*((D*p)./sqrt((D*p).^2 + eps^2)));
            reg.hessVec = @(p, V) beta*(D'*((eps^2./((D*p).^2 + eps^2).^(3/2)).*(D*V)));
        case 'lasso'
            reg.R = @(p) beta*sum(sqrt((p - p0).^2 + eps^2) - eps);
            reg.grad = @(p) beta*(p - p0)./sqrt((p - p0).^2 + eps^2);
            reg.hessVec = @(p, V) beta*(eps^2./((p - p0).^2 + eps^2).^(3/2)).*V;
        otherwise
            error('Unknown regularization type: %s', opts.type);
    end
end
