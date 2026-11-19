% Number of forward and adjoint solves of post.prob, minus nFwd0 and nAdj0 if given. Structs
% without an adjoint counter (as in the Gaussian tests) count as 0 adjoint solves.
function [nFwd, nAdj] = solveCounts(post, nFwd0, nAdj0)
    if nargin < 2, nFwd0 = 0; nAdj0 = 0; end
    nFwd = post.prob.nFwd - nFwd0;
    nAdj = -nAdj0;
    if isprop(post.prob, 'nAdj') || isfield(post.prob, 'nAdj')
        nAdj = post.prob.nAdj - nAdj0;
    end
end
