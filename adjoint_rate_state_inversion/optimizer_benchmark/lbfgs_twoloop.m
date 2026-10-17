% Two-loop recursion: returns H*q where H approximates the inverse Hessian from pairs (S, Y).
% Pairs with y's <= 0 are skipped.
function r = lbfgs_twoloop(q, S, Y)
    keep = sum(S.*Y, 1) > 0;
    S = S(:, keep); Y = Y(:, keep);
    L = size(S, 2);
    if L == 0
        r = q; return
    end
    rho = 1./sum(S.*Y, 1);
    a = zeros(L, 1);
    for i = L:-1:1
        a(i) = rho(i)*(S(:,i)'*q);
        q = q - a(i)*Y(:,i);
    end
    gamma = (S(:,L)'*Y(:,L))/(Y(:,L)'*Y(:,L));
    r = gamma*q;
    for i = 1:L
        b = rho(i)*(Y(:,i)'*r);
        r = r + (a(i) - b)*S(:,i);
    end
end
