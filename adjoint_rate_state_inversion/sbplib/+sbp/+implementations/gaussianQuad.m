% GAUSSIANQUAD computes the nth order Gaussian quadrature points and weights
% [r,w] = gaussianQuad(n)
%
% inputs:
%   int n: quadrature order
%
% outputs:
%   doulbe vector r: quadrature points
%   double vector w: quadrature weights
%
% Based on routines given in
%   @BOOK{HesthavenWarburton2008,
%     title = {Nodal Discontinuous {G}alerkin Methods: {A}lgorithms, Analysis, and
%    Applications},
%     publisher = {Springer},
%     year = {2008},
%     author = {Hesthaven, Jan S. and Warburton, Tim},
%     volume = {54},
%     series = {Texts in Applied Mathematics},
%     doi = {10.1007/978-0-387-72067-8}
%   }
function [r,w] = gaussianQuad(n)
  if(n==0)
    r = 0;
    w = 2;
  else
    h = 2*(0:n-1);
    Je = 2./(h + 2).*(1:n).^2./sqrt((h+1).*(h+3));
    A = spdiags([Je,0;0,Je]',[-1 1],n+1,n+1);
    [V,D] = eig(full(A));
    r = diag(D);
    w = 2*(V(1,:)').^2;
  end
end
