% LEGENDREPOLYNOMIALS Computes Legendre polynomials P_i(r)
% of order i from 0 to n using recurrence.
%
% input:
%   double vector r: points to evaluate the polynomials at (-1 <= r <= 1)
%   int n: maximum polynomial order
%   bool do_normalized: true if the polynomials should be normalized.
%
% output:
%   double matrix Ps: Legendre polynomials stored as [P_0, P_1, ..., P_n]
%                     where P_i is the i-th order polynomial stored as a
%                     column vector.
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

function Ps = legendrePolynomials(r,n,do_normalized)
  if do_normalized % Scaling, and recurrence formula from Hesthaven, Warburton
    gamma_i = @(i) 2/(2*i+1);
    a = @(i) sqrt(i^2/((2*i+1)*(2*i-1)));
    recurrence = @(i,P_cur,P_prev) 1/(a(i+1))*(r.*P_cur - a(i)*P_prev);
  else % Bonnet's recurrence formula
    gamma_i = @(i) 1;
    recurrence = @(i,P_cur,P_prev) 1/(i+1)*((2*i+1)*r.*P_cur - i*P_prev);
  end
  Ps = zeros(length(r),n+1);
  r = r(:);
  % 0th polynomial
  if(n >= 0)
    Ps(:,1) = 1/sqrt(gamma_i(0))*ones(size(r));
  end
  % 1st polynomial
  if(n>=1)
    Ps(:,2) = 1/sqrt(gamma_i(1))*r;
  end
  % Use recurrance relation for remaining polynomial orders
  for i=1:n-1
    col = i+1; %column index is offset by 1
    Ps(:,col+1) = recurrence(i,Ps(:,col),Ps(:,col-1));
  end
end
