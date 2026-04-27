% Returns the 1D periodic grid where the last grid point at x_r
% is omitted (it is the same as the the first grid point at x_l)
function [x,h] = get_periodic_grid(x_l,x_r,m)
    L = x_r-x_l;
    h = L/m;
    x = linspace(x_l,x_r,m+1)';
    x = x(1:end-1);
end