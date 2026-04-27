% Computes the grid points x and grid spacing h used by the boundary optimized SBP operators 
% with improved boundary accuracy, presented in 
% 'Boundary optimized diagonal-norm SBP operators - Mattsson, Almquist, van der Weide 2018'.
%
% lim - cell array with domain limits
% N - Number of grid points
% order - order of accuracy of sbp operator.
function [x,h] = upwindAccurateBoundaryOptimizedGrid(lim,N,order)
    assert(iscell(lim) && numel(lim) == 2,'The limit should be cell array with 2 elements.');
    L = lim{2} - lim{1};
    assert(L>0,'Limits must be given in increasing order.');
    %%%% Non-equidistant grid points %%%%%
    xb = boundaryPoints(order);
    m = length(xb)-1; % Number of non-equidistant points
    assert(N-2*(m+1)>=0,'Not enough grid points to contain the boundary region. Requires at least %d points.',2*(m+1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Compute h %%%%%%%%%%
    h = L/(2*xb(end) + N-1-2*m);
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Define grid %%%%%%%%
    x = h*[xb; linspace(xb(end)+1,L/h-xb(end)-1,N-2*(m+1))'; L/h-flip(xb) ];
    x = x + lim{1};
    %%%%%%%%%%%%%%%%%%%%%%%%%
end
function xb = boundaryPoints(order)
    switch order
        case 2
            x0 =  0.0000000000000e+00;
            x1 =  1.0158066888926694105e+00;
            x2 =  2.0158066888926694105e+00;
            xb = [x0 x1 x2]';
        case 3
            x0 =  0.0000000000000e+00;
            x1 =  0.78643399666738163986;
            x2 =  1.78643399666738163986;
            xb = [x0 x1 x2]';
        case 4
            x0 =  0.0000000000000e+00;
            x1 =  0.67776056835278144950;
            x2 =  1.67776056835278144950;
            xb = [x0 x1 x2]';
        case 5
            x0 = 0.0000000000000e+00;
            x1 = 0.62048159180300444152;
            x2 = 1.62048159180300444152;
            xb = [x0 x1 x2]';
        case 6
            x0 =  0.0000000000000e+00;
            x1 =  0.46488458913890803995;
            x2 =  x1 + 0.87082799992290223798;
            xb = [x0 x1 x2]';
        case 7
            x0 =  0.0000000000000e+00;
            x1 =  0.47968913918690932952;
            x2 =  x1 + 0.89563398599582400787;
            xb = [x0 x1 x2]';
        case 8
            x0 =  0.0000000000000e+00;
            x1 =  0.46931054871854486577;
            x2 =  x1 + 0.90168072751351422227;
            xb = [x0 x1 x2]';
        case 9
            x0 =  0.0000000000000e+00;
            x1 =  0.46457087199875349857;
            x2 =  x1 + 0.89474171339471133859;
            xb = [x0 x1 x2]';
        otherwise
            error('Invalid operator order %d.',order);
    end
end