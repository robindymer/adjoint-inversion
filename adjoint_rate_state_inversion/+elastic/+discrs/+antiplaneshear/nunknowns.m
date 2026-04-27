function [n, ns, nsip, nsim] = nunknowns(domain, mbGrid, mbDiffOp, bc_method, ic_method)
    
    % ---- Determine total number of unknowns -------%
    n = mbGrid.nPoints;
    switch bc_method
    case 'standard' % No additional variables should be tracked
        ns = 0; 
    case 'erickson2022' % hardcoded to be imposed on both boundaries if used
        ns = 2;
    end

    switch ic_method
    case 'standard'  % No additional variables should be tracked
        nsip = 0;
        nsim = 0;
    case 'erickson2022' % 2 additional variables (one for each side of the interface)
        bidm = domain.boundaryGroups.interface{1};
        bidp = domain.boundaryGroups.interface{2};
        em = mbDiffOp.getBoundaryOperator('e', bidm);
        ep = mbDiffOp.getBoundaryOperator('e', bidp);
        [~, nsim] = size(em);
        [~, nsip] = size(ep);
    end     

end