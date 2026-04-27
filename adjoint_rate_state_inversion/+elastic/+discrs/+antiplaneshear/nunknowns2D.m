function [n, ns, nsip, nsim] = nunknowns2D(domain, mbGrid, mbDiffOp, bc)
    
    % ---- Determine total number of unknowns -------%
    n = mbGrid.nPoints;
    
    % Boundaries without data are imposed using Erickson BCs
    % and requires tracking additional variables.
    ns = 0; 
    for i = 1:length(bc.ids)
        if isempty(bc.data{i})
            bid = bc.ids{i};
            ei = mbDiffOp.getBoundaryOperator('e', bid);
            [~, ni] = size(ei);
            ns = ns + ni;   
        end
    end

    bidm = domain.boundaryGroups.fault_minus;
    bidp = domain.boundaryGroups.fault_plus;
    em = mbDiffOp.getBoundaryOperator('e', bidm);
    ep = mbDiffOp.getBoundaryOperator('e', bidp);
    [~, nsim] = size(em);
    [~, nsip] = size(ep);
end