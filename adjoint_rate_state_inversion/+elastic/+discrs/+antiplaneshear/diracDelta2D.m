function dirac_deltas = diracDelta2D(sources, mbGrid, order, opSet)
    dirac_deltas = [];
    moment_conditions = order;
    smoothness_conditions = 0;
    if ~isempty(sources)
        ns = length(sources.x);
        dirac_deltas = sparse(mbGrid.N(),ns);
        % Create dirac delta functions
        for i = 1:ns
            x_i  = sources.x{i};
            blockId = sources.blockIds(i);
            grd = mbGrid.grids{blockId};
            delta_fun_local = diracDiscrCurve(x_i, grd, moment_conditions, smoothness_conditions, order, opSet);
            delta_fun = mbGrid.expandFunc(delta_fun_local, blockId);
            dirac_deltas(:,i) = delta_fun;
        end
    end
end