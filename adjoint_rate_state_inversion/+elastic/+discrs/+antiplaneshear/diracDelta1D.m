function dirac_deltas = diracDelta1D(sources, mbGrid, mbDiffOp, moment_conditions)
    dirac_deltas = [];
    smoothness_conditions = 0;
    if ~isempty(sources)
        ns = length(sources.x);
        dirac_deltas = cell(ns,1);
        % Create dirac delta functions
        for i = 1:ns
            x_i  = sources.x(i);
            blockId = sources.blockIds(i);
            grd = mbGrid.grids{blockId};
            H = mbDiffOp.diffOps{i}.H;
            delta_fun_local = diracDiscr(grd, x_i, moment_conditions, smoothness_conditions, H);
            delta_fun = mbGrid.expandFunc(delta_fun_local, blockId);
            dirac_deltas{i} = sparse(delta_fun);
        end
    end
end