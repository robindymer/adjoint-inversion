classdef SupergridLayerSW < multiblock.DefCurvilinear
    properties
        p
        m
        ratio
        blocks
    end

    methods
        function obj = SupergridLayerSW(p, h0x, h0y, m, ratio)
            default_arg('m',51);
            default_arg('ratio', 2^6);

            % Same grid size in start of layer as in DOI, even with different m:s
            measureDOIx = h0x*(m-1);
            measureDOIy = h0y*(m-1);

            % Structs for supergrid parameters
            xPars = struct;
            yPars = struct;

            xPars.side = 'l';
            yPars.side = 'l';
            xPars.measuresAtEnds = {measureDOIx, measureDOIx*ratio};
            yPars.measuresAtEnds = {measureDOIy, measureDOIy*ratio};
            xPars.x0 = p(1);
            yPars.x0 = p(2);

            % Create block
            blocks = {elastic.superGridBlock(xPars, yPars)};

            blockNames = {'B1'};
            conn = cell(1,1);

            boundaryGroups = struct();
            boundaryGroups.farfield = multiblock.BoundaryGroup({{1,'w'}, {1,'s'}});

            boundaryGroups.W = multiblock.BoundaryGroup({{1,'w'}});
            boundaryGroups.S = multiblock.BoundaryGroup({{1,'s'}});

            obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blockNames);

            % Save properties
            obj.p = p;
            obj.m = m;
            obj.ratio = ratio;
            obj.blocks = blocks;
        end

        % m is a dummy here
        function ms = getGridSizes(obj, m)
            ms = {[obj.m, obj.m]};
        end

    end
end
