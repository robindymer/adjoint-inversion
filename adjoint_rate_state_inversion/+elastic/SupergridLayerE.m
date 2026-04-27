classdef SupergridLayerE < multiblock.DefCurvilinear
    properties
        p
        m
        ratio
        blocks
    end

    methods
        function obj = SupergridLayerE(p, h0, m, ratio, curves)
            default_arg('m',51);
            default_arg('ratio', 2^6);
            nBlocks = length(p)-1;
            default_arg('curves', []);

            % Same grid size in start of layer as in DOI, even with different m:s
            measureDOI = h0*(m-1);

            % Structs for supergrid parameters
            xPars = struct;
            yPars = struct;

            xPars.side = 'r';
            yPars.side = 'r';
            xPars.measuresAtEnds = {measureDOI, measureDOI*ratio};

            % Create blocks
            blocks = cell(nBlocks, 1);
            for i = 1:nBlocks
                xPars.x0 = p{i}(1);
                yPars.x0 = p{i}(2);

                measureTangential = p{i+1}(2) - p{i}(2);
                yPars.measuresAtEnds = {measureTangential, measureTangential};

                if isempty(curves)
                    blocks{i} = elastic.superGridBlock(xPars, yPars);
                else
                    lines = struct;
                    lines.L = curves{i}.reverse();
                    blocks{i} = elastic.superGridBlock(xPars, yPars, lines);
                end
            end

            blockNames = cell(nBlocks, 1);
            for i = 1:nBlocks
                blockNames{i} = sprintf('B%d',i);
            end

            conn = cell(nBlocks,nBlocks);
            for i = 1:nBlocks-1
                conn{i,i+1} = {'n', 's'};
            end

            boundaryGroups = struct();

            farfieldBG = cell(nBlocks,1);
            for i = 1:nBlocks
                farfieldBG{i} = {i,'e'};
            end

            boundaryGroups.farfield = multiblock.BoundaryGroup(farfieldBG);

            boundaryGroups.E = multiblock.BoundaryGroup(farfieldBG);
            boundaryGroups.S = multiblock.BoundaryGroup({{1,'s'}});
            boundaryGroups.N = multiblock.BoundaryGroup({{nBlocks,'n'}});

            obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blockNames);

            % Save properties
            obj.p = p;
            obj.m = m;
            obj.ratio = ratio;
            obj.blocks = blocks;
        end

        % mTangential -- vector of number of points in tangential direction, for each block.
        function ms = getGridSizes(obj, mTangential)

            % If mTangential is a scalar, use that number of grid points for all blocks
            if length(mTangential) == 1 && obj.nBlocks > 1
                mTangential = mTangential*ones(obj.nBlocks, 1);
            end

            ms = cell(1, obj.nBlocks);
            for i = 1:obj.nBlocks
                ms{i}(2) = mTangential(i);
                ms{i}(1) = obj.m;
            end
        end

    end
end
