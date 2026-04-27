classdef SupergridLayerS < multiblock.DefCurvilinear
    properties
        p
        m
        ratio
        blocks
    end

    methods
        function obj = SupergridLayerS(p, h0, m, ratio, curves)
            default_arg('m',51);
            default_arg('ratio', 2^6);
            default_arg('curves', []);

            nBlocks = length(p)-1;

            % Same grid size in start of layer as in DOI, even with different m:s
            measureDOI = h0*(m-1);

            % Structs for supergrid parameters
            xPars = struct;
            yPars = struct;

            xPars.side = 'r';
            yPars.side = 'l';
            yPars.measuresAtEnds = {measureDOI, measureDOI*ratio};

            % Create blocks
            blocks = cell(nBlocks, 1);
            for i = 1:nBlocks
                xPars.x0 = p{i}(1);
                yPars.x0 = p{i}(2);

                measureTangential = p{i+1}(1) - p{i}(1);
                xPars.measuresAtEnds = {measureTangential, measureTangential};

                if isempty(curves)
                    blocks{i} = elastic.superGridBlock(xPars, yPars);
                else
                    lines = struct;
                    lines.T = curves{i}.reverse();
                    blocks{i} = elastic.superGridBlock(xPars, yPars, lines);
                end
            end

            blockNames = cell(nBlocks, 1);
            for i = 1:nBlocks
                blockNames{i} = sprintf('B%d',i);
            end

            conn = cell(nBlocks,nBlocks);
            for i = 1:nBlocks-1
                conn{i,i+1} = {'e', 'w'};
            end

            boundaryGroups = struct();

            farfieldBG = cell(nBlocks,1);
            for i = 1:nBlocks
                farfieldBG{i} = {i,'s'};
            end

            boundaryGroups.farfield = multiblock.BoundaryGroup(farfieldBG);

            boundaryGroups.W = multiblock.BoundaryGroup({{1,'w'}});
            boundaryGroups.S = multiblock.BoundaryGroup(farfieldBG);
            boundaryGroups.E = multiblock.BoundaryGroup({{nBlocks,'e'}});

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
                ms{i}(1) = mTangential(i);
                ms{i}(2) = obj.m;
            end
        end

    end
end
