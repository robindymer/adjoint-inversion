classdef SupergridLayerWE < multiblock.DefCurvilinear
    properties
        m
        ratio
        SupergridDefW
        SupergridDefNW
        SupergridDefN
        SupergridDefNE
        SupergridDefE
    end

    methods
        % If curves are supplied, curves.W should be a cell array of tangential curves for each block.
        % Same for curves.E.
        function obj = SupergridLayerWE(pW, pE, h0xW, h0xE, m, ratio, curves)
            default_arg('m',51);
            default_arg('ratio', 2^6);

            curves_default = struct;
            curves_default.W = [];
            curves_default.E = [];
            default_struct('curves', curves_default);

            nBlocksW = length(pW)-1;
            nBlocksE = length(pE)-1;
            nBlocksVec = [nBlocksW, nBlocksE];
            nBlocksCumulative = cumsum(nBlocksVec);
            nBlocks = sum(nBlocksVec);

            % Ensure same grid size in start of layer as in DOI, even with different m:s
            measureDOIWx = h0xW*(m-1);
            measureDOIEx = h0xE*(m-1);

            SupergridDefW = elastic.SupergridLayerW(pW, h0xW, m, ratio, curves.W);
            SupergridDefE = elastic.SupergridLayerE(pE, h0xE, m, ratio, curves.E);

            % Blocks
            bW = SupergridDefW.blocks;
            bE = SupergridDefE.blocks;
            blocks = {bW{:}, bE{:}};

            %------ Connections --------------
            conn = cell(nBlocks, nBlocks);

            % West layer
            for i = 1:nBlocksW-1
                conn{i,i+1} = {'n','s'};
            end

            % East layer
            nStart = nBlocksCumulative(1);
            for i = 1:nBlocksE-1
                conn{nStart+i,nStart+i+1} = {'n','s'};
            end
            % --------------------------------

            %------ Boundary groups ----------
            boundaryGroups = struct;

            W = cell(1, nBlocksW);
            for i = 1:nBlocksW
                W{i} = {i,'w'};
            end
            boundaryGroups.W = multiblock.BoundaryGroup(W);

            E = cell(1, nBlocksE);
            for i = 1:nBlocksE
                E{i} = {nBlocks+1-i,'e'};
            end
            boundaryGroups.E = multiblock.BoundaryGroup(E);

            N = {{nBlocksW,'n'}, {nBlocks,'n'}};
            boundaryGroups.N = multiblock.BoundaryGroup(N);

            S = {{1,'s'}, {nBlocksW+1,'s'}};
            boundaryGroups.S = multiblock.BoundaryGroup(S);

            boundaryGroups.farfield = multiblock.BoundaryGroup([W,E]);
            % --------------------------------

            blockNames = cell(nBlocks, 1);
            for i = 1:nBlocks
                blockNames{i} = sprintf('B%d',i);
            end

            obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blockNames);

            % Save properties
            obj.m = m;
            obj.ratio = ratio;

            obj.SupergridDefW = SupergridDefW;
            obj.SupergridDefE = SupergridDefE;
        end

        % mTangential -- vector of number of points in tangential direction, for each block.
        function ms = getGridSizes(obj, mTangentialW, mTangentialE)

            msW = obj.SupergridDefW.getGridSizes(mTangentialW);
            msE = obj.SupergridDefE.getGridSizes(mTangentialE);

            ms = [msW, msE];
        end

    end
end
