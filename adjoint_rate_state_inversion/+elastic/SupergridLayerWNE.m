classdef SupergridLayerWNE < multiblock.DefCurvilinear
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
        function obj = SupergridLayerWNE(pW, pN, pE, h0xW, h0yN, h0xE, m, ratio)
            default_arg('m',51);
            default_arg('ratio', 2^6);

            nBlocksW = length(pW)-1;
            nBlocksN = length(pN)-1;
            nBlocksE = length(pE)-1;
            nBlocksNW = 1;
            nBlocksNE = 1;
            nBlocksVec = [nBlocksW, nBlocksNW, nBlocksN, nBlocksNE, nBlocksE];
            nBlocksCumulative = cumsum(nBlocksVec);
            nBlocks = sum(nBlocksVec);


            % Ensure same grid size in start of layer as in DOI, even with different m:s
            measureDOIWx = h0xW*(m-1);
            measureDOINy = h0yN*(m-1);
            measureDOIEx = h0xE*(m-1);

            SupergridDefW = elastic.SupergridLayerW(pW, h0xW, m, ratio);
            SupergridDefN = elastic.SupergridLayerN(pN, h0yN, m, ratio);
            SupergridDefE = elastic.SupergridLayerE(pE, h0xE, m, ratio);
            SupergridDefNW = elastic.SupergridLayerNW(pN{1}, h0xW, h0yN, m, ratio);
            SupergridDefNE = elastic.SupergridLayerNE(pN{end}, h0xE, h0yN, m, ratio);

            % Blocks
            bW = SupergridDefW.blocks;
            bNW = SupergridDefNW.blocks;
            bN = SupergridDefN.blocks;
            bNE = SupergridDefNE.blocks;
            bE = SupergridDefE.blocks;
            blocks = {bW{:}, bNW{:}, bN{:}, bNE{:}, bE{:}};

            %------ Connections --------------
            conn = cell(nBlocks, nBlocks);

            % West layer
            for i = 1:nBlocksW-1
                conn{i,i+1} = {'n','s'};
            end

            % West to Northwest
            conn{nBlocksW, nBlocksW+1} = {'n','s'};

            % Northwest to North
            conn{nBlocksW+1, nBlocksW+2} = {'e','w'};

            % North layer
            nStart = nBlocksCumulative(2);
            for i = 1:nBlocksN-1
                conn{nStart+i,nStart+i+1} = {'e','w'};
            end

            % North to Northeast
            conn{nStart+nBlocksN,nStart+nBlocksN+1} = {'e','w'};

            % Northeast to East
            nStart = nBlocksCumulative(4);
            conn{nStart,nBlocks} = {'s','n'};

            % East layer
            nStart = nBlocksCumulative(4);
            for i = 1:nBlocksE-1
                conn{nStart+i,nStart+i+1} = {'n','s'};
            end
            % --------------------------------

            %------ Boundary groups ----------
            boundaryGroups = struct;

            W = cell(1, nBlocksW+1);
            for i = 1:nBlocksW+1
                W{i} = {i,'w'};
            end
            boundaryGroups.W = multiblock.BoundaryGroup(W);

            N = cell(1, nBlocksN+2);
            nStart = nBlocksCumulative(1);
            for i = 1:nBlocksN+2
                N{i} = {nStart+i,'n'};
            end
            boundaryGroups.N = multiblock.BoundaryGroup(N);

            E = cell(1, nBlocksE+1);
            for i = 1:nBlocksE+1
                E{i} = {nBlocks+1-i,'e'};
            end
            boundaryGroups.E = multiblock.BoundaryGroup(E);

            S = {{1,'s'}, {nBlocksCumulative(4)+1,'s'}};
            boundaryGroups.S = multiblock.BoundaryGroup(S);

            boundaryGroups.farfield = multiblock.BoundaryGroup([W,N,E]);
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
            obj.SupergridDefN = SupergridDefN;
            obj.SupergridDefE = SupergridDefE;
            obj.SupergridDefNW = SupergridDefNW;
            obj.SupergridDefNE = SupergridDefNE;
        end

        % mTangential -- vector of number of points in tangential direction, for each block.
        function ms = getGridSizes(obj, mTangentialW, mTangentialN, mTangentialE)

            msW = obj.SupergridDefW.getGridSizes(mTangentialW);
            msN = obj.SupergridDefN.getGridSizes(mTangentialN);
            msE = obj.SupergridDefE.getGridSizes(mTangentialE);

            msNW = obj.SupergridDefNW.getGridSizes([]);
            msNE = obj.SupergridDefNE.getGridSizes([]);

            ms = [msW, msNW, msN, msNE, msE];
        end

    end
end
