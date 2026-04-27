classdef SupergridLayerWSE < multiblock.DefCurvilinear
    properties
        m
        ratio
        SupergridDefW
        SupergridDefSW
        SupergridDefS
        SupergridDefSE
        SupergridDefE
    end

    methods
        function obj = SupergridLayerWSE(pW, pS, pE, h0xW, h0yS, h0xE, m, ratio)
            default_arg('m',51);
            default_arg('ratio', 2^6);

            nBlocksW = length(pW)-1;
            nBlocksS = length(pS)-1;
            nBlocksE = length(pE)-1;
            nBlocksSW = 1;
            nBlocksSE = 1;
            nBlocksVec = [nBlocksW, nBlocksSW, nBlocksS, nBlocksSE, nBlocksE];
            nBlocksCumulative = cumsum(nBlocksVec);
            nBlocks = sum(nBlocksVec);


            % Ensure same grid size in start of layer as in DOI, even with different m:s
            measureDOIWx = h0xW*(m-1);
            measureDOISy = h0yS*(m-1);
            measureDOIEx = h0xE*(m-1);

            SupergridDefW = elastic.SupergridLayerW(pW, h0xW, m, ratio);
            SupergridDefS = elastic.SupergridLayerS(pS, h0yS, m, ratio);
            SupergridDefE = elastic.SupergridLayerE(pE, h0xE, m, ratio);
            SupergridDefSW = elastic.SupergridLayerSW(pS{1}, h0xW, h0yS, m, ratio);
            SupergridDefSE = elastic.SupergridLayerSE(pS{end}, h0xE, h0yS, m, ratio);

            % Blocks
            bW = SupergridDefW.blocks;
            bSW = SupergridDefSW.blocks;
            bS = SupergridDefS.blocks;
            bSE = SupergridDefSE.blocks;
            bE = SupergridDefE.blocks;
            blocks = {bW{:}, bSW{:}, bS{:}, bSE{:}, bE{:}};

            %------ Connections --------------
            conn = cell(nBlocks, nBlocks);

            % West layer
            for i = 1:nBlocksW-1
                conn{i,i+1} = {'n','s'};
            end

            % West to Southwest
            conn{1, nBlocksW+1} = {'s','n'};

            % Southwest to South
            conn{nBlocksW+1, nBlocksW+2} = {'e','w'};

            % South layer
            nStart = nBlocksCumulative(2);
            for i = 1:nBlocksS-1
                conn{nStart+i,nStart+i+1} = {'e','w'};
            end

            % South to Southeast
            conn{nStart+nBlocksS,nStart+nBlocksS+1} = {'e','w'};

            % Southeast to East
            nStart = nBlocksCumulative(4);
            conn{nStart,nStart+1} = {'n','s'};

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

            S = cell(1, nBlocksS+2);
            nStart = nBlocksCumulative(1);
            for i = 1:nBlocksS+2
                S{i} = {nStart+i,'s'};
            end
            boundaryGroups.S = multiblock.BoundaryGroup(S);

            E = cell(1, nBlocksE+1);
            for i = 1:nBlocksE+1
                E{i} = {nBlocks+1-i,'e'};
            end
            boundaryGroups.E = multiblock.BoundaryGroup(E);

            N = {{nBlocksW,'n'}, {nBlocks,'n'}};
            boundaryGroups.N = multiblock.BoundaryGroup(N);

            boundaryGroups.farfield = multiblock.BoundaryGroup([W,S,E]);
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
            obj.SupergridDefS = SupergridDefS;
            obj.SupergridDefE = SupergridDefE;
            obj.SupergridDefSW = SupergridDefSW;
            obj.SupergridDefSE = SupergridDefSE;
        end

        % mTangential -- vector of number of points in tangential direction, for each block.
        function ms = getGridSizes(obj, mTangentialW, mTangentialS, mTangentialE)

            msW = obj.SupergridDefW.getGridSizes(mTangentialW);
            msS = obj.SupergridDefS.getGridSizes(mTangentialS);
            msE = obj.SupergridDefE.getGridSizes(mTangentialE);

            msSW = obj.SupergridDefSW.getGridSizes([]);
            msSE = obj.SupergridDefSE.getGridSizes([]);

            ms = [msW, msSW, msS, msSE, msE];
        end

    end
end
