classdef SupergridLayerWNES < multiblock.DefCurvilinear
    properties
        m
        ratio
        SupergridDefW
        SupergridDefNW
        SupergridDefN
        SupergridDefNE
        SupergridDefE
        SupergridDefSE
        SupergridDefS
        SupergridDefSW
    end

    methods

        % If curves are supplied, curves.W should be a cell array of tangential curves for each block.
        % Same for curves.E, etc.
        function obj = SupergridLayerWNES(pW, pN, pE, pS, h0xW, h0yN, h0xE, h0yS, m, ratio, curves)
            default_arg('m',51);
            default_arg('ratio', 2^6);

            curves_default = struct;
            curves_default.W = [];
            curves_default.E = [];
            curves_default.S = [];
            curves_default.N = [];
            default_struct('curves', curves_default);

            nBlocksW = length(pW)-1;
            nBlocksN = length(pN)-1;
            nBlocksE = length(pE)-1;
            nBlocksS = length(pS)-1;
            nBlocksNW = 1;
            nBlocksNE = 1;
            nBlocksSW = 1;
            nBlocksSE = 1;
            nBlocksVec = [nBlocksW, nBlocksNW, nBlocksN, nBlocksNE, nBlocksE, ...
                          nBlocksSE, nBlocksS, nBlocksSW];
            nBlocksCumulative = cumsum(nBlocksVec);
            nBlocks = sum(nBlocksVec);

            % Ensure same grid size in start of layer as in DOI, even with different m:s
            measureDOIWx = h0xW*(m-1);
            measureDOINy = h0yN*(m-1);
            measureDOIEx = h0xE*(m-1);
            measureDOISy = h0yS*(m-1);

            SupergridDefW = elastic.SupergridLayerW(pW, h0xW, m, ratio, curves.W);
            SupergridDefN = elastic.SupergridLayerN(pN, h0yN, m, ratio, curves.N);
            SupergridDefE = elastic.SupergridLayerE(pE, h0xE, m, ratio, curves.E);
            SupergridDefS = elastic.SupergridLayerS(pS, h0yS, m, ratio, curves.S);
            SupergridDefNW = elastic.SupergridLayerNW(pN{1}, h0xW, h0yN, m, ratio);
            SupergridDefNE = elastic.SupergridLayerNE(pN{end}, h0xE, h0yN, m, ratio);
            SupergridDefSW = elastic.SupergridLayerSW(pS{1}, h0xW, h0yS, m, ratio);
            SupergridDefSE = elastic.SupergridLayerSE(pS{end}, h0xE, h0yS, m, ratio);

            % Blocks
            bW = SupergridDefW.blocks;
            bNW = SupergridDefNW.blocks;
            bN = SupergridDefN.blocks;
            bNE = SupergridDefNE.blocks;
            bE = SupergridDefE.blocks;
            bSE = SupergridDefSE.blocks;
            bS = SupergridDefS.blocks;
            bSW = SupergridDefSW.blocks;
            blocks = {bW{:}, bNW{:}, bN{:}, bNE{:}, bE{:}, bSE{:}, bS{:}, bSW{:}};

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
            nNE = nBlocksCumulative(4);
            nE = nBlocksCumulative(5);
            conn{nNE,nE} = {'s','n'};

            % East layer
            nStart = nBlocksCumulative(4);
            for i = 1:nBlocksE-1
                conn{nStart+i,nStart+i+1} = {'n','s'};
            end

            % East to Southeast
            nE = nBlocksCumulative(4) + 1;
            nSE = nE + nBlocksE;
            conn{nE, nSE} = {'s', 'n'};

            % Southeast to South
            conn{nSE, nBlocks-1} = {'w', 'e'};

            % South layer
            nStart = nBlocksCumulative(6);
            for i = 1:nBlocksS-1
                conn{nStart+i,nStart+i+1} = {'e','w'};
            end

            % South to Southwest
            nS = nBlocksCumulative(6) + 1;
            nSW = nBlocks;
            conn{nS, nSW} = {'w', 'e'};

            % West to Southwest
            nW = 1;
            nSW = nBlocks;
            conn{nW, nSW} = {'s', 'n'};
            % --------------------------------

            %------ Boundary groups ----------
            boundaryGroups = struct;

            W = cell(1, nBlocksW+2);
            nStart = 0;
            W{1} = {nBlocks, 'w'};
            for i = 1:nBlocksW+1
                W{i+1} = {nStart+i,'w'};
            end
            boundaryGroups.W = multiblock.BoundaryGroup(W);

            N = cell(1, nBlocksN+2);
            nStart = nBlocksCumulative(1);
            for i = 1:nBlocksN+2
                N{i} = {nStart+i,'n'};
            end
            boundaryGroups.N = multiblock.BoundaryGroup(N);

            E = cell(1, nBlocksE+2);
            nStart = nBlocksCumulative(3);
            for i = 1:nBlocksE+2
                E{i} = {nStart+i,'e'};
            end
            boundaryGroups.E = multiblock.BoundaryGroup(E);

            S = cell(1, nBlocksS+2);
            nStart = nBlocksCumulative(5);
            for i = 1:nBlocksS+2
                S{i} = {nStart+i,'s'};
            end
            boundaryGroups.S = multiblock.BoundaryGroup(S);

            boundaryGroups.farfield = multiblock.BoundaryGroup([W,N,E,S]);
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
            obj.SupergridDefS = SupergridDefS;
            obj.SupergridDefNW = SupergridDefNW;
            obj.SupergridDefNE = SupergridDefNE;
            obj.SupergridDefSW = SupergridDefSW;
            obj.SupergridDefSE = SupergridDefSE;
        end

        % mTangential -- vector of number of points in tangential direction, for each block.
        function ms = getGridSizes(obj, mTangentialW, mTangentialN, mTangentialE, mTangentialS)

            msW = obj.SupergridDefW.getGridSizes(mTangentialW);
            msN = obj.SupergridDefN.getGridSizes(mTangentialN);
            msE = obj.SupergridDefE.getGridSizes(mTangentialE);
            msS = obj.SupergridDefS.getGridSizes(mTangentialS);

            msNW = obj.SupergridDefNW.getGridSizes([]);
            msNE = obj.SupergridDefNE.getGridSizes([]);
            msSW = obj.SupergridDefSW.getGridSizes([]);
            msSE = obj.SupergridDefSE.getGridSizes([]);

            ms = [msW, msNW, msN, msNE, msE, msSE, msS, msSW];
        end

    end
end
