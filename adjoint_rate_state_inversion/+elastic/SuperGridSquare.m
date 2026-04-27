classdef SuperGridSquare < multiblock.DefCurvilinear
    properties
        L
        DOI_IDs
        W_IDs, E_IDs, S_IDs, N_IDs
        SG
    end

    methods
        % L     - length of domain of interest
        % SG    - struct of supergrid parameters
        function obj = SuperGridSquare(L, SG)
            default_arg('L', 1);

            defaultSG = struct;
            defaultSG.ratio = 2^6;
            defaultSG.mDOI = 51;
            defaultSG.mSG = 51;
            default_struct('SG', defaultSG);

            xl = -L/2;
            xr = L/2;

            yl = -L/2;
            yr = L/2;

            % Blow-up in domain of interest:
            measureDOI = L;

            % Same grid size in start of layer as in DOI, even with different m:s
            mRatio = (SG.mDOI-1) / (SG.mSG-1);
            measureSG = measureDOI / mRatio;

            % Structs for supergrid parameters
            xPars = struct;
            xPars.x0 = xl;
            xPars.side = 'r';
            xPars.measuresAtEnds = {measureDOI, measureDOI};

            yPars = xPars;
            yPars.x0 = yl;

            % BLock config:

            %   LT ---- CT ---- RT
            %   |       |       |
            %   LC ---- CC ---- RC
            %   |       |       |
            %   LB ---- CB ---- RB

            % Center
            xParsC = xPars;
            yParsC = yPars;

            % Right
            xParsR = xPars;
            xParsR.x0 = xr;
            xParsR.measuresAtEnds = {measureSG, SG.ratio*measureSG};

            yParsR = yPars;
            yParsR.x0 = yr;
            yParsR.measuresAtEnds = {measureSG, SG.ratio*measureSG};

            % Left
            xParsL = xPars;
            xParsL.x0 = xl;
            xParsL.measuresAtEnds = {measureSG, SG.ratio*measureSG};
            xParsL.side = 'l';

            yParsL = yPars;
            yParsL.x0 = yl;
            yParsL.measuresAtEnds = {measureSG, SG.ratio*measureSG};
            yParsL.side = 'l';

            % Create block Ti:s
            LT = elastic.superGridBlock(xParsL, yParsR);
            CT = elastic.superGridBlock(xParsC, yParsR);
            RT = elastic.superGridBlock(xParsR, yParsR);

            LC = elastic.superGridBlock(xParsL, yParsC);
            CC = elastic.superGridBlock(xParsC, yParsC);
            RC = elastic.superGridBlock(xParsR, yParsC);

            LB = elastic.superGridBlock(xParsL, yParsL);
            CB = elastic.superGridBlock(xParsC, yParsL);
            RB = elastic.superGridBlock(xParsR, yParsL);

            blocks = {LT, CT, RT, LC, CC, RC, LB, CB, RB};
            blockNames = {'LT', 'CT', 'RT', 'LC', 'CC', 'RC', 'LB', 'CB', 'RB'};

            conn = cell(9,9);
            conn{1,2} = {'e', 'w'};
            conn{1,4} = {'s', 'n'};
            conn{2,3} = {'e', 'w'};
            conn{2,5} = {'s', 'n'};
            conn{3,6} = {'s', 'n'};

            conn{4,5} = {'e', 'w'};
            conn{4,7} = {'s', 'n'};
            conn{5,6} = {'e', 'w'};
            conn{5,8} = {'s', 'n'};
            conn{6,9} = {'s', 'n'};

            conn{7,8} = {'e', 'w'};
            conn{8,9} = {'e', 'w'};

            boundaryGroups = struct();
            boundaryGroups.E = multiblock.BoundaryGroup({{3,'e'}, {6,'e'}, {9, 'e'}});
            boundaryGroups.W = multiblock.BoundaryGroup({{1,'w'}, {4,'w'}, {7, 'w'}});
            boundaryGroups.N = multiblock.BoundaryGroup({{1,'n'}, {2,'n'}, {3, 'n'}});
            boundaryGroups.S = multiblock.BoundaryGroup({{7,'s'}, {8,'s'}, {9, 's'}});
            boundaryGroups.all = multiblock.BoundaryGroup({ {3,'e'}, {6,'e'}, {9, 'e'},...
                                                            {1,'w'}, {4,'w'}, {7, 'w'},...
                                                            {1,'n'}, {2,'n'}, {3, 'n'},...
                                                            {7,'s'}, {8,'s'}, {9, 's'} });

            obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blockNames);

            % Save properties
            obj.L = L;
            obj.DOI_IDs = [5];
            obj.W_IDs = [1,4,7];
            obj.E_IDs = [3,6,9];
            obj.N_IDs = 1:3;
            obj.S_IDs = 7:9;
            obj.SG = SG;
        end

        function ms = getGridSizes(obj, m)

            mSG = obj.SG.mSG;

            LT = [mSG, mSG];
            CT = [m, mSG];
            RT = [mSG, mSG];

            LC = [mSG, m];
            CC = [m, m];
            RC = [mSG, m];

            LB = [mSG, mSG];
            CB = [m, mSG];
            RB = [mSG, mSG];

            ms = {LT, CT, RT, LC, CC, RC, LB, CB, RB};
        end

    end
end