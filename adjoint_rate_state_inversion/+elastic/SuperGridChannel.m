classdef SuperGridChannel < multiblock.DefCurvilinear
    properties
        L
        DOI_IDs
        W_IDs, E_IDs, S_IDs, N_IDs
        SG
    end

    methods
        % L     - length of domain of interest
        % SG    - struct of supergrid parameters
        function obj = SuperGridChannel(L, SG)
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
            %   L ---- C ---- R

            % Center
            xParsC = xPars;
            yParsC = yPars;

            % Right
            xParsR = xPars;
            xParsR.x0 = xr;
            xParsR.measuresAtEnds = {measureSG, SG.ratio*measureSG};

            % Left
            xParsL = xPars;
            xParsL.x0 = xl;
            xParsL.measuresAtEnds = {measureSG, SG.ratio*measureSG};
            xParsL.side = 'l';


            % Create block Ti:s
            LB = elastic.superGridBlock(xParsL, yParsC);
            CB = elastic.superGridBlock(xParsC, yParsC);
            RB = elastic.superGridBlock(xParsR, yParsC);

            blocks = {LB, CB, RB};
            blockNames = {'L', 'C', 'R'};

            conn = cell(3,3);
            conn{1,2} = {'e', 'w'};
            conn{2,3} = {'e', 'w'};

            boundaryGroups = struct();
            boundaryGroups.E = multiblock.BoundaryGroup({{3,'e'}});
            boundaryGroups.W = multiblock.BoundaryGroup({{1,'w'}});
            boundaryGroups.N = multiblock.BoundaryGroup({{1,'n'}, {2,'n'}, {3, 'n'}});
            boundaryGroups.S = multiblock.BoundaryGroup({{1,'s'}, {2,'s'}, {3, 's'}});
            boundaryGroups.all = multiblock.BoundaryGroup({ {3,'e'},...
                                                            {1,'w'},...
                                                            {1,'n'}, {2,'n'}, {3, 'n'},...
                                                            {1,'s'}, {2,'s'}, {3, 's'} });

            obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blockNames);

            % Save properties
            obj.L = L;
            obj.DOI_IDs = [2];
            obj.W_IDs = [1];
            obj.E_IDs = [3];
            obj.S_IDs = [];
            obj.N_IDs = [];
            obj.SG = SG;
        end

        function ms = getGridSizes(obj, m)

            mSG = obj.SG.mSG;

            L = [mSG, m];
            C = [m, m];
            R = [mSG, m];

            ms = {L, C, R};
        end

    end
end