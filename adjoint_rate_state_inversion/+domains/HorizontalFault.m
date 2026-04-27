classdef HorizontalFault < multiblock.DefCurvilinear 

	properties
		xlim
		ylim
	end

	methods
		function obj = HorizontalFault(xlim, ylim)
			default_arg('xlim', [-10, 10]);
			default_arg('ylim', [-10, 0, 10]);

			xl = xlim(1);
			xr = xlim(2);
			yb = ylim(1);
            y0 = ylim(2);
			yt = ylim(3);

            % ---- Block 1 (bot) -------------------------

            % Corners
            NW = [xl; y0];
            NE = [xr; y0];
            SW = [xl; yb];
            SE = [xr; yb];

            W = parametrization.Curve.line(NW, SW);
            E = parametrization.Curve.line(SE, NE);
            S = parametrization.Curve.line(SW, SE);
            N = parametrization.Curve.line(NE, NW);

            B1 = parametrization.Ti(S, E, N, W);

            % ---- Block 2 (top) -------------------------

            % Corners
            NW = [xl; yt];
            NE = [xr; yt];
            SW = [xl; y0];
            SE = [xr; y0];

            W = parametrization.Curve.line(NW, SW);
            E = parametrization.Curve.line(SE, NE);
            S = parametrization.Curve.line(SW, SE);
            N = parametrization.Curve.line(NE, NW);

            B2 = parametrization.Ti(S, E, N, W);
            
            %-----------------------------------------------
            blocks = {B1, B2};
            blockNames = {'\Omega^-', '\Omega^+'};

            conn = cell(2,2);

            boundaryGroups = struct;

            boundaryGroups.surface = multiblock.BoundaryGroup({{2,'n'}});
            boundaryGroups.bottom = multiblock.BoundaryGroup({{1,'s'}});

            boundaryGroups.left = multiblock.BoundaryGroup({{1,'w'},{2,'w'}});
            boundaryGroups.right = multiblock.BoundaryGroup({{1,'e'},{2,'e'}});
            boundaryGroups.sides = multiblock.joinBoundaryGroups(boundaryGroups.left, boundaryGroups.right);

            boundaryGroups.faultBot = multiblock.BoundaryGroup({{1,'n'}});
            boundaryGroups.faultTop = multiblock.BoundaryGroup({{2,'s'}});
            boundaryGroups.fault = multiblock.joinBoundaryGroups(boundaryGroups.faultBot, boundaryGroups.faultTop);

            boundaryGroups.fault_minus = boundaryGroups.faultBot;
            boundaryGroups.fault_plus = boundaryGroups.faultTop;


			obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blockNames);

			obj.xlim = xlim;
			obj.ylim = ylim;
		end

		function ms = getGridSizes(obj, mFault)
            % mFault = number of points along fault

            faultLength = obj.xlim(2) - obj.xlim(1);

            % Heigts in blocks 1 and 2
            H1 = obj.ylim(2) - obj.ylim(1);
            H2 = obj.ylim(3) - obj.ylim(2);

            hFault = faultLength/(mFault - 1);

            % Same spacings in y-directions as on fault
            my1 = round(H1/hFault) + 1;
            my2 = round(H2/hFault) + 1;

            ms = {[mFault, my1], [mFault, my2]};

        end

	end

end
