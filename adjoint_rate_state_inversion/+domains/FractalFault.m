classdef FractalFault < multiblock.DefCurvilinear 

	properties
		xlim
		ylim
        hx
	end

	methods
		function obj = FractalFault(xlim, ylim, m)
			default_arg('xlim', [-10, 10]);
            default_arg('ylim', [-10, 10]);
            default_arg('m', 1001);
			
            [x,y] = domains.fractal_profile(xlim, m);
            if isrow(x)
                x = x';
            end
            if isrow(y)
                y = y';
            end
			xl = x(1);
			xr = x(end);
            yb = ylim(1);
            yt = ylim(2);

            % TODO: Determine optimal number of grid points in y-direction
            % ymin = min(y);
            % ymax = max(y);
            
            parametrizedFault = parametrization.dataSpline([x,y]);

            % ---- Block 1 (bot) -------------------------

            % Sides
            NW = parametrizedFault(0);
            NE = parametrizedFault(1);
            SW = [xl; yb];
            SE = [xr; yb];

            W = parametrization.Curve.line(NW, SW);
            E = parametrization.Curve.line(SE, NE);
            S = parametrization.Curve.line(SW, SE);
            N = parametrizedFault.reverse();

            B1 = parametrization.Ti(S, E, N, W);

            % ---- Block 2 (top) -------------------------

            % Corners
            NW = [xl; yt];
            NE = [xr; yt];
            SW = parametrizedFault(0);
            SE = parametrizedFault(1);

            W = parametrization.Curve.line(NW, SW);
            E = parametrization.Curve.line(SE, NE);
            S = parametrizedFault;
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
            H1 = 0 - obj.ylim(1);
            H2 = obj.ylim(2) - 0;

            hFault = faultLength/(mFault - 1);

            % Same spacings in y-directions as on fault
            my1 = round(H1/hFault) + 1;
            my2 = round(H2/hFault) + 1;

            ms = {[mFault, my1], [mFault, my2]};

        end

	end

end
