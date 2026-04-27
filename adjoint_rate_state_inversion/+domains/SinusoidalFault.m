classdef SinusoidalFault < multiblock.DefCurvilinear 

	properties
		xlim
		ylim
	end

	methods
		function obj = SinusoidalFault(xlim, ylim)
			default_arg('xlim', [-1, 0, 1]);
			default_arg('ylim', [-1, 0]);

			xm = xlim(1);
			x0 = xlim(2);
			xp = xlim(3);
			yb = ylim(1);
			yt = ylim(2);

            % ---- Fault shape ---------
            L = x0 - xm;
            H = yt - yb;
            A = L/10;
            wvl = H/2;
            k = elastic.wavelengthToWavenumber(wvl, 1);

            faultShapeFun = @(y) A*sin(k*y);
            g = @(y) [x0 + faultShapeFun(y); y]; 

            % Reparametrize curve
            gFault = @(t) g(yb + t*(yt-yb));
            CFault = parametrization.Curve(gFault);
			% ---- Block 1 (left) -------------------------

            % Corners
            NW = [xm; yt];
            NE = [x0; yt];
            SW = [xm; yb];
            SE = [x0; yb];

            W = parametrization.Curve.line(NW, SW);
            S = parametrization.Curve.line(SW, SE);
            N = parametrization.Curve.line(NE, NW);
            E = CFault;

            B1 = parametrization.Ti(S, E, N, W);

            % ---- Block 2 (right) -------------------------

            % Corners
            NW = [x0; yt];
            NE = [xp; yt];
            SW = [x0; yb];
            SE = [xp; yb];

            W = CFault.reverse();
            E = parametrization.Curve.line(SE, NE);
            S = parametrization.Curve.line(SW, SE);
            N = parametrization.Curve.line(NE, NW);

            B2 = parametrization.Ti(S, E, N, W);
            %-----------------------------------------------

            blocks = {B1, B2};
            blockNames = {'left', 'right'};

            conn = cell(2,2);

            boundaryGroups = struct;

            boundaryGroups.surface = multiblock.BoundaryGroup({{1,'n'},{2,'n'}});
            boundaryGroups.bottom = multiblock.BoundaryGroup({{1,'s'},{2,'s'}});

            boundaryGroups.left = multiblock.BoundaryGroup({{1,'w'}});
            boundaryGroups.right = multiblock.BoundaryGroup({{2,'e'}});
            boundaryGroups.sides = multiblock.joinBoundaryGroups(boundaryGroups.left, boundaryGroups.right);

            boundaryGroups.faultLeft = multiblock.BoundaryGroup({{1,'e'}});
            boundaryGroups.faultRight = multiblock.BoundaryGroup({{2,'w'}});
            boundaryGroups.fault = multiblock.joinBoundaryGroups(boundaryGroups.faultLeft, boundaryGroups.faultRight);

            boundaryGroups.fault_minus = boundaryGroups.faultLeft;
            boundaryGroups.fault_plus = boundaryGroups.faultRight;


			obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blockNames);

			obj.xlim = xlim;
			obj.ylim = ylim;
		end

		function ms = getGridSizes(obj, mFault)
            % mFault = number of points along fault

            faultLength = obj.ylim(2) - obj.ylim(1);

            % Widths in blocks 1 and 2
            W1 = obj.xlim(2) - obj.xlim(1);
            W2 = obj.xlim(3) - obj.xlim(2);

            hFault = faultLength/(mFault - 1);

            % Same spacings in x-directions as on fault
            mx1 = round(W1/hFault) + 1;
            mx2 = round(W2/hFault) + 1;

            ms = {[mx1, mFault], [mx2, mFault]};

        end

	end

end
