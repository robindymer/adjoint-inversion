classdef D1UpwindOptimal < sbp.OpSet
    properties
        Dp, Dm % SBP operator approximating first derivative
        H % Norm matrix
        HI % H^-1
        e_l % Left boundary operator
        e_r % Right boundary operator
        m % Number of grid points.
        h % Step size
        x % grid
        borrowing % Struct with borrowing limits for different norm matrices
    end

    methods
        function obj = D1UpwindOptimal(m, lim, order)

            [x, h] = sbp.grid.upwindAccurateBoundaryOptimizedGrid(lim, m ,order);
            switch order
                case 2
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_l, obj.e_r] = sbp.implementations.d1_upwind_optimal_2(m,h);
                case 3
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_l, obj.e_r] = sbp.implementations.d1_upwind_optimal_3(m,h);
                case 4
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_l, obj.e_r] = sbp.implementations.d1_upwind_optimal_4(m,h);
                case 5
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_l, obj.e_r] = sbp.implementations.d1_upwind_optimal_5(m,h);
                case 6
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_l, obj.e_r] = sbp.implementations.d1_upwind_optimal_6(m,h);
                case 7
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_l, obj.e_r] = sbp.implementations.d1_upwind_optimal_7(m,h);
                case 8
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_l, obj.e_r] = sbp.implementations.d1_upwind_optimal_8(m,h);
                case 9
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_l, obj.e_r] = sbp.implementations.d1_upwind_optimal_9(m,h);
                otherwise
                    error('No such order implemented')
            end
            obj.m = m;
            obj.x = x;
            obj.h = h;
            obj.borrowing = [];
        end

        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end
    end

end