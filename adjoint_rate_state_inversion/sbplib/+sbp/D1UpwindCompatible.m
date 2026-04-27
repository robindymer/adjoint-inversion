classdef D1UpwindCompatible < sbp.OpSet
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
        function obj = D1UpwindCompatible(m,lim,order)

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
            obj.x = linspace(x_l,x_r,m)';

            ops = sbp.D2Standard(m, lim, order);
            D1 = ops.D1;
            H = ops.H;

            obj.H = H;
            obj.HI = inv(H);
            obj.e_l = ops.e_l;
            obj.e_r = ops.e_r;

            switch order
                case 2
                    ops = sbp.D2Standard(m, lim, 2);
                    Dp = D1 + obj.h^1*1/2*(H\ops.M);
                    Dm = D1 - obj.h^1*1/2*(H\ops.M);
                case 4
                    ops = sbp.D4Variable(m, lim, 2);
                    Dp = D1 - obj.h^3*1/12*(H\ops.M4);
                    Dm = D1 + obj.h^3*1/12*(H\ops.M4);
                otherwise
                    error('Invalid operator order %d.',order);
            end

            obj.Dp = Dp;
            obj.Dm = Dm;

            obj.m = m;
        	obj.borrowing = [];

        end

        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end
    end


end





