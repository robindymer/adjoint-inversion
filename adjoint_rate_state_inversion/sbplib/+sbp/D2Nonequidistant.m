classdef D2Nonequidistant < sbp.OpSet
    % Implements the boundary optimized variable coefficient
    % second derivative.
    %
    % The boundary closure uses the first and last rows of the
    % boundary-optimized D1 operator, i.e. the operators are
    % fully compatible.
    properties
        H % Norm matrix
        HI % H^-1
        D1 % SBP operator approximating first derivative
        D2 % SBP operator approximating second derivative
        DI % Dissipation operator
        e_l % Left boundary operator
        e_r % Right boundary operator
        d1_l % Left boundary first derivative
        d1_r % Right boundary first derivative
        m % Number of grid points.
        h % Step size
        x % grid
        M % Norm matrix, second derivative
        borrowing % Struct with borrowing limits for different norm matrices
        options % Struct holding options used to create the operator
    end

    methods

        function obj = D2Nonequidistant(m, lim, order, options)
            % m - number of gridpoints
            % lim - cell array holding the limits of the domain
            % order - order of the operator
            % options - struct holding options used to construct the operator
            %           struct.stencil_width:   {'minimal', 'nonminimal', 'wide'}
            %               minimal: minimal compatible stencil width (default)
            %               nonminimal: a few additional stencil points compared to minimal
            %               wide: wide stencil obtained by applying D1 twice
            %           struct.AD: {'op', 'upwind'}
            %               'op': order-preserving AD (preserving interior stencil order) (default)
            %               'upwind': upwind AD (order-1 upwind interior stencil)
            %           struct.variable_coeffs: {true, false}
            %               true: obj.D2 is a function handle D2(c) returning a matrix 
            %                     for coefficient vector c (default)
            %               false: obj.D2 is a matrix.
            default_arg('options', struct);
            default_field(options,'stencil_width','minimal');
            default_field(options,'AD','op');
            default_field(options,'variable_coeffs',true);
            [x, h] = sbp.grid.accurateBoundaryOptimizedGrid(lim, m, order);
            switch order
                case 4
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = sbp.implementations.d2_noneq_variable_4(m, h, options);
                case 6
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = sbp.implementations.d2_noneq_variable_6(m, h, options);
                case 8
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = sbp.implementations.d2_noneq_variable_8(m, h, options);
                case 10
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = sbp.implementations.d2_noneq_variable_10(m, h, options);
                case 12
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = sbp.implementations.d2_noneq_variable_12(m, h, options);
                otherwise
                    error('Invalid operator order %d.', order);
            end

            if ~options.variable_coeffs
                obj.D2 = obj.D2(ones(m,1));
            end

            % Boundary operators
            obj.e_l = sparse(m, 1); obj.e_l(1) = 1;
            obj.e_r = sparse(m, 1); obj.e_r(m) = 1;
            obj.d1_l = (obj.e_l' * obj.D1)';
            obj.d1_r = (obj.e_r' * obj.D1)';

            obj.M = [];

            % Borrowing coefficients
            obj.borrowing.H11 = obj.H(1, 1) / h; % First element in H/h,
            obj.borrowing.M.d1 = obj.H(1, 1) / h; % First element in H/h is borrowing also for M
            obj.borrowing.R.delta_D = inf;

            % grid data
            obj.x = x;
            obj.h = h;
            obj.m = m;

            % misc
            obj.options = options;
        end

        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end

    end

end
