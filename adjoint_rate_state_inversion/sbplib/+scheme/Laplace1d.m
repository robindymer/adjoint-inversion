classdef Laplace1d < scheme.Scheme
    properties
        grid
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        H % Discrete norm
        M % Derivative norm
        a

        D2
        Hi
        e_l
        e_r
        d_l
        d_r
        % Borrowing constants
        theta_M
        theta_R
        theta_H
    end

    methods
        function obj = Laplace1d(grid, order, a, opset)
            default_arg('a', 1);

            assertType(grid, 'grid.Cartesian');

            ops = opset(grid.size(), grid.lim{1}, order);
            
            if isa(ops.D2,'function_handle')
                ops.D2 = ops.D2(ones(grid.size,1));
            end

            obj.D2 = sparse(ops.D2);
            obj.H =  sparse(ops.H);
            obj.Hi = sparse(ops.HI);
            obj.M =  sparse(ops.M);
            obj.e_l = sparse(ops.e_l);
            obj.e_r = sparse(ops.e_r);
            obj.d_l = -sparse(ops.d1_l);
            obj.d_r = sparse(ops.d1_r);

            obj.grid = grid;
            obj.order = order;

            obj.a = a;
            obj.D = a*obj.D2;

            h = ops.h;
            obj.theta_M = h*ops.borrowing.M.d1;
            obj.theta_R = h*ops.borrowing.R.delta_D;
            obj.theta_H = h*ops.borrowing.H11;
        end


        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary,type,data)
            default_arg('type','neumann');
            default_arg('data',0);
            
            e = obj.getBoundaryOperator('e', boundary);
            d = obj.getBoundaryOperator('d', boundary);
            
            Hi = obj.Hi;
            a = obj.a;
            th_R = obj.theta_R;
            th_H = obj.theta_H;

            switch type
                % Dirichlet boundary condition

                case {'D','d','dirichlet'}
                    tuning = 1.0;

                    tau = tuning*(1/th_R + 1/th_H);

                    closure = a*Hi*d*e' ...
                             -a*Hi*e*tau*e';

                    penalty = -a*Hi*d ...
                              +a*Hi*e*tau;

                % Neumann boundary condition
                case {'N','n','neumann'}
                    tau = -e;

                    closure = obj.a*obj.Hi*tau*d';
                    penalty = -obj.a*obj.Hi*tau;

                % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj, boundary, neighbour_scheme, neighbour_boundary, type)
            default_arg('type',struct);
            default_field(type,'tuning',1);
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            e_u = obj.getBoundaryOperator('e', boundary);
            d_u = obj.getBoundaryOperator('d', boundary);
            b_b_u = 1;

            e_v = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);
            d_v = neighbour_scheme.getBoundaryOperator('d', neighbour_boundary);
            b_b_v = 1;

            % Penalty strenght
            tuning = type.tuning;
            th_H_u = obj.theta_H;
            th_R_u = obj.theta_R;
            th_H_v = neighbour_scheme.theta_H;
            th_R_v = neighbour_scheme.theta_R;
            tau = 1/4*tuning*(b_b_u*(1/th_R_u + 1/th_H_u) + b_b_v*(1/th_R_v + 1/th_H_v));
           
            % Operators/coefficients that are only required from this side
            Hi = obj.Hi;
            a = obj.a;

            closure =   1/2*a*Hi*d_u*b_b_u*e_u' ...
                        -1/2*a*Hi*e_u*b_b_u*d_u' ...
                        -a*Hi*e_u*tau*e_u';

            penalty =   -1/2*a*Hi*d_u*b_b_u*e_v' ...
                        -1/2*a*Hi*e_u*b_b_v*d_v' ...
                        +a*Hi*e_u*tau*e_v';
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string
        % boundary  -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(op, {'e', 'd'})
            assertIsMember(boundary, {'l', 'r'})

            o = obj.([op, '_', boundary]);
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        % Note: for 1d diffOps, the boundary quadrature is the scalar 1.
        function H_b = getBoundaryQuadrature(obj, boundary)
            assertIsMember(boundary, {'l', 'r'})

            H_b = 1;
        end

        % Returns the boundary sign. The right boundary is considered the positive boundary
        % boundary -- string
        function s = getBoundarySign(obj, boundary)
            assertIsMember(boundary, {'l', 'r'})

            switch boundary
                case {'r'}
                    s = 1;
                case {'l'}
                    s = -1;
            end
        end

        function N = size(obj)
            N = obj.grid.size();
        end

    end
end
