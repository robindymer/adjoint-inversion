classdef Laplace1dVariable < scheme.Scheme
    properties
        grid
        order % Order accuracy for the approximation

        D % scheme operator
        H % Discrete norm

        % Variable coefficients a(b u_x)_x
        a
        b

        Hi
        e_l
        e_r
        d_l
        d_r
        gamm  % borrowing constant (includes h)
    end

    methods
        function obj = Laplace1dVariable(g, order, a, b, opset)
            default_arg('a', @(x) 0*x + 1);
            default_arg('b', @(x) 0*x + 1);
            default_arg('opset', @sbp.D2Variable);
            assertType(g, 'grid.Cartesian');

            ops = opset(g.size(), g.lim{1}, order);
            obj.H =  sparse(ops.H);
            obj.Hi = sparse(ops.HI);
            obj.e_l = sparse(ops.e_l);
            obj.e_r = sparse(ops.e_r);
            obj.d_l = sparse(ops.d1_l);
            obj.d_r = sparse(ops.d1_r);

            obj.grid = g;
            obj.order = order;

            a = grid.evalOn(g, a);
            b = grid.evalOn(g, b);

            A = spdiag(a);
            B = spdiag(b);

            obj.a = A;
            obj.b = B;

            obj.D = A*ops.D2(b);
            
            % Borrowing constant
            switch toString(opset)
            case {'sbp.D2Variable','sbp.D4Variable'} % Compatible operators
                obj.gamm = g.h*ops.borrowing.M.d1;
            case {'sbp.D2VariableCompatible','sbp.D2Nonequidistant'} % Fully compatible operators
                obj.gamm = g.h*ops.borrowing.H11;
            end
        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        function [closure, penalty] = boundary_condition(obj,boundary,type,tuning)
            default_arg('type','neumann');
            default_arg('tuning',1); % Tuning used for Dirichlet conditions

            e = obj.getBoundaryOperator('e', boundary);
            d = obj.getBoundaryOperator('d', boundary);
            s = obj.getBoundarySign(boundary);

            d_n = s*d; % Normal derivative;
            b_b = e'*obj.b*e; % Coefficient on boundary
            Hi = obj.Hi;
            a = obj.a;

            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    tau1 = -tuning/obj.gamm;
                    tau2 =  1;

                    tau = tau1*e + tau2*d_n;

                    closure = a*Hi*tau*b_b*e';
                    penalty = -a*Hi*tau*b_b;

                % Neumann boundary condition
                case {'N','n','neumann'}
                    tau = -e;
                    closure = a*Hi*tau*b_b*d_n';
                    penalty = -a*Hi*tau*b_b;
                
                % Traction boundary condition
                case {'T','t','traction'}
                    tau = -e;
                    closure = a*Hi*tau*b_b*d_n';
                    penalty = -a*Hi*tau;

                % From Erickson 2022: tau = -\alpha u_t
                case {'erickson2022'}
                    % Traction (normal derivative) flux
                    closure = -e*b_b*d_n';

                    % Displacement (u) flux
                    closure = closure + d_n*b_b*e';

                    % Multiply by a*Hi
                    closure = a*Hi*closure;

                    % Two penalties needed. Set to zero here and handled in discr.
                    penalty = 0*e;

                % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj, boundary, neighbour_scheme, neighbour_boundary, type, tuning)
            default_arg('type',[]);
            default_arg('tuning',1.1); % Penalty parameter tuning
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            e_u = obj.getBoundaryOperator('e', boundary);
            d_u = obj.getBoundaryOperator('d', boundary);
            s_u = obj.getBoundarySign(boundary);
            d_n_u = s_u*d_u; % Normal derivative

            e_v = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);
            d_v = neighbour_scheme.getBoundaryOperator('d', neighbour_boundary);
            s_v = neighbour_scheme.getBoundarySign(neighbour_boundary);
            d_n_v = s_v*d_v; % Normal derivative

            b_u = e_u'*obj.b*e_u;
            b_v = e_v'*neighbour_scheme.b*e_v;

            gamm_u = obj.gamm;
            gamm_v = neighbour_scheme.gamm;

            closure = zeros(size(obj.D));
            penalty = zeros(length(obj.D), length(b_v));

            % Continuity of bu_x
            closure = closure - 1/2*e_u*b_u*d_n_u';
            penalty = penalty - 1/2*e_u*b_v*d_n_v';

            % Continuity of u (symmetrizing term)
            closure = closure + 1/2*d_n_u*b_u*e_u';
            penalty = penalty - 1/2*d_n_u*b_u*e_v';

            % Continuity of u (symmetric term)
            tau = 1/4*(b_u/gamm_u + b_v/gamm_v)*tuning;
            closure = closure - e_u*tau*e_u';
            penalty = penalty + e_u*tau*e_v';

            % Multiply by Hi and a.
            closure = obj.Hi*obj.a*closure;
            penalty = obj.Hi*obj.a*penalty;

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
