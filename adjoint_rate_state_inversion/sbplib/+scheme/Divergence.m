classdef Divergence < scheme.Scheme

% Approximates the divergence
% Interface and boundary condition methods are just dummies

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order of accuracy for the approximation

        D
        D1
        H
    end

    methods

        function obj = Divergence(g, order, opSet)
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});

            dim = 2;

            m = g.size();
            m_tot = g.N();

            h = g.scaling();
            lim = g.lim;
            if isempty(lim)
                x = g.x;
                lim = cell(length(x),1);
                for i = 1:length(x)
                    lim{i} = {min(x{i}), max(x{i})};
                end
            end

            % 1D operators
            ops = cell(dim,1);
            for i = 1:dim
                ops{i} = opSet{i}(m(i), lim{i}, order);
            end

            I = cell(dim,1);
            D1 = cell(dim,1);

            for i = 1:dim
                I{i} = speye(m(i));
                D1{i} = ops{i}.D1;
            end

            %====== Assemble full operators ========

            % D1
            obj.D1{1} = kron(D1{1},I{2});
            obj.D1{2} = kron(I{1},D1{2});

            I_dim = speye(dim, dim);

            % E{i}^T picks out component i.
            E = cell(dim,1);
            I = speye(m_tot,m_tot);
            for i = 1:dim
                e = sparse(dim,1);
                e(i) = 1;
                E{i} = kron(I,e);
            end

            Div = sparse(m_tot, dim*m_tot);
            for i = 1:dim
                Div = Div + obj.D1{i}*E{i}';
            end
            obj.D = Div;
            obj.H = [];
            %=========================================%'

        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       bc                  is a cell array of component and bc type, e.g. {1, 'd'} for Dirichlet condition
        %                           on the first component. Can also be e.g.
        %                           {'normal', 'd'} or {'tangential', 't'} for conditions on
        %                           tangential/normal component.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, bc)
            error('Not implemented')
        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- tuning:           penalty strength, defaults to 1.2
        %          -- interpolation:    type of interpolation, default 'none'
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            error('Not implemented')
        end

        function N = size(obj)
            N = obj.dim*prod(obj.m);
        end
    end
end
