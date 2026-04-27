classdef Staggered < grid.Structured
    properties
        gridGroups  % Cell array of grid groups, each group is a cell array
        nGroups % Number of grid groups
        h      % Interior grid spacing
        d      % Number of dimensions
        logic  % Grid in logical domain, if any.
    end

    methods

        % Accepts multiple grids and combines them into a staggered grid
        % Each grid entry is a cell array of grids that store the same field
        function obj = Staggered(d, varargin)
            default_arg('d', 2);

            obj.d = d;

            obj.nGroups = length(varargin);
            obj.gridGroups = cell(obj.nGroups, 1);
            for i = 1:obj.nGroups
                obj.gridGroups{i} = varargin{i};
            end

            obj.h = [];
            obj.logic = [];
        end

        % N returns the number of points in the first grid group
        function o = N(obj)
            o = 0;
            gs = obj.gridGroups{1};
            for i = 1:length(gs)
                o = o+gs{i}.N();
            end
        end

        % D returns the spatial dimension of the grid
        function o = D(obj)
            o = obj.d;
        end

        % size returns a reference size
        function m = size(obj)
            m = obj.gridGroups{1}{1};
        end

        % points returns an n x 1 vector containing the coordinates for the first grid group.
        function X = points(obj)
            X = [];
            gs = obj.gridGroups{1};
            for i = 1:length(gs)
                X = [X; gs{i}.points()];
            end
        end

        % matrices returns a cell array with coordinates in matrix form.
        % For 2d case these will have to be transposed to work with plotting routines.
        function X = matrices(obj)
            error('grid:Staggered1d:matrices', 'Not implemented')
        end

        function h = scaling(obj)
            if isempty(obj.h)
                error('grid:Staggered1d:NoScalingSet', 'No scaling set')
            end

            h = obj.h;
        end

        % Restricts the grid function gf on obj to the subgrid g.
        % Only works for even multiples
        function gf = restrictFunc(obj, gf, g)
            error('grid:Staggered1d:NotImplemented','This method does not exist yet')
        end

        % Projects the grid function gf on obj to the grid g.
        function gf = projectFunc(obj, gf, g)
            error('grid:Staggered1d:NotImplemented','This method does not exist yet')
        end

        % Return the names of all boundaries in this grid.
        function bs = getBoundaryNames(obj)
            switch obj.d()
                case 1
                    bs = {'l', 'r'};
                case 2
                    bs = {'w', 'e', 's', 'n'};
                case 3
                    bs = {'w', 'e', 's', 'n', 'd', 'u'};
                otherwise
                    error('not implemented');
            end
        end

        % Return coordinates for the given boundary
        % gridGroup (scalar)    - grid group to return coordinates for
        % subGrids (array)      - specifies which grids in the grid group to include (default: all grids in the grid group)
        function X = getBoundary(obj, name, gridGroup, subGrids)

            default_arg('gridGroup' , 1);
            grids = obj.gridGroups{gridGroup};
            default_arg('subGrids' , 1:numel(grids));

            X = [];
            for i = subGrids
                X = [X; grids{i}.getBoundary(name)];
            end
        end

    end
end