classdef DefCurvilinear < multiblock.Definition
    properties
        nBlocks
        blockMaps % Maps from logical blocks to physical blocks build from transfinite interpolation
        blockNames
        connections % Cell array specifying connections between blocks
        boundaryGroups % Structure of boundaryGroups
    end

    methods
        % Defines a multiblock setup for transfinite interpolation blocks
        % TODO: How to bring in plotting of points?
        function obj = DefCurvilinear(blockMaps, connections, boundaryGroups, blockNames)
            default_arg('boundaryGroups', struct());
            default_arg('blockNames',{});

            nBlocks = length(blockMaps);

            obj.nBlocks = nBlocks;

            obj.blockMaps = blockMaps;

            assert(all(size(connections) == [nBlocks, nBlocks]));
            obj.connections = connections;


            if isempty(blockNames)
                obj.blockNames = cell(1, nBlocks);
                for i = 1:length(blockMaps)
                    obj.blockNames{i} = sprintf('%d', i);
                end
            else
                assert(length(blockNames) == nBlocks);
                obj.blockNames = blockNames;
            end

            obj.boundaryGroups = boundaryGroups;
        end

        % Returns a multiblock.Grid given some parameters
        % ms: cell array of [mx, my] vectors
        % Currently defaults to an equidistant curvilinear grid if varargin is empty.
        % If varargin is non-empty, the first argument should supply the grid type, followed by
        % additional arguments required to construct the grid.
        % Grid types:
        %          'equidist' - equidistant curvilinear grid
        %                       Additional argumets: none
        %          'boundaryopt' - boundary optimized grid based on boundary
        %                          optimized SBP operators
        %                          Additional arguments: order, stencil option
        function g = getGrid(obj, ms, varargin)
            % If a scalar is passed, defer to getGridSizes implemented by subclass
            % TODO: This forces the interface of subclasses.
            % Should ms be included in varargin? Figure out bow to do it properly
            if ~iscell(ms) && length(ms) == 1
                ms = obj.getGridSizes(ms);
            end
            if isempty(varargin) || strcmp(varargin{1},'equidist')
                gridgenerator = @(blockMap,m) grid.equidistantCurvilinear(blockMap, m);
            elseif strcmp(varargin{1},'boundaryopt')
                 order = varargin{2};
                 stenciloption = varargin{3};
                 gridgenerator = @(blockMap,m) grid.boundaryOptimizedCurvilinear(blockMap,m,{0,1},{0,1},...
                     order,stenciloption);
            elseif strcmp(varargin{1},'general')
                gridgenerator = @(blockMap,m) grid.generalCurvilinear(blockMap,m,{0,1},{0,1},varargin{2:end});
            else
                error('No grid type supplied!');
            end
            grids = cell(1, obj.nBlocks);
            for i = 1:obj.nBlocks
                grids{i} = gridgenerator(obj.blockMaps{i}.S, ms{i});
            end

            g = multiblock.Grid(grids, obj.connections, obj.boundaryGroups);
        end

        function g = getLebedevGrid(obj, varargin)
            ms = obj.getGridSizes(varargin{:});

            grids = cell(1, obj.nBlocks);
            for i = 1:obj.nBlocks
                % grids{i} = grid.equidistantCurvilinear(obj.blockMaps{i}.S, ms{i});
                grids{i} = grid.lebedev2dCurvilinear(obj.blockMaps{i}.S, ms{i});
            end

            g = multiblock.Grid(grids, obj.connections, obj.boundaryGroups);
        end

        function h = show(obj, label, gridLines, varargin)
            default_arg('label', 'name')
            default_arg('gridLines', false);

            h = [];
            if isempty('label') && ~gridLines
                for i = 1:obj.nBlocks
                    h = [h, obj.blockMaps{i}.show(2,2)];
                end
                axis equal
                return
            end

            if gridLines
                ms = obj.getGridSizes(varargin{:});
                for i = 1:obj.nBlocks
                    h = [h, obj.blockMaps{i}.show(ms{i}(1),ms{i}(2))];
                end
            end


            switch label
                case 'name'
                    labels = obj.blockNames;
                case 'id'
                    labels = {};
                    for i = 1:obj.nBlocks
                        labels{i} = num2str(i);
                    end
                case 'none'
                    axis equal
                    return
            end

            for i = 1:obj.nBlocks
                parametrization.Ti.label(obj.blockMaps{i}, labels{i});
            end

            axis equal
        end
    end

    methods (Abstract)
        % Returns the grid size of each block in a cell array
        % The input parameters are determined by the subclass
        ms = getGridSizes(obj, varargin)
        % end
    end

end


