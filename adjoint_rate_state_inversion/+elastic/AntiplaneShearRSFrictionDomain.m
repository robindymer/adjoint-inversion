classdef AntiplaneShearRSFrictionDomain < multiblock.Definition
    properties
    
    nBlocks
    xlims
    blockNames
    connections
    boundaryGroups % Structure of boundaryGroups

    end


    methods
        function obj = AntiplaneShearRSFrictionDomain(x,blockNames)
            default_arg('blockNames',[]);
            N = 2; % Number of blocks
            assert(length(x) == N+1, 'x should have length 3, i.e. x = [x_l, x_i, x_r]');
            assert(issorted(x), 'The elements of x seem to be in the wrong order')
            
            x_l = x(1);
            x_i = x(2);
            x_r = x(3);
            
        
            % Block limits
            xlims = cell(N,1);
            xlims{1} = {x_l, x_i};
            xlims{2} = {x_i, x_r};
            
            % The blocks are connected through a frictional interface, and thus cant
            % be treated using the normal interface coupling routines. 
            % Therefore connections `conns` are left empty.
            conn = cell(N,N);
            
            % Block names (id number as default)
            if isempty(blockNames)
                obj.blockNames = cell(1, N);
                for i = 1:N
                    obj.blockNames{i} = sprintf('%d', i);
                end
            else
                assert(length(blockNames) == N);
                obj.blockNames = blockNames;
            end

            % Boundary groups
            boundaryGroups = struct();
            left = {1,'l'};
            right = {2,'r'};
            interface_m = {1,'r'};
            interface_p = {2,'l'};
            boundaryGroups.left = multiblock.BoundaryGroup({left});
            boundaryGroups.right = multiblock.BoundaryGroup({right});
            boundaryGroups.outer = multiblock.BoundaryGroup({left, right});
            boundaryGroups.interface = multiblock.BoundaryGroup({interface_m, interface_p});
            
            obj.nBlocks = N;
            obj.xlims = xlims;
            obj.connections = conn;
            obj.boundaryGroups = boundaryGroups;
        end


        % Returns a multiblock.Grid given some parameters
        % ms: vector of grid points [m_m, m_p] for the two blocks
        % For same m in both blocks, just input scalar m.
        % Currently defaults to equidistant grid if varargin is empty.
        % If varargin is non-empty, the first argument should supply the grid type, followed by
        % additional arguments required to construct the grid.
        % Grid types:
        %          'equidist' - equidistant grid
        %                       Additional argumets: none
        %          'boundaryopt' - boundary optimized grid based on boundary
        %                          optimized SBP operators
        %                          Additional arguments: order, stencil option
        % Example: g = getGrid() - the local blocks are 101 equidistant grids.
        %          g = getGrid(ms,) - block i is an equidistant grid with size given by ms{i}.
        %          g = getGrid(ms,'equidist') - block i is an equidistant grid with size given by ms{i}.
        %          g = getGrid(ms,'boundaryopt',4,'acc') - block i is a Cartesian grid with size given by ms{i}
        %              and nodes placed according to the boundary optimized accurate 4th order SBP operator.
        function g = getGrid(obj, ms, varargin)
            default_arg('ms',[101])

            % Extend ms if input is a scalar
            if (numel(ms) == 1)
                ms = [ms, ms];
            end
            if isempty(varargin) || strcmp(varargin{1},'equidist')
               gridgenerator = @(m,xlim) grid.equidistant(m, xlim);
            elseif strcmp(varargin{1},'boundaryopt')
                order = varargin{2};
                stenciloption = varargin{3};
                gridgenerator = @(m,xlim) grid.boundaryOptimized(m, xlim, order, stenciloption);
            else
                error('No grid type supplied!');
            end
            grids = cell(1, obj.nBlocks);
            for i = 1:obj.nBlocks
                grids{i} = gridgenerator(ms(i), obj.xlims{i});
            end
            g = multiblock.Grid(grids, obj.connections, obj.boundaryGroups);
        end

        % label is the type of label used for plotting,
        % default is block name, 'id' show the index for each block.
        function show(obj, label, gridLines, varargin)
            default_arg('label', 'name')
            default_arg('gridLines', false);
        
            x_l = obj.xlims{1}{1};
            x_i = obj.xlims{1}{2};
            x_r = obj.xlims{2}{2};
            
            if isempty(label) && ~gridLines
                x = [x_l x_i x_r];
                plot(x,0*x,'.','MarkerSize',8);
                hold on;
                xline(x_i,'--','linewidth',1.5);
                hold off;
                axis equal
                return
            end

            if gridLines
                if isempty(varargin)
                    m = 10;
                else
                    m = varargin{1};
                end
                g = getGrid(obj, m, varargin{2:end});
                x = g.points();
                plot(x,0*x,'.','MarkerSize',8);
                hold on;
                xline(x_i,'--','linewidth',1.5);
                hold off;
            end
            
            if ~isempty(label)
                switch label
                case 'name'
                    label1 = sprintf('$\\Omega^%s$',obj.blockNames{1});
                    label2 = sprintf('$\\Omega^%s$',obj.blockNames{2});
                    labels = {label1, label2};          
                case 'id'
                    labels = {'$\Omega^1$', '$\Omega^2$'};
                end
                xm_m = (x_l + x_i)/2;
                xm_p = (x_i + x_r)/2;
                text(xm_m,0.2,labels{1},'interpreter','latex');
                text(xm_p,0.2,labels{2},'interpreter','latex');
            end
            axis equal
        end
    end
end
