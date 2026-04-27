classdef StaggeredSurface < handle
    properties
        grid
        surfs
        gridGroup
        subGrid

        ZData
        CData

    end

    methods
        function obj = StaggeredSurface(g, gf, gridGroup, subGrid)

            default_arg('gridGroup', 1);
            default_arg('subGrid', 1);

            obj.grid = g;
            obj.gridGroup = gridGroup;
            obj.subGrid = subGrid;

            % coords = obj.grid.points();
            % X = obj.grid.funcToPlotMatrices(coords(:,1));
            % Y = obj.grid.funcToPlotMatrices(coords(:,2));
            % V = obj.grid.funcToPlotMatrices(gf);
            X = {};
            Y = {};
            V = {};

            holdState = ishold();
            hold on

            surfs = cell(1, obj.grid.nBlocks);
            gfIndex = 1;
            for i = 1:g.nBlocks()

                gi = g.grids{i}.gridGroups{gridGroup}{subGrid};

                X{i} = grid.funcToPlotMatrix(gi, gi.coords(:,1));
                Y{i} = grid.funcToPlotMatrix(gi, gi.coords(:,2));

                Ni = gi.N();
                gf_i = gf(gfIndex:gfIndex+Ni-1);
                V{i} = grid.funcToPlotMatrix(gi, gf_i);

                surfs{i} = surf(X{i}, Y{i}, V{i});
                gfIndex = gfIndex + Ni;
            end

            if holdState == false
                hold off
            end

            obj.surfs = [surfs{:}];

            obj.ZData = gf;
            obj.CData = gf;
        end

        function set(obj, propertyName, propertyValue)
            set(obj.surfs, propertyName, propertyValue);
        end

        function obj = set.ZData(obj, gf)
            obj.ZData = gf;

            % V = obj.grid.funcToPlotMatrices(gf);
            gfIndex = 1;
            for i = 1:obj.grid.nBlocks()
                gi = obj.grid.grids{i}.gridGroups{obj.gridGroup}{obj.subGrid};
                Ni = gi.N();
                gf_i = gf(gfIndex:gfIndex+Ni-1);
                Vi = grid.funcToPlotMatrix(gi, gf_i);
                obj.surfs(i).ZData = Vi;

                gfIndex = gfIndex + Ni;
            end
        end

        function obj = set.CData(obj, gf)
            obj.CData = gf;

            % V = obj.grid.funcToPlotMatrices(gf);
            gfIndex = 1;
            for i = 1:obj.grid.nBlocks()
                gi = obj.grid.grids{i}.gridGroups{obj.gridGroup}{obj.subGrid};
                Ni = gi.N();
                gf_i = gf(gfIndex:gfIndex+Ni-1);
                Vi = grid.funcToPlotMatrix(gi, gf_i);
                obj.surfs(i).CData = Vi;

                gfIndex = gfIndex + Ni;
            end
        end
    end
end
