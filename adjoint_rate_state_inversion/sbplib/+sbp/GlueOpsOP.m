% Order-preserving (OP) glue operators.
%
% Let ^* denote the adjoint. These operators satsify
%
% Iuv2.good = Iv2u.bad^*
% Iv2u.good = Iu2v.bad^*
%
% The .bad operators have the same order of accuracy as the operators
% by Mattsson and Carpenter (MC) in InterpOpsMC, i.e. order p,
% if the interior stencil is order 2p. The .good operators are
% one order more accurate, i.e. order p+1.
%
% For PDEs of second order in space, the OP operators allow for the same
% convergence rate as with conforming interfaces, which is an improvement
% by one order compared what is possible with the MC operators.
classdef GlueOpsOP < sbp.InterpOps
    properties

        % Structs of interpolation operators, fields .good and .bad
        Iu2v
        Iv2u
    end

    methods
        % m_u, m_v         --   number of grid points along the interface
        % order_u, order_v --   order of accuracy in the different blocks
        function obj = GlueOpsOP(m_u, m_v, order_u, order_v, opset_u, opset_v, order_preserving)
            default_arg('opset_u',@sbp.D2Variable);
            default_arg('opset_v',@sbp.D2Variable);
            default_arg('order_preserving', true);
            Iu2v = struct;
            Iv2u = struct;
            function [g, optype] = interp_grid_and_optype(m, order, opset)
                switch toString(opset)
                case {'sbp.D2Standard','sbp.D2Variable','D2VariableCompatible','sbp.D4Standard','sbp.D4Variable'}
                    optype = 'trad';
                    g = grid.equidistant(m, {0,1});
                    
                case {'D1Nonequidistant','sbp.D2Nonequidistant'}
                    optype = 'bopt';
                    g = grid.boundaryOptimized(m, {0,1}, order);
                end
            end

            [g_u, optype_u] = interp_grid_and_optype(m_u, order_u, opset_u);
            [g_v, optype_v] = interp_grid_and_optype(m_v, order_v, opset_v);

            % TODO: The h-dependency should be removed from the projection operator implementation. The resulting
            % interpolation operators should be independent of h.
            if order_preserving
                [Iu2v.good, Iv2u.bad] = sbp.implementations.projOpsGridToGrid(optype_u, g_u.points, g_u.h, order_u, 'fd2glue',...
                                                                            optype_v, g_v.points, g_v.h, order_v, 'glue2fd');
                [Iu2v.bad, Iv2u.good] = sbp.implementations.projOpsGridToGrid(optype_u, g_u.points, g_u.h, order_u, 'glue2fd',...
                                                                            optype_v, g_v.points, g_v.h, order_v, 'fd2glue');
            else
                [Iu2v.good, Iv2u.bad] = sbp.implementations.projOpsGridToGrid(optype_u, g_u.points, g_u.h, order_u, 'none',...
                                                                              optype_v, g_v.points, g_v.h, order_v, 'none');
                Iu2v.bad = Iu2v.good;
                Iv2u.good = Iv2u.bad;
            end                                                                              
            obj.Iu2v = Iu2v;
            obj.Iv2u = Iv2u;
        end

        function str = string(obj)
            str = [class(obj)];
        end

    end
end
