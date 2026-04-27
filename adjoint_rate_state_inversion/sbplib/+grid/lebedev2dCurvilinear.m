% Creates a curvilinear 2d lebedev2d grid
% over the logical domain xi_lim, eta_lim, ...
% If all limits are ommited they are set to {0,1}.
% Examples:
%   g = grid.lebedev2dCurvilinear(mapping, [m_xi, m_eta])
%   g = grid.lebedev2dCurvilinear(mapping, [m_xi, m_eta], xi_lim, eta_lim)
%   g = grid.lebedev2dCurvilinear(mapping, [10, 15], {0,1}, {0,1})
function g = lebedev2dCurvilinear(mapping, m, varargin)
    if isempty(varargin)
        varargin = repmat({{0,1}}, [1 length(m)]);
    end

    if length(m) ~= length(varargin)
        error('grid:lebedev2d:NonMatchingParameters','The number of provided dimensions do not match.')
    end

    for i = 1:length(m)
        if ~iscell(varargin{i}) || numel(varargin{i}) ~= 2
           error('grid:lebedev2d:InvalidLimits','The limits should be cell arrays with 2 elements.');
        end

        if varargin{i}{1} > varargin{i}{2}
            error('grid:lebedev2d:InvalidLimits','The elements of the limit must be increasing.');
        end
    end

    g_logic = grid.lebedev2d(m, varargin{:});

    gu1_logic = g_logic.gridGroups{1}{1};
    gu2_logic = g_logic.gridGroups{1}{2};
    gs1_logic = g_logic.gridGroups{2}{1};
    gs2_logic = g_logic.gridGroups{2}{2};

    gu1 = grid.Curvilinear(mapping, gu1_logic.x{1}, gu1_logic.x{2});
    gu2 = grid.Curvilinear(mapping, gu2_logic.x{1}, gu2_logic.x{2});
    gs1 = grid.Curvilinear(mapping, gs1_logic.x{1}, gs1_logic.x{2});
    gs2 = grid.Curvilinear(mapping, gs2_logic.x{1}, gs2_logic.x{2});

    gu = {gu1, gu2};
    gs = {gs1, gs2};

    dim = 2;
    g = grid.Staggered(dim, gu, gs);

    g.logic = g_logic;
end