function S = addFunctionHandles(S_cell)
	S = [];
	for i = 1:numel(S_cell)
        data = S_cell{i};
        if ~isempty(data)
            if isempty(S)
                S = data;
            else
                S = @(t) S(t) + data(t);
            end
        end
    end

end
