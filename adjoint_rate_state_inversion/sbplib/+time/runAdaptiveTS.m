function ts = runAdaptiveTS(ts, T, plotCallback)
	default_arg('plotCallback', []);

	[~, t] = ts.getV();
    next_plot_time = t;
	while t < T
		% Make sure we don't overshoot final time
		if t + ts.k > T
			ts.k = T-t;
		end
		ts.step();
        [~, t] = ts.getV();
		if ~isempty(plotCallback)
			next_plot_time = plotCallback(next_plot_time);
        end
	end
end