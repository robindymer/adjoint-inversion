function marmousiMovies(downSampling, saveMovie)
	default_arg('downSampling', 40);
	default_arg('saveMovie', false);

	fontsize = 16;

	filenameGoodAdj = sprintf('goodAdjDS%d', downSampling);
	filenameGoodFwd = sprintf('goodFwdDS%d', downSampling);

	filenameGood = {filenameGoodAdj, filenameGoodFwd};
	% climsPressure = {5e-12, 2};
	climsPressure = {5e-13, 1e-1};

	filenameBad = sprintf('badAdjDS%d', downSampling);

	%---- Elastic acoustic ------
	for k = 1:numel(filenameGood)
		filename = filenameGood{k};
		if saveMovie
			movieName = [filename 'PressurePaddedRickerA2.avi'];
			writerObj = VideoWriter(movieName);
			open(writerObj);
		end

		s = loadFile([filename, '.mat']);
		U = s.u;
		Ut = s.ut;
		T = s.t;
		E_elastic = s.E_elastic;
		E_acoustic = s.E_acoustic;
		g_elastic = s.g.elastic;
		g_acoustic = s.g.acoustic;
		U_elastic = E_elastic'*U;
		Ut_acoustic = E_acoustic'*Ut;
		LAMBDA = s.LAMBDA;
		RHO_acoustic = s.RHO_acoustic;
		Div = s.Div;

		%---- Pressure ------
		P_el = -LAMBDA*Div*U_elastic;
		P_ac = -RHO_acoustic*Ut_acoustic;
		% ------------------

		%--- Setup plot ---
		fhandle = setupPlot();
		Sur_el = multiblock.Surface(g_elastic, P_el(:,1));
		shading interp
		hold on;
		Sur_ac = multiblock.Surface(g_acoustic, P_ac(:,1));
		h = gca;
	    adjustPlot(fhandle, fontsize)
	    % colorbar
	    % -----------------

	    [~, N] = size(U);
		for i = 1:N;
			t = T(i);
			% if t > Tend
			% 	break;
			% end
			u = U(:,i);
			updatePlotElAc(h, Sur_el, Sur_ac, P_el(:,i), P_ac(:,i), t, climsPressure{k});
			pause(0.1)
			drawnow;

			if saveMovie
				%===== Add frame to movie ===%
		        frame = getframe(gcf);
		        writeVideo(writerObj,frame);
		        %===========================%
		    end
		end

		if saveMovie
			%== Close movie object ==%
			close(writerObj);
			%========================%
		end
	end
	%------------------------

	%---- Elastic ----------
	filename = filenameBad;
	if saveMovie
		movieName = [filename 'PressurePaddedRickerA2.avi'];
		writerObj = VideoWriter(movieName);
		open(writerObj);
	end

	s = loadFile([filename '.mat']);
	U = s.u;
	T = s.t;
	g = s.g;
	Div = s.Div;
	LAMBDA = s.LAMBDA;

	%---- Pressure ------
	P = -LAMBDA*Div*U;
	% ------------------

	%--- Setup plot ---
	fhandle = setupPlot();
	Sur = multiblock.Surface(g, P(:,1));
	shading interp
	h = gca;
    adjustPlot(fhandle, fontsize)
    % colorbar
    % -----------------

    [~, N] = size(U);
	for i = 1:N;
		t = T(i);
		% if t > Tend
		% 	break;
		% end
		u = U(:,i);
		updatePlotEl(h, Sur, P(:,i), t, climsPressure{1});
		pause(0.1)
		drawnow;

		if saveMovie
			%===== Add frame to movie ===%
	        frame = getframe(gcf);
	        writeVideo(writerObj,frame);
	        %===========================%
	    end
	end

	if saveMovie
		%== Close movie object ==%
		close(writerObj);
		%========================%
	end
	%------------------------


end

function updatePlotElAc(h, Sur_el, Sur_ac, p_el, p_ac, t, cmax)
	Sur_el.ZData = p_el;
    Sur_el.CData = p_el;
    Sur_ac.ZData = p_ac;
    Sur_ac.CData = p_ac;
    caxis(h, [-cmax, cmax]);
    title(['t = ' sprintf('%4.2f',t)]);
    xlim([0, 17000])
    ylim([-3500, 0])
end

function updatePlotEl(h, Sur, p, t, cmax)
	Sur.ZData = p;
    Sur.CData = p;
    caxis(h, [-cmax, cmax]);
    title(['t = ' sprintf('%4.2f',t)]);
    xlim([0, 17000])
    ylim([-3500, 0])
end

function f = setupPlot()
	scrsz = get(0,'ScreenSize');
    f = figure('Position',[0.05*scrsz(3) 0.05*scrsz(4) 0.95*scrsz(3) 0.65*scrsz(4)]);

    tint = 0.0;
    bgcolor = [1, 1, 1];
    set(gcf, 'Color', bgcolor );
    set(gca, 'Color', bgcolor, 'Xcolor',[tint,tint,tint], ...
        'Ycolor',[tint,tint,tint]);
end

function adjustPlot(fhandle, fontsize)
	xlabel('x')
    ylabel('y')
    shading interp
    setFontSize(fhandle, fontsize);
    axis equal
    xlim([0, 17000])
    ylim([-3500, 0])
end

function s = loadFile(filename)
	load(filename);
	s = saveData;
end




