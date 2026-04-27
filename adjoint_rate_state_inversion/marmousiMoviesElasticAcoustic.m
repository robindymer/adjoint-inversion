function marmousiMoviesElasticAcoustic(filename, downSampling, saveMovie,...
									climsP, climsVx, climsVy)
	default_arg('filename', 'elAc')
	default_arg('downSampling', 40);
	default_arg('saveMovie', false);

	default_arg('climsP', 6e-4);
	default_arg('climsVx', 1e-9);
	default_arg('climsVy', 1e-9);

	Pstr = 'Pressure';
	Vxstr = 'v_x';
	Vystr = 'v_y';

	fontsize = 16;

	filename = [filename, 'DS', num2str(downSampling)];

	%---- Elastic acoustic ------
	if saveMovie
		movieName = [filename 'Pressure.avi'];
		wObjPressure = VideoWriter(movieName);
		open(wObjPressure);

		movieName = [filename 'Vx.avi'];
		wObjVx = VideoWriter(movieName);
		open(wObjVx);

		movieName = [filename 'Vy.avi'];
		wObjVy = VideoWriter(movieName);
		open(wObjVy);
	end

	s = loadFile([filename, '.mat']);
	T = s.t;
	g_elastic = s.g.elastic;
	g_acoustic = s.g.acoustic;

	gradOp = multiblock.DiffOp(@scheme.Gradient, g_acoustic, 4);
	grad = gradOp.D;
	v_acoustic = grad*s.u_acoustic;

	%--- Setup plots ---
	fhandleP = setupPlot();
	cmap;
	Sur_el_P = multiblock.Surface(g_elastic, s.p_elastic(:,1));
	shading interp
	hold on;
	Sur_ac_P = multiblock.Surface(g_acoustic, s.p_acoustic(:,1));
	h_P = gca;
    adjustPlot(fhandleP, fontsize)
    colorbar

    fhandleVx = setupPlot();
	cmap;
	Sur_el_Vx = multiblock.Surface(g_elastic, s.ut_elastic(1:2:end-1,1));
	shading interp
	hold on;
	Sur_ac_Vx = multiblock.Surface(g_acoustic, v_acoustic(1:2:end-1,1));
	h_Vx = gca;
    adjustPlot(fhandleVx, fontsize)
    colorbar

    fhandleVy = setupPlot();
	cmap;
	Sur_el_Vy = multiblock.Surface(g_elastic, s.ut_elastic(2:2:end,1));
	shading interp
	hold on;
	Sur_ac_Vy = multiblock.Surface(g_acoustic, v_acoustic(2:2:end,1));
	h_Vy = gca;
    adjustPlot(fhandleVy, fontsize)
    colorbar
    % -----------------

    [~, N] = size(s.p_elastic);
	for i = 1:N;
		t = T(i);
		% if t > Tend
		% 	break;
		% end
		updatePlotElAc(h_P, Sur_el_P, Sur_ac_P, s.p_elastic(:,i), s.p_acoustic(:,i), t, climsP, Pstr);
		drawnow;
		updatePlotElAc(h_Vx, Sur_el_Vx, Sur_ac_Vx, s.ut_elastic(1:2:end-1,i), v_acoustic(1:2:end-1,i), t, climsVx, Vxstr);
		drawnow;
		updatePlotElAc(h_Vy, Sur_el_Vy, Sur_ac_Vy, s.ut_elastic(2:2:end,i), v_acoustic(2:2:end,i), t, climsVy, Vystr);
		drawnow;

		if saveMovie
			%===== Add frame to movie ===%
	        frame = getframe(fhandleP);
	        writeVideo(wObjPressure,frame);

	        frame = getframe(fhandleVx);
	        writeVideo(wObjVx,frame);

	        frame = getframe(fhandleVy);
	        writeVideo(wObjVy,frame);
	        %===========================%
	    end
	end

	if saveMovie
		%== Close movie objects ==%
		close(wObjPressure);
		close(wObjVx);
		close(wObjVy);
		%========================%
	end
	%------------------------

end

function updatePlotElAc(h, Sur_el, Sur_ac, p_el, p_ac, t, cmax, fieldStr)
	Sur_el.ZData = p_el/cmax;
    Sur_el.CData = p_el/cmax;
    Sur_ac.ZData = p_ac/cmax;
    Sur_ac.CData = p_ac/cmax;
    cbh = colorbar;
    caxis(h, [-1, 1]);
    set(cbh, 'YTick', -1:1);
    title(h, [fieldStr '   t = ' sprintf('%4.2f',t)]);
    % xlim([0, 17000])
    % ylim([-3500, 0])

    xlim([9000, 12000])
    ylim([-2000, 0])
    % axis equal
end

function f = setupPlot()
	scrsz = get(0,'ScreenSize');
    f = figure('Position',[0.05*scrsz(3) 0.05*scrsz(4) 0.75*scrsz(3) 2/3*0.75*scrsz(3)]);

    tint = 0.0;
    bgcolor = [1, 1, 1];
    set(gcf, 'Color', bgcolor );
    set(gca, 'Color', bgcolor, 'Xcolor',[tint,tint,tint], ...
        'Ycolor',[tint,tint,tint]);
end

function adjustPlot(fhandle, fontsize)
	xlabel('distance (km)')
    ylabel('depth (km)')
    shading interp
    setFontSize(fhandle, fontsize);

    % axis equal
    xlim([9000, 12000])
    ylim([-2000, 0])
    yticks([-2000:400:0]);
    yticklabels({'2', '1.6' '1.2' '0.8' '0.4' '0'})
    xticks([9000:1000:12000]);
    xticklabels({'9', '10', '11', '12'})
    line([9000,12000],[-449.6,-449.6],[10,10],'Color', 'black', 'linewidth',1);
end

function s = loadFile(filename)
	load(filename);
	s = saveData;
end




