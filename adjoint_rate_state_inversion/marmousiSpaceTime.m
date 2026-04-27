function marmousiSpaceTime(filename, downSampling, saveFig,...
									climsP, climsVx, climsVy)
	default_arg('filename', 'elAc')
	default_arg('downSampling', 40);
	default_arg('saveFig', false);

	default_arg('climsP', 6e-4);
	default_arg('climsVx', 1e-9);
	default_arg('climsVy', 1e-9);

	Pstr = 'Pressure';
	Vxstr = 'v_x';
	Vystr = 'v_y';

	fontsize = 16;

	if saveFig
		switch filename
		case 'CmapElSeafloor10p7'
			filenameP = ['marmousiSpaceTime/elasticSeafloor_' 'Pressure'];
			filenameVx = ['marmousiSpaceTime/elasticSeafloor_' 'Vx'];
			filenameVy = ['marmousiSpaceTime/elasticSeafloor_' 'Vy'];
		otherwise
			filenameP = ['marmousiSpaceTime/' filename 'Pressure'];
			filenameVx = ['marmousiSpaceTime/' filename 'Vx'];
			filenameVy = ['marmousiSpaceTime/' filename 'Vy'];
		end
	end

	filename = [filename, 'DS', num2str(downSampling)];

	s = loadFile([filename, '.mat']);

	T = s.t;
	g_elastic = s.g.elastic;

	Tfinal = 3;
	p = g_elastic.points();
	X = g_elastic.funcToPlotMatrices(p(:,1));
	X = X{1};
	x = X(end,:);
    t = Tfinal-s.t;
    [X,T] = meshgrid(x,t);

	[~, N] = size(s.p_elastic);
	vx = [];
	vy = [];
	p = [];
	for i = 1:N;
		t = T(i);
		p_i = g_elastic.funcToPlotMatrices(s.p_elastic(:,i));
		p_i = p_i{1};
		p = [p_i(end, :); p];

		vx_i = g_elastic.funcToPlotMatrices(s.ut_elastic(1:2:end-1,i));
		vx_i = vx_i{1};
		vx = [vx_i(end, :); vx];

		vy_i = g_elastic.funcToPlotMatrices(s.ut_elastic(2:2:end,i));
		vy_i = vy_i{1};
		vy = [vy_i(end, :); vy];
	end

	%---Plot ---
	fhandleP = setupPlot();
	plotRecordSection(fhandleP, X, T, p, climsP, fontsize, Pstr);
	if saveFig
		print(filenameP, '-depsc');
	end

	fhandleVx = setupPlot();
	plotRecordSection(fhandleVx, X, T, vx, climsVx, fontsize, Vxstr);
	if saveFig
		print(filenameVx, '-depsc');
	end

	fhandleVy = setupPlot();
	plotRecordSection(fhandleVy, X, T, vy, climsVy, fontsize, Vystr);
	if saveFig
		print(filenameVy, '-depsc');
	end


 %    fhandleVx = setupPlot();
	% cmap;
	% Sur_el_Vx = multiblock.Surface(g_elastic, s.ut_elastic(1:2:end-1,1));
	% shading interp
	% hold on;
	% Sur_ac_Vx = multiblock.Surface(g_acoustic, v_acoustic(1:2:end-1,1));
	% h_Vx = gca;
 %    adjustPlot(fhandleVx, fontsize)
 %    colorbar

 %    fhandleVy = setupPlot();
	% cmap;
	% Sur_el_Vy = multiblock.Surface(g_elastic, s.ut_elastic(2:2:end,1));
	% shading interp
	% hold on;
	% Sur_ac_Vy = multiblock.Surface(g_acoustic, v_acoustic(2:2:end,1));
	% h_Vy = gca;
 %    adjustPlot(fhandleVy, fontsize)
 %    colorbar

end


function f = setupPlot()
	scrsz = get(0,'ScreenSize');
    f = figure('Position',[0.05*scrsz(3) 0.05*scrsz(4) 0.75*scrsz(3) 2/3*0.75*scrsz(3)]);

    tint = 0.0;
    bgcolor = [1, 1, 1];
    set(gcf, 'Color', bgcolor );
    set(gca, 'Color', bgcolor, 'Xcolor',[tint,tint,tint], ...
        'Ycolor',[tint,tint,tint]);
    grid off;
end

function s = loadFile(filename)
	load(filename);
	s = saveData;
end

function plotRecordSection(fhandle, X, T, v, cmax, fontsize, titleStr)
	cmap;
	surf(X,-T,v/cmax);
	xlabel('offset (km)')
    ylabel('time (s)')
    view(0,90);
    shading interp
    grid off
    setFontSize(fhandle, fontsize);

    cbh = colorbar;
    h = gca;
    caxis(h, [-1, 1]);
    set(cbh, 'YTick', -1:1);

    % axis equal
    xlim([9000, 12000])
    ylim([-3, 0])
    yticks(-3:0);
    yticklabels({'3' '2' '1', '0'})
    xticks(9200:500:11700);
    xticklabels({'-1.5' '-1' '-0.5' '0' '0.5' '1'})
    title(titleStr);
end




