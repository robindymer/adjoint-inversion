% Optimize optObjFunc using the parameters specified in parSetFunc
%
% optObjFunc: 	Handle to an AdjointOptimization class
% parSetFunc: 		Handle to a parameter definition function
function optimize(optObjFunc, parSetFunc, invPars, parSetOpts, order, optMethod, maxIter, saveMovie, movieName)
	default_arg('movieName','Movie');
	default_arg('saveMovie', false);
	default_arg('maxIter', 20);
	default_arg('optMethod', 'goldensection');
	default_arg('order', 2);
	default_arg('parSetOpts', struct);

	% Plot settings
	fontsize = 16;

	[parset, trueParset] = parSetFunc(parSetOpts);
	% nReceivers = numel(parset.receivers.x);
	% nSources = numel(parset.sources.x);
	dim = parset.dim;


	% Run with true parameters and save data
	trueOptObj = optObjFunc(trueParset, invPars, order);
	trueOptObj.runForward(false);
	trueReceiverData = trueOptObj.forwardReceiverRecordings;
	trueTimeIntegrationData = trueOptObj.forwardTimeIntegrationData;
	k = trueOptObj.k;
	trueParset = trueOptObj.pars;
	clear trueOptObj;

	% Create optimization object
	optObj = optObjFunc(parset, invPars, order, k);
	optObj.setSyntheticReceiverData(trueReceiverData,trueTimeIntegrationData);
	parset = optObj.pars;

	% Compute initial misfit
	optObj.runForward(false);
	M0 = optObj.computeMisfit();

	% Setup movie file
	if saveMovie
		movieName = [movieName '.avi'];
		writerObj = VideoWriter(movieName);
		writerObj.FrameRate = 2;
		open(writerObj);
	end

	M = zeros(maxIter,1);

	tic
	switch optMethod
	case 'goldensection'

		% Setup figure
		[upPar, upRec, upMis, figure_handle] = optObj.setupPlot(trueParset, parset, M, M0, maxIter);

		% Golden section parameters
		tol = 1e-2;

		for i = 1:maxIter
			% Compute gradient
			grad = optObj.computeGradient();
			M(i) = optObj.computeMisfit();

			upMis(M);
			upRec(optObj);
			upPar(optObj);

			subplot(3,1,3);
			title(sprintf('iter: %d', i-1), 'fontsize', fontsize);
			drawnow;

			if saveMovie
				%===== Add frame to movie ===%
		        frame = getframe(gcf);
		        writeVideo(writerObj,frame);
		        %===========================%
		    end

			% One golden section iteration
			relnorm = optObj.gradientNorm(grad);
			L = 0.25*1/relnorm;
			optObj = adjopt.goldenSection(optObj, grad, L, tol);
		end

		optObj.runForward();
		M(end) = optObj.computeMisfit();
		upMis(M);

	% LBFGS via MATLABs fmincon
	case 'lbfgs'
		options = optimoptions('fmincon');
		options.SpecifyObjectiveGradient = true;
		options.Algorithm = 'interior-point';
		options.HessianApproximation = 'lbfgs';
		options.MaxIterations = maxIter;

		% Set up misfit history
		history.fval = [];
		%options.OutputFcn = @outputMisfit;
        options.PlotFcn = {@optimplotx,@optimplotfval};
		options.Display = 'iter-detailed';

		fun = @(parsVec) adjopt.vectorGradient(parsVec, optObj);
		x0 = optObj.parametersNativeToVector( optObj.pars );

		% fmincon arguments
		A = [];
		b = [];
		Aeq = [];
		beq = [];
		lb = [];
		ub = [];
		nonlcon = [];

		% Setup movie figure
		%[upPar, upRec, upMis, figure_handle] = optObj.setupPlot(trueParset, parset, M, M0, maxIter);

		% Optimize
		disp('');
		str = util.replace_string('','Iteration %d',0);
		[parsVec, ~, flag, output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
		str = util.replace_string(str,'');

		% Plot final
		%optObj.setupPlot(trueParset, optObj.pars, history.fval, M0, length(history.fval));
        parsVec

	end

	toc
	
	% 'output' function that tells fmincon to store the objective function at each iteration
	function stop = outputMisfit(x, optimValues, state)
		stop = false;

		% Store misfit (fval = value of objective function)
		fval = optimValues.fval;
		history.fval = [history.fval; fval];
		it = optimValues.iteration;

		if it >= 1
			M(it) = fval;
			upMis(M);
			upRec(optObj);
			upPar(optObj);
			drawnow;
		end

		if saveMovie
	        frame = getframe(gcf);
	        writeVideo(writerObj,frame);
		end

		% Print iteration number
		str = util.replace_string(str,'Iteration %d', it);

	end

	if saveMovie
		%== Close movie object ==%
		close(writerObj);
		%========================%
	end

end


