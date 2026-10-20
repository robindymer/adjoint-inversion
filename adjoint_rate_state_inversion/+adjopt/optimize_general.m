% Optimize optObjFunc using the parameters specified in parSetFunc
%
% optObjFunc: 	Handle to an AdjointOptimization class
% parSetFunc: 		Handle to a parameter definition function
% NOTE: Previously called optimize_lbfgs, now we call it optimize_general
function [parsVec, misFit, flag, output, history, initialPars, truePars, parset] = optimize_general(optObjFunc, parSetFunc, invPars, parSetOpts, order, maxIter, loadDataPath, loadOptPath, opt_routine_str)
    default_arg('maxIter', 20);
    default_arg('order', 4);
    default_arg('parSetOpts', struct);
    default_arg('loadDataPath', '');
    default_arg('loadOptPath', '');
    default_arg('opt_routine_str','lbfgs');

    % Plot settings
    if ~strcmp(loadDataPath,'') % Load data and create optimization object
        [parset, ~] = parSetFunc(parSetOpts);
        fprintf('Loading data from %s\n',loadDataPath);
        S = load(loadDataPath);
        nr = numel(parset.receivers.x);
        receiverRecordings = cell(nr,1);
        for i = 1:nr
            for j = 1:numel(S.receiverData.positions)
                if S.receiverData.positions{j} == parset.receivers.x{i}
                    break
                end
            end
            receiverRecordings{i} = S.receiverData.recordings{j};
        end
        % Create optimization object and set data
        optObj = optObjFunc(parset, invPars, order);
        optObj.receiverData = receiverRecordings;

        initialPars = optObj.pars;
        truePars = [];
    else % Run with true parameters, save data and create optimization object
        [parset, trueParset] = parSetFunc(parSetOpts);
        disp('Generating synthetic data')
        trueOptObj = optObjFunc(trueParset, invPars, order);
        trueOptObj.runForward(false);
        trueReceiverData = trueOptObj.forwardReceiverRecordings;
        trueTimeIntegrationData = trueOptObj.forwardTimeIntegrationData;
        k = trueOptObj.k;
        truePars = trueOptObj.pars;
        clear trueOptObj;

        % Create optimization object and set data
        optObj = optObjFunc(parset, invPars, order, k);
        optObj.setSyntheticReceiverData(trueReceiverData, trueTimeIntegrationData);
        initialPars = optObj.pars;
    end

    % Set up misfit history. history.t is the wall time (s) since the optimizer started,
    % so solvers can be compared by time rather than by iterations.
    history.x = [];
    history.fval = [];
    history.t = [];
    
    fun = @(parsVec) adjopt.vectorGradient(parsVec, optObj);
    
    % initialize solution and set history of loading
    if ~strcmp(loadOptPath,'')
        IC = load(loadOptPath);
        history = IC.history;
        x0 = IC.history.x(:,end);
        [~, histIter] = size(IC.history.x);
        if isfield(history, 't')
            maxIter = maxIter - (histIter-1); % One entry per iteration, including iteration 0
        else
            % Older histories also stored fmincon's 'init' and 'done' states and have no time
            maxIter = maxIter - (histIter-3);
            history.t = nan(histIter, 1);
        end
        % Set loaded parameters in optimization object
        optObj.setParameters(optObj.parametersVectorToNative(x0));
        optObj.updateForwardDiscr();
    else
        x0 = optObj.parametersNativeToVector( optObj.pars );
    end
    tOpt = tic; % Start of optimization, for history.t

    if strcmp(opt_routine_str,'lbfgs')
        % LBFGS via MATLABs fmincon
        options = optimoptions('fmincon');
        options.SpecifyObjectiveGradient = true;
        options.Algorithm = 'interior-point';
        options.HessianApproximation = 'lbfgs';
        options.MaxIterations = maxIter;
        options.OutputFcn = @outputMisfit;
        options.Display = 'iter-detailed';
        options.OptimalityTolerance = 1e-10;
        %options.PlotFcn = {@optimplotx,@optimplotfval};

        % fmincon arguments
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        % No bound: the interior-point barrier for a > 0 stalled L-BFGS (notes/optimizer_comparison.md)
        lb = [];
        ub = [];
        nonlcon = [];

        % Optimize
        tic
        [parsVec, misFit, flag, output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        toc
    elseif strcmp(opt_routine_str,'bfgs')
        % BFGS via MATLABs fminunc
        options = optimoptions('fminunc');
        options.SpecifyObjectiveGradient = true;
        options.Algorithm = 'quasi-newton';
        options.HessianApproximation = 'bfgs';
        options.MaxIterations = maxIter;
        options.OutputFcn = @outputMisfit;
        options.Display = 'iter-detailed';
        options.OptimalityTolerance = 1e-10;
        options.PlotFcn = {@optimplotx,@optimplotfval};

        % Optimize
        tic
        [parsVec, misFit, flag, output] = fminunc(fun,x0,options);
        toc
    elseif strcmp(opt_routine_str,'CG')
        % TODO: Implement CG optimization routine
        % LBFGS via MATLABs fmincon
        options = optimoptions('fmincon');
        options.SpecifyObjectiveGradient = true;
        options.Algorithm = 'interior-point';
        % options.HessianApproximation = 'lbfgs';
        options.SubproblemAlgorithm = 'cg';
        options.HessianMultiplyFcn = @(x, lambda, v) adjopt.hessianVector(x, v, optObj);
        options.MaxIterations = maxIter;
        options.OutputFcn = @outputMisfit;
        options.Display = 'iter-detailed';
        options.OptimalityTolerance = 1e-10;
        %options.PlotFcn = {@optimplotx,@optimplotfval};

        % fmincon arguments
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = []; % No bound, as for 'lbfgs'
        ub = [];
        nonlcon = [];

        % Optimize
        tic
        [parsVec, misFit, flag, output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        toc
    elseif any(strcmp(opt_routine_str,{'newton-tr','newton-tr-noprec'}))
        % Trust-region Newton-CG with Hessian-vector products and L-BFGS preconditioning.
        % 'newton-tr-noprec' is the same without preconditioner, i.e. plain Steihaug CG.
        hessVec = @(x, V) adjopt.hessianVector(x, V, optObj);
        tropts.maxIter = maxIter;
        tropts.precond = strcmp(opt_routine_str,'newton-tr');
        tropts.OutputFcn = @(x, f, g, info) outputMisfit(x, struct('fval', f), 'iter');
        % Rate-and-state parameters are positive.
        lb = zeros(size(x0));

        % Optimize
        tic
        [parsVec, misFit, flag, output] = adjopt.trustRegionNewtonCG(fun, hessVec, x0, lb, tropts);
        toc
    elseif strcmp(opt_routine_str,'fminunc-tr')
        % MATLAB's trust-region Newton-CG (Steihaug) with Hessian-vector products, no bounds.
        % Hinfo (third output of the objective) is the current x, passed to HessianMultiplyFcn.
        options = optimoptions('fminunc');
        options.Algorithm = 'trust-region';
        options.SpecifyObjectiveGradient = true;
        options.HessianMultiplyFcn = @(Hinfo, Y) adjopt.hessianVector(Hinfo, Y, optObj);
        options.MaxIterations = maxIter;
        options.OutputFcn = @outputMisfit;
        options.Display = 'iter-detailed';
        options.OptimalityTolerance = 1e-10;

        % Optimize
        tic
        [parsVec, misFit, flag, output] = fminunc(@misfitGradientHinfo, x0, options);
        toc
    elseif strcmp(opt_routine_str,'newton-cg-ls')
        % Line-search Newton-CG with Armijo backtracking by M. Heinkenschloss (matlab_impl_constr),
        % used as is. It calls fval, grad, Hessvec, xnew and xprod by name, so the adapters in
        % heinkenschloss_adapter/ are put on the path only for this solve.
        baseDir = fileparts(fileparts(mfilename('fullpath')));
        adapterPaths = {fullfile(fileparts(baseDir), 'matlab_impl_constr', 'optimization'), ...
                        fullfile(baseDir, 'heinkenschloss_adapter')};
        addpath(adapterPaths{:});
        cleanupPaths = onCleanup(@() rmpath(adapterPaths{:}));

        usr_par.optObj = optObj;
        usr_par.outputFcn = @(x, f) outputMisfit(x, struct('fval', f), 'iter');
        ncgopts.iprint = 1;
        ncgopts.maxit = maxIter;
        ncgopts.gtol = 1e-4;
        ncgopts.stol = 1e-14;

        % Optimize
        tic
        [parsVec, iter, flag] = newton_cg(x0, ncgopts, usr_par);
        toc
        misFit = history.fval(end);
        output = struct('iterations', iter);
    else
        error('Unknown optimization routine: %s', opt_routine_str);
    end

    
    % 'output' function that tells fmincon to store the objective function at each iteration.
    % Only 'iter' is stored (it includes iteration 0), not 'init' or 'done', which repeat iterates.
    function stop = outputMisfit(x, optimValues, state)
        stop = false;
        if ~strcmp(state, 'iter')
            return
        end

        % Store misfit (fval = value of objective function)
        fval = optimValues.fval;
        history.x = [history.x, x];
        history.fval = [history.fval; fval];
        history.t = [history.t; toc(tOpt)];
    end

    % Objective for fminunc 'trust-region', which passes its third output to HessianMultiplyFcn
    function [f, g, Hinfo] = misfitGradientHinfo(x)
        [f, g] = fun(x);
        Hinfo = x;
    end
end


