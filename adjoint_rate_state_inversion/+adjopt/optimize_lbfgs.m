% Optimize optObjFunc using the parameters specified in parSetFunc
%
% optObjFunc: 	Handle to an AdjointOptimization class
% parSetFunc: 		Handle to a parameter definition function
function [parsVec, misFit, flag, output, history, initialPars, truePars, parset] = optimize_lbfgs(optObjFunc, parSetFunc, invPars, parSetOpts, order, maxIter, loadDataPath, loadOptPath, lbfgs)
    default_arg('maxIter', 20);
    default_arg('order', 4);
    default_arg('parSetOpts', struct);
    default_arg('loadDataPath', '');
    default_arg('loadOptPath', '');
    default_arg('lbfgs', true);

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

    % Set up misfit history
    history.x = [];
    history.fval = [];
    
    fun = @(parsVec) adjopt.vectorGradient(parsVec, optObj);
    
    % initialize solution and set history of loading
    if ~strcmp(loadOptPath,'')
        IC = load(loadOptPath);
        history = IC.history;
        x0 = IC.history.x(:,end);
        [~, histIter] = size(IC.history.x);
        maxIter = maxIter - (histIter-3);
        % Set loaded parameters in optimization object
        optObj.setParameters(optObj.parametersVectorToNative(x0));
        optObj.updateForwardDiscr();
    else
        x0 = optObj.parametersNativeToVector( optObj.pars );
    end

    if lbfgs
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
        lb = [];
        ub = [];
        nonlcon = [];

        % Optimize
        tic
        [parsVec, misFit, flag, output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        toc

    else
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
    end

    
    % 'output' function that tells fmincon to store the objective function at each iteration
    function stop = outputMisfit(x, optimValues, state)
        stop = false;

        % Store misfit (fval = value of objective function)
        fval = optimValues.fval;
        history.x = [history.x, x];
        history.fval = [history.fval; fval];
    end
end


