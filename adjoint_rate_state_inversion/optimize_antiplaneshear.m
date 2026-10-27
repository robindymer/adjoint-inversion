function [parsVec, misFit, flag, output, history, initialPars, truePars, parSet] = optimize_antiplaneshear(maxIter, m, m_p, misfitType, optPars, val, scalingStudy, loadData, horizontalFault, loadOptPath, opt_routine_str, receiverSpacing, regOpts)
    
    addSubpaths;
    default_arg('scalingStudy',true);
    default_arg('loadData',false);
    default_arg('horizontalFault',false);
    default_arg('loadOptPath','');
    % lbfgs or CG or bfgs
    % lbfgs and bfgs are quasi-newton methods, CG is a conjugate gradient method
    % lbfgs usually better than bfgs for many unknowns, but require same information (gradient)
    % CG require gradient and hessian-vector product
    % newton-tr is trust-region Newton-CG (newton-tr-noprec without preconditioner), requires gradient and hessian-vector product
    % fminunc-tr is MATLAB's trust-region Newton-CG, newton-cg-ls is M. Heinkenschloss' line-search
    % Newton-CG (matlab_impl_constr). Both require gradient and hessian-vector product
    default_arg('opt_routine_str','lbfgs');
    % 2: N_rec = 88 receivers, 1: N_rec = 325 receivers (see pars.rsFriction2DFractalFaultInversion)
    default_arg('receiverSpacing',2);
    % Regularization added to the misfit, see adjopt.regularization. E.g. struct('type','tv','beta',1)
    default_arg('regOpts',struct('type','none'));
    order = 8;
    if any(strcmp(opt_routine_str, {'CG', 'newton-tr', 'newton-tr-noprec', 'fminunc-tr', 'newton-cg-ls'})) % Requires Hessian-vector products
        optObj = @AntiplaneShear2DRSFrictionOptHessian;
    else
        optObj = @AntiplaneShear2DRSFrictionOpt;
    end
    
    if loadData
         % Generate this data by running fractal_fault_generate_data.m with the 
         % pars.rsFriction2DFractalFaultForward parameter set.
        assert(isequal(misfitType,'velocity'));
        data_str = 'highres';
        loadDataPath = 'mat/fractal_fault_m1001_mp1001/receiverData.mat';
    else
        loadDataPath = '';
        data_str = 'synth';
    end

    if horizontalFault
        fault_str = 'horiz';
        parSetFun = @pars.rsFriction2DHorizontalFaultInversion;
    else
        fault_str = 'rough';
        parSetFun = @pars.rsFriction2DFractalFaultInversion;
    end

    if scalingStudy
        parSetOpts.initialGuessScalings = val;
        parSetOpts.initialGuessValues = [];
    else
        parSetOpts.initialGuessScalings = [];
        parSetOpts.initialGuessValues = val;
    end
    parSetOpts.misfitType = misfitType;
    parSetOpts.inversionParameters = optPars;
    parSetOpts.m = m;
    parSetOpts.m_p = m_p;
    parSetOpts.receiverSpacing = receiverSpacing;

    % if lbfgs
    %     opt_str = 'lbfgs';
    % else
    %     opt_str = 'bfgs';
    % end
    if scalingStudy
        scalings = parSetOpts.initialGuessScalings;
        str_ending = 'scaled';
        for i = 1:length(scalings)
            parname = optPars{i};
            s = scalings(i);
            str_ending = erase(sprintf('%s_%s%1.2f',str_ending,parname,s),'.');
        end
    else
        vals = parSetOpts.initialGuessValues;
        str_ending = 'vals';
        for i = 1:length(vals)
            parname = optPars{i};
            p = vals(i);
            str_ending = erase(sprintf('%s_%s%1.4f',str_ending,parname,p),'.');
        end
    end
    if receiverSpacing ~= 2
        str_ending = sprintf('%s_rs%d',str_ending,receiverSpacing);
    end
    if ~strcmp(regOpts.type,'none')
        % %.0e, since the file name may not contain '.' (see the file counter below)
        str_ending = sprintf('%s_%s_beta%.0e',str_ending,regOpts.type,regOpts.beta);
        if any(strcmp(regOpts.type,{'tv','lasso'})) && isfield(regOpts,'eps')
            str_ending = sprintf('%s_eps%.0e',str_ending,regOpts.eps);
        end
    end
    filename = sprintf('mat/inversion/%s_%s_maxiter%d_m%d_mp%d_%s_%s_%s.mat',fault_str,opt_routine_str,maxIter,m,m_p,misfitType,data_str,str_ending);
    if isfile(filename)
        filecounter = 1;
        split_str = split(filename,'.');
        newfilename = [split_str{1},'_',num2str(filecounter),'.',split_str{2}];
        while isfile(newfilename)
            filecounter = filecounter+1;
            newfilename = [split_str{1},'_',num2str(filecounter),'.',split_str{2}];
        end
        filename = newfilename;
    end
    mkdir('mat/inversion/');
    % history is saved here after every iteration, and the file is removed once the run has finished
    checkpointPath = strrep(filename, '.mat', '_checkpoint.mat');
    [parsVec, misFit, flag, output, history, initialPars, truePars, parSet] = adjopt.optimize_general(optObj, parSetFun, optPars, parSetOpts, order, maxIter, loadDataPath, loadOptPath, opt_routine_str, checkpointPath, regOpts);

    save(filename, '-v7.3','-nocompression');
    delete(checkpointPath);
end