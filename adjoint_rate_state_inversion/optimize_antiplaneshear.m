function [parsVec, misFit, flag, output, history, initialPars, truePars, parSet] = optimize_antiplaneshear(maxIter, m, m_p, misfitType, optPars, val, scalingStudy, loadData, horizontalFault, loadOptPath, lbfgs)
    
    addSubpaths;
    default_arg('scalingStudy',true);
    default_arg('loadData',false);
    default_arg('horizontalFault',false);
    default_arg('loadOptPath','');
    default_arg('lbfgs',true);
    order = 8;
    optObj = @AntiplaneShear2DRSFrictionOpt;
    
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
    [parsVec, misFit, flag, output, history, initialPars, truePars, parSet] = adjopt.optimize_lbfgs(optObj, parSetFun, optPars, parSetOpts, order, maxIter, loadDataPath, loadOptPath, lbfgs);

    if lbfgs
        opt_str = 'lbfgs';
    else
        opt_str = 'bfgs';
    end

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
    filename = sprintf('mat/inversion/%s_%s_maxiter%d_m%d_mp%d_%s_%s_%s.mat',fault_str,opt_str,maxIter,m,m_p,misfitType,data_str,str_ending);
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
    save(filename, '-v7.3','-nocompression');
end