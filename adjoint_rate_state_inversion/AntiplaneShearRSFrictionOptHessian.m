% Adjoint optimization for anti-plane shear in 1D with a frictional
% fault interface

classdef AntiplaneShearRSFrictionOptHessian < adjopt.AdjointOptimization

properties
    % Sources
    adjointReceiverRecordings 	%Recorded during LATEST adjoint simulation (at forward source coords).
    %Stored with FORWARD time convention, like everything else.
    secondOrderAdjointReceiverRecordings
    secondOrderForwardReceiverRecordings
    
    % Receivers
    receiverData 				%True data, recorded by "seismometers"
    forwardReceiverRecordings	%Recorded during LATEST forward simulation
    misfitType
    
    % Discretizations
    forwardDiscr
    adjointDiscr
    secondOrderForwardDiscr
    secondOrderAdjointDiscr
    
    % Time stepping options
    tsOpts
    
    %Data recorded during LATEST forward simulation
    forwardFaultVariables	   
    forwardTimeIntegrationData 
    
    %Data recorded during LATEST adjoint simulation
    adjointFaultVariables	    
    adjointTimeIntegrationData

    %Data recorded during LATEST second order forward simulation
    secondOrderForwardFaultVariables
    secondOrderForwardTimeIntegrationData

    %Data recorded during LATEST second order forward simulation
    secondOrderAdjointFaultVariables
    secondOrderAdjointTimeIntegrationData
    
    pars			%Current parameters.
    T 				%Final time
    H               % Space quadrature
    Ht              % Time quadrature for LATEST forward simulation
    
    m				%Grid points
    order
    dim
    k
    
    inversionPars
    F_p
    G_p
    eps_pertubations % Parameter pertubation, required for hessian vector computation
    supportedPars
    
end

methods
    
    function obj = AntiplaneShearRSFrictionOptHessian(pars, inversionPars, order, k)
        % Use convention that this represents H_ab, so first a pertubation and then b
        default_arg('inversionPars', {'a', 'b'});
        default_arg('order', 4);
        default_arg('k', []);
        
        % Unpack parameter struct
        opset = pars.opset;
        m = pars.m;
        xlims = pars.xlims;
        material = pars.material;
        bc = pars.bc;
        friction = pars.friction;
        sources = pars.sources;
        receivers = pars.receivers;
        secondOrderSources = pars.secondOrderSources;
        secondOrderReceivers = pars.secondOrderReceivers;
        initialconditions = pars.initialconditions;
        interp_data = pars.interpolate_data; % TODO: ???????? Maybe not necessary right now
        
        tsOpts = pars.tsOpts;
        T = pars.T;
        
        % Create forward discretization object
        forwardDiscr = elastic.AntiplaneShearRSFrictionFwdDiscr(opset, m, xlims, order, material, bc, friction, sources, initialconditions);
        disp(initialconditions)
        % Initialize adjoint discretization object. Does not set any friction/receiver data.
        adjointDiscr = elastic.AntiplaneShearRSFrictionAdjDiscr(opset, m, xlims, order, material, bc, friction, receivers, interp_data);
        disp(interp_data)
        % Same thing here, do not set friction and second order source / reciever data
        secondOrderForwardDiscr = elastic.AntiplaneShearRSFrictionSecondOrderFwdDiscr(opset, m, xlims, order, material, bc, friction, secondOrderSources, interp_data);
        secondOrderAdjointDiscr = elastic.AntiplaneShearRSFrictionSecondOrderAdjDiscr(opset, m, xlims, order, material, bc, friction, secondOrderReceivers, interp_data);

        pars.a = pars.friction.params.a;
        pars.b = pars.friction.params.b;

        obj.eps_pertubations = struct();
        obj.F_p = cell(numel(inversionPars), 1);
        obj.G_p = cell(numel(inversionPars), 1);
        selectDerivativeFun = str2func('pars.rsGetPartialDerivativeFun');
        for i = 1:numel(inversionPars)
            parName = inversionPars{i};
            [obj.F_p{i}, obj.G_p{i}] = selectDerivativeFun(friction.funs, parName);

            epsName = sprintf('eps_%s', parName);
            if isfield(pars.friction.params, epsName)
                obj.eps_pertubations.(parName) = pars.friction.params.(epsName);
            end
        end
        
        % Set instance variables
        obj.dim = 1;
        obj.forwardDiscr = forwardDiscr;
        obj.adjointDiscr = adjointDiscr; % Note: Must be updated with data prior to calling runAdjoint
        obj.secondOrderForwardDiscr = secondOrderForwardDiscr;
        obj.secondOrderAdjointDiscr = secondOrderAdjointDiscr;
        
        obj.tsOpts = tsOpts;
        
        obj.forwardReceiverRecordings = [];
        obj.adjointReceiverRecordings = [];
        obj.secondOrderForwardReceiverRecordings = [];
        obj.secondOrderAdjointReceiverRecordings = [];
        obj.receiverData = [];
        obj.misfitType = pars.misfitType;
        
        obj.forwardFaultVariables = obj.faultVariablesStruct();
        obj.adjointFaultVariables = obj.faultVariablesStruct();
        obj.secondOrderForwardFaultVariables = obj.faultVariablesStruct();
        obj.secondOrderAdjointFaultVariables = obj.faultVariablesStruct();
        
        obj.forwardTimeIntegrationData = obj.timeIntegrationStruct();
        obj.adjointTimeIntegrationData = obj.timeIntegrationStruct();
        obj.secondOrderForwardTimeIntegrationData = obj.timeIntegrationStruct();
        obj.secondOrderAdjointTimeIntegrationData = obj.timeIntegrationStruct();

        
        obj.m = m;
        obj.order = order;
        obj.T = T;
        obj.pars = pars;
        obj.inversionPars  = inversionPars;
        obj.H = forwardDiscr.H;
        obj.Ht = [];
        obj.k = k;
        obj.supportedPars = {'a', 'b'};
    end
    
    function runForward(obj, plotFlag, T, saveOpts, progressBar)
        default_arg('plotFlag', false);
        default_arg('T', obj.T);
        default_arg('saveOpts', []);
        default_arg('progressBar', []);
        
        discr = obj.forwardDiscr;
        % Receiver dirac deltas for measuring misfit
        receiverDeltas = obj.adjointDiscr.dirac_deltas;
        % Time stepping options
        fwdTsOpts.method = obj.tsOpts.forwardMethod;
        fwdTsOpts.T = T;
        fwdTsOpts.k = obj.tsOpts.k;
        fwdTsOpts.cont_time = true;

        if fwdTsOpts.method.adaptive
            [receiverData, faultData, timeData] = obj.runSimulationAdaptive(discr, fwdTsOpts,...
            receiverDeltas, plotFlag, saveOpts, progressBar);
        else
            [receiverData, faultData, timeData] = obj.runSimulation(discr, fwdTsOpts,...
            receiverDeltas, plotFlag, saveOpts, progressBar);
        end
        obj.forwardReceiverRecordings = receiverData;
        obj.forwardFaultVariables = faultData;
        obj.forwardTimeIntegrationData = timeData;
        obj.Ht = timeData.Ht;
        if ~iscolumn(obj.Ht)
            obj.Ht = transpose(obj.Ht);
        end
    end
    
    function runAdjoint(obj, plotFlag, T, saveOpts, progressBar)
        default_arg('plotFlag', false);
        default_arg('T', obj.T);
        default_arg('saveOpts', []);
        default_arg('progressBar', []);
        
        discr = obj.adjointDiscr;
        % Receiver dirac deltas for measuring misfit
        receiverDeltas = obj.forwardDiscr.dirac_deltas;
        
        % Time stepping options
        adjTsOpts.method = obj.tsOpts.adjointMethod;
        adjTsOpts.T = T;
        if obj.adjointDiscr.interpolate_data
            adjTsOpts.cont_time = true;
        else
            adjTsOpts.cont_time = false;
        end

        if adjTsOpts.method.adaptive 
            assert(obj.adjointDiscr.interpolate_data, 'Interpolate data must be set to true.');
            adjTsOpts.k = obj.tsOpts.k;
            [receiverData, faultData, timeData] = obj.runSimulationAdaptive(discr, adjTsOpts,...
            receiverDeltas, plotFlag, saveOpts, progressBar);
        else  % Time integrate using a standard RK method, but with time steps taken from forward problem
            adjTsOpts.k = fliplr(obj.forwardTimeIntegrationData.k); % Reverse timestep vector
            [receiverData, faultData, timeData] = obj.runSimulation(discr, adjTsOpts,...
            receiverDeltas, plotFlag, saveOpts, progressBar);
        end
        
        % TBD: Is this needed?									
        nReceivers = numel(receiverDeltas);
        for i = 1:nReceivers
            receiverData{i} = rot90(receiverData{i});
        end
        obj.adjointReceiverRecordings = receiverData;
        obj.adjointFaultVariables = faultData;
        obj.adjointTimeIntegrationData = timeData;
    end


    function runSecondOrderForward(obj, plotFlag, T, saveOpts, progressBar)
        default_arg('plotFlag', false);
        default_arg('T', obj.T);
        default_arg('saveOpts', []);
        default_arg('progressBar', []);
        
        discr = obj.secondOrderForwardDiscr;
        % Receiver dirac deltas for measuring misfit
        receiverDeltas = obj.secondOrderAdjointDiscr.dirac_deltas;
        
        % Time stepping options
        fwdTsOpts.method = obj.tsOpts.forwardMethod;
        fwdTsOpts.T = T;
        if obj.secondOrderForwardDiscr.interpolate_data
            fwdTsOpts.cont_time = true;
        else
            fwdTsOpts.cont_time = false;
        end

        if fwdTsOpts.method.adaptive
            assert(obj.secondOrderForwardDiscr.interpolate_data, 'Interpolate data must be set to true.');
            fwdTsOpts.k = obj.tsOpts.k;
            [receiverData, faultData, timeData] = obj.runSimulationAdaptive(discr, fwdTsOpts,...
            receiverDeltas, plotFlag, saveOpts, progressBar);
        else
            % Time integrate using a standard RK method, but with time steps taken from forward problem
            fwdTsOpts.k = obj.forwardTimeIntegrationData.k;
            [receiverData, faultData, timeData] = obj.runSimulation(discr, fwdTsOpts,...
            receiverDeltas, plotFlag, saveOpts, progressBar);
        end
        
        obj.secondOrderForwardReceiverRecordings = receiverData;
        obj.secondOrderForwardFaultVariables = faultData;
        obj.secondOrderForwardTimeIntegrationData = timeData;
        obj.Ht = timeData.Ht;
        if ~iscolumn(obj.Ht)
            obj.Ht = transpose(obj.Ht);
        end
    end

    function runSecondOrderAdjoint(obj, plotFlag, T, saveOpts, progressBar)
        default_arg('plotFlag', false);
        default_arg('T', obj.T);
        default_arg('saveOpts', []);
        default_arg('progressBar', []);
        
        discr = obj.secondOrderAdjointDiscr;
        % Receiver dirac deltas for measuring misfit
        receiverDeltas = obj.secondOrderForwardDiscr.dirac_deltas;
        
        % Time stepping options
        adjTsOpts.method = obj.tsOpts.adjointMethod;
        adjTsOpts.T = T;
        if obj.secondOrderAdjointDiscr.interpolate_data
            adjTsOpts.cont_time = true;
        else
            adjTsOpts.cont_time = false;
        end

        if adjTsOpts.method.adaptive 
            assert(obj.secondOrderAdjointDiscr.interpolate_data, 'Interpolate data must be set to true.');
            adjTsOpts.k = obj.tsOpts.k;
            [receiverData, faultData, timeData] = obj.runSimulationAdaptive(discr, adjTsOpts,...
            receiverDeltas, plotFlag, saveOpts, progressBar);
        else  % Time integrate using a standard RK method, but with time steps taken from forward problem
            adjTsOpts.k = fliplr(obj.secondOrderForwardTimeIntegrationData.k); % Reverse timestep vector
            [receiverData, faultData, timeData] = obj.runSimulation(discr, adjTsOpts,...
            receiverDeltas, plotFlag, saveOpts, progressBar);
        end
        
        % TBD: Is this needed?									
        nReceivers = numel(receiverDeltas);
        for i = 1:nReceivers
            receiverData{i} = rot90(receiverData{i});
        end
        obj.secondOrderAdjointReceiverRecordings = receiverData;
        obj.secondOrderAdjointFaultVariables = faultData;
        obj.secondOrderAdjointTimeIntegrationData = timeData;
    end
    
    function [receiverRecordings, faultData, timeData] = runSimulation(obj, discr, tsOpts,...
        receiverDeltas, plotFlag, saveOpts, progressBar)
        default_arg('saveOpts', []);
        default_arg('progressBar', false)
        saveData = obj.saveStruct();
        
        % Setup timestepper ts, vector of timesteps ks, and number of timesteps N needed to
        % reach final time.
        if isempty(tsOpts.k) || numel(tsOpts.k) == 1 % no timestep or single timestep specified
            % Construct ts aligned to final time T
            [ts, N] = discr.getTimestepper(tsOpts.method, tsOpts.T, tsOpts.k);
            ks = ts.k*ones(N,1); % Timestep vector
        else
            ts = discr.getTimestepper(tsOpts.method);
            N = numel(tsOpts.k);
            ks = tsOpts.k;
        end
        [~, ~, ~, nStages] = ts.getTableau();
        
        % For storing data at receivers.
        switch obj.misfitType
        case 'displacement'
            % Record u
            Erec = discr.E.u;
        case 'velocity'
            % Record v
            Erec = discr.E.v;
        end
        
        % Initialize data structs
        nReceivers = numel(receiverDeltas);
        receiverRecordings = cell(nReceivers, obj.dim);
        for i = 1:nReceivers
            receiverRecordings{i} = zeros(1, N*nStages);
        end
        
        faultData = obj.faultVariablesStruct();
        faultData.V = zeros(1, N*nStages); 
        faultData.Psi = zeros(1, N*nStages); 
        
        
        timeData = obj.timeIntegrationStruct();
        timeData.nStages = nStages;
        timeData.T = zeros(1, N*nStages);
        timeData.Ht = zeros(1, N*nStages);
        timeData.k = ks;
        
        % Setup plot
        if plotFlag
            [update, fh] = discr.setupPlot();
            fh();
            r = discr.getTimeSnapshot(0);
            update(r);
            
            if ~isempty(saveOpts) && isfield(saveOpts, 'videoFilename') && ~isempty(saveOpts.videoFilename)
                vidObj = VideoWriter(saveOpts.videoFilename);
                vidObj.FrameRate = 1;
                open(vidObj);
                writeVideo(vidObj, getframe(gcf));
            else
                vidObj = [];
            end
        else
            vidObj = [];
        end
        
        % Initialize progress bar
        if progressBar
            s = util.replace_string('','   %d %%',0);
        end
        
        for i = 1:N
            ts.k = ks(i); % Set timestep and step forward.
            ts.step();
            
            [~, ~, W, T] = ts.getV;
            i_stages = (i-1)*nStages+1 : i*nStages;
            
            % Update progress bars
            if progressBar
                s = util.replace_string(s,'   %.2f %%',i/N*100);
            end
            
            % Store at receivers
            for j = 1:nReceivers		
                receiverRecordings{j}(i_stages) = (obj.H*receiverDeltas{j})'*Erec*W;
            end
            
            % Store fault data and time integration data for the iteration
            switch discr.friction.method
            case 'standard'
                faultData.V(:, i_stages) = discr.fault_jump(discr.E.v*W);
            case 'erickson2022'
                if tsOpts.cont_time
                    for is = 1:nStages
                        faultData.V(:, i_stages(is)) = discr.V_star(T(is), W(:,is));
                    end
                else
                    for is = 1:nStages
                        faultData.V(:, i_stages(is)) = discr.V_star(i_stages(is), W(:,is));
                    end
                end
            end
            faultData.Psi(:, i_stages) = discr.E.Psi*W;
            timeData.T(i_stages) = T;
            timeData.Ht(i_stages) = ts.getTimeStepQuadrature();
            
            % Plot
            if( plotFlag && (mod(i,100) == 0) )
                fh();
                update(discr.getTimeSnapshot(ts));
                drawnow;
                if ~isempty(vidObj)
                    frame = getframe(gcf);
                    writeVideo(vidObj, imresize(frame.cdata, [vidObj.Height, vidObj.Width]));
                end
            end
            
            % Save solution if requested
            if ~isempty(saveOpts) && isfield(saveOpts, 'frameSpacing') && (mod(i, saveOpts.frameSpacing) == 1)
                [w, t] = ts.getV;
                saveData.t = [saveData.t, t];
                saveData.u = [saveData.u, discr.E.u*w];
                saveData.v = [saveData.v, discr.E.v*w];
                saveData.Psi = [saveData.Psi, discr.E.Psi*w];
            end
        end
        
        % End progress bar
        if progressBar
            s = util.replace_string(s,'');
        end
        
        if ~isempty(saveOpts) && isfield(saveOpts, 'filename') && ~isempty(saveOpts.filename)
            saveData.grid = obj.forwardDiscr.grid;
            save(saveOpts.filename, 'saveData', '-v7.3','-nocompression');
        end
        
        if plotFlag && ~isempty(vidObj)
            close(vidObj);
        end
    end
    
    function [receiverRecordings, faultData, timeData] = runSimulationAdaptive(obj, discr, tsOpts,...
        receiverDeltas, plotFlag, saveOpts, progressBar)
        default_arg('saveOpts', []);
        default_arg('progressBar', false)
        saveData = obj.saveStruct();
        
        [ts, N] = discr.getTimestepper(tsOpts.method, tsOpts.T, tsOpts.k);
        [~, b, ~, nStages] = ts.getTableau();
        
        % For storing data at receivers.
        switch obj.misfitType
        case 'displacement'
            % Record u
            Erec = discr.E.u;
        case 'velocity'
            % Record v
            Erec = discr.E.v;
        end
        
        % Initialize data structs
        buffsz = 1e5; % Initial guess of number of timelevels.
        nReceivers = numel(receiverDeltas);
        receiverRecordings = cell(nReceivers, obj.dim);
        for i = 1:nReceivers
            receiverRecordings{i} = zeros(1,buffsz*nStages);
        end
        
        faultData = obj.faultVariablesStruct();
        faultData.V = zeros(1, buffsz*nStages); 
        faultData.Psi = zeros(1, buffsz*nStages); 
        
        timeData = obj.timeIntegrationStruct();
        timeData.T = zeros(1, buffsz*nStages);
        timeData.Ht = zeros(1, buffsz*nStages);
        timeData.nStages = nStages;

        % Setup plot
        if plotFlag
            [update, fh] = discr.setupPlot();
            fh();
            r = discr.getTimeSnapshot(0);
            update(r);
            
            if ~isempty(saveOpts) && isfield(saveOpts, 'videoFilename') && ~isempty(saveOpts.videoFilename)
                vidObj = VideoWriter(saveOpts.videoFilename);
                vidObj.FrameRate = 1;
                open(vidObj);
                writeVideo(vidObj, getframe(gcf));
            else
                vidObj = [];
            end
        else
            vidObj = [];
        end
        
        % Initialize progress bar
        if progressBar
            s = util.replace_string('','   %d %%',0);
        end
        
        t = ts.t;
        i = 1; % Iteration counter
        Nalloc = 1;
        while t < tsOpts.T
            
            if t + ts.k > tsOpts.T % Make sure we don't overshoot final time
                ts.k = tsOpts.T-t;
            end

            % Increase size of data vectors if needed
            if i > Nalloc*buffsz
                % Add additional vector of size buffsz
                tmp = zeros(1,buffsz);
                Tmp = zeros(1,buffsz*nStages);
                
                for j = 1:nReceivers
                    receiverRecordings{j} = [receiverRecordings{j}, Tmp];
                end
                 
                timeData.k = [timeData.k, tmp];
                timeData.T = [timeData.T, Tmp];
                timeData.Ht = [timeData.Ht, Tmp];

                faultData.V = [faultData.V, Tmp];
                faultData.Psi = [faultData.Psi, Tmp];
                
                Nalloc = Nalloc+1;
                clear tmp;
                clear Tmp;
            end                
            
            ts.step();
            
            % Update progress bars
            if progressBar
                s = util.replace_string(s,'   %.2f %%',i/N*100);
            end
            
            [~, t, W, T] = ts.getV; % Solution and times for all stages in previous step
            i_stages = (i-1)*nStages+1 : i*nStages; % Global indices for stage rates
            
            % Store at receivers
            for j = 1:nReceivers			
                receiverRecordings{j}(i_stages) = (obj.H*receiverDeltas{j})'*Erec*W;
            end
            
            
            % Store fault data and time integration data for the iteration
            k = T(end) - T(1); % Extract the timestep used for the stages;
            switch discr.friction.method
            case 'standard'
                faultData.V(:, i_stages) = discr.fault_jump(discr.E.v*W);
            case 'erickson2022'
                for is = 1:nStages
                    faultData.V(:, i_stages(is)) = discr.V_star(T(is), W(:,is));
                end
            end
            faultData.Psi(:, i_stages) = discr.E.Psi*W;
            timeData.T(i_stages) = T;
            timeData.Ht(i_stages) = k*b; % Compute timestep quadrature for current timestep
            timeData.k(i) = k;
            
            % Plot
            if( plotFlag && (mod(i,100) == 0) )
                fh();
                update(discr.getTimeSnapshot(ts));
                drawnow;
                if ~isempty(vidObj)
                    frame = getframe(gcf);
                    writeVideo(vidObj, imresize(frame.cdata, [vidObj.Height, vidObj.Width]));
                end
            end
            
            % Save solution if requested
            if ~isempty(saveOpts) && isfield(saveOpts, 'frameSpacing') && (mod(i, saveOpts.frameSpacing) == 1)
                [w, t] = ts.getV;
                saveData.t = [saveData.t, t];
                saveData.u = [saveData.u, discr.E.u*w];
                saveData.v = [saveData.v, discr.E.v*w];
                saveData.Psi = [saveData.Psi, discr.E.Psi*w];
            end
            % Update iteration counter
            i = i+1; 
        end

        % Prune vectors
        for j = 1:nReceivers
            receiverRecordings{j} = receiverRecordings{j}(1:i_stages(end));
        end
        timeData.k = timeData.k(1:i-1);
        timeData.T = timeData.T(1:i_stages(end));
        timeData.Ht = timeData.Ht(1:i_stages(end));

        faultData.V = faultData.V(1:i_stages(end));
        faultData.Psi = faultData.Psi(1:i_stages(end));
        
        % End progress bar
        if progressBar
            s = util.replace_string(s,'');
        end
        
        if ~isempty(saveOpts) && isfield(saveOpts, 'filename') && ~isempty(saveOpts.filename)
            saveData.grid = obj.forwardDiscr.grid;
            save(saveOpts.filename, 'saveData', '-v7.3','-nocompression');
        end
        
        if plotFlag && ~isempty(vidObj)
            close(vidObj);
        end
    end
    
    % Updates the inversion parameters, pars = pars + steplength*direction;
    % pars 			-- parameter set in native format
    % direction 	-- direction in native gradient format
    % steplength 	-- scalar
    function updateParameters(obj, pars, direction, steplength)
        pars = obj.parametersNativeToVector(pars);
        direction = obj.gradientNativeToVector(direction);
        newPars = pars + steplength*direction;
        newPars = obj.parametersVectorToNative(newPars);
        obj.setParameters(newPars);
    end
    
    
    % Overwrites inversion parameters in obj.pars by those specified in pars
    % pars -- parameters in native format
    function setParameters(obj, pars)
        obj.pars.a = pars.a;
    end
    
    % Converts inversion parameters from native format to one column vector
    % parsNative -- parameters in native fomat
    function parsVec = parametersNativeToVector(obj, parsNative)
        parsVec = parsNative.a;
    end
    
    % Converts vector of inversion parameters to the native format
    % parsVec: column vector of inversion parameters
    function parsNative = parametersVectorToNative(obj, parsVec)
        parsNative.a = parsVec;
    end
    
    % Converts gradient from native format to one column vector
    % gradNative -- gradient in native fomat
    function gradVec = gradientNativeToVector(obj, gradNative)
        gradVec = gradNative;
    end
    
    % Converts gradient vector to native format
    % gradVec -- gradient vector
    function gradNative = gradientVectorToNative(obj, gradVec)
        gradNative = gradVec;
    end
    
    % Computes the difference between to gradients in native format.
    % grad1 -- gradient in native format
    % grad2 -- gradient in native format
    function difference = compareGradients(obj, grad1, grad2)
        grad1 = obj.gradientNativeToVector(grad1);
        grad2 = obj.gradientNativeToVector(grad2);
        difference = abs(grad1 - grad2);
        difference = obj.gradientVectorToNative(difference);
    end
    
    function updateForwardDiscr(obj, pars)
        default_arg('pars',obj.pars)
        % TODO: Based on obj.inversionPars, select which parameters to update

        a = pars.a;
        b = pars.b;
        % Update parameters in forward discr
        obj.forwardDiscr.friction.params.a = a;
        obj.forwardDiscr.friction.params.b = b;

        % Update friction functions
        obj.forwardDiscr.setFaultTraction();
        obj.forwardDiscr.setStateEvolution();
    end
    
    % TODO: Consider splitting into 2 functions dep on if data is interpolated.
    % Update adjoint sources, friction functions
    function updateAdjointDiscr(obj)
        % TODO: Based on obj.inversionPars, select which parameters to update
        discr = obj.adjointDiscr;

        % Update parameter values
        a = obj.pars.a;
        b = obj.pars.b;
        discr.friction.params.a = a;
        discr.friction.params.b = b;
        
        %% 
        % Update adjoint sources with misfit residual
        %

        % Make sure that the adjoint sources are the receivers with the correct data
        data = obj.receiverData;
        approx = obj.forwardReceiverRecordings;
        nReceivers = numel(discr.dirac_deltas);
        interp_data = discr.interpolate_data;
        
        % Solve for negative adjoint velocity potential by reversing sign of forcing
        misfit_residual = cell(nReceivers, obj.dim);
                
        T_fwd = obj.forwardTimeIntegrationData.T;
        nStages = obj.forwardTimeIntegrationData.nStages;

        if interp_data
            % Remove stages from data for interpolation
            t_fwd = obj.removeStagedData(T_fwd,nStages);
            for i = 1:nReceivers
                residual = approx{i} - data{i}(T_fwd); % Residual in receiver measurements including substages
                switch obj.misfitType
                case 'displacement'
                    R_data = obj.removeStagedData(obj.integrateResidual(residual),nStages);
                    R_pp = spline(t_fwd, R_data);
                    R = @(t) ppval(R_pp, t);
                    res_data_i = @(t) -R(obj.T-t);
                case 'velocity'
                    r_data = obj.removeStagedData(residual,nStages);
                    r_pp = spline(t_fwd, r_data);
                    r = @(t) ppval(r_pp, t);
                    res_data_i = @(t) r(obj.T-t); % Reverse in time
                otherwise
                    error('misfit %s not implemented', obj.misfitType);
                end
                misfit_residual{i} = res_data_i;
            end
        else
            for i = 1:nReceivers
                residual = approx{i} - data{i}(T_fwd); % Residual in receiver measurements including substages
                switch obj.misfitType
                case 'displacement'
                    R = obj.integrateResidual(residual);
                    res_data_i = -R; % Minus sign due to how the forcing enters the PDE
                case 'velocity'
                    res_data_i = residual;
                otherwise
                    error('misfit %s not implemented', obj.misfitType);
                end
                misfit_residual{i} = fliplr(res_data_i); % Reverse in time
            end
        end

        %% 
        % Update adjoint friction functions with data from forward solve
        %


        % Forward solve data
        V = obj.forwardFaultVariables.V;
        Psi = obj.forwardFaultVariables.Psi;

        % Compute coefficients needed for adjoint friction functions
        F_V = discr.friction.funs.tau_V(V, Psi, a);
        F_Psi = discr.friction.funs.tau_Psi(V, Psi, a);
        G_V = discr.friction.funs.g_V(V, Psi, a, b);
        G_Psi = discr.friction.funs.g_Psi(V, Psi, a, b);

        if interp_data
            % Update source data and functions
            discr.sources.funs = misfit_residual; % Source data is now a function handle
            discr.setInterpPointSources();
            
            % Update fault data and functions
            % Need to remove stages for the interpolation.
            discr.friction.data.tau_V = obj.removeStagedData(F_V, nStages);
            discr.friction.data.tau_Psi = obj.removeStagedData(F_Psi, nStages);
            discr.friction.data.g_V = obj.removeStagedData(G_V, nStages);
            discr.friction.data.g_Psi = obj.removeStagedData(G_Psi, nStages);
            discr.friction.data.t = obj.T-obj.removeStagedData(T_fwd,nStages); % Reverse in time

            discr.setInterpFaultTraction();
            discr.setInterpStateEvolution();
        else
            % Update source data and functions
            discr.sources.data = misfit_residual;
            discr.setPointSources();

            % Update fault data and functions
            discr.friction.data.tau_V = fliplr(F_V);
            discr.friction.data.g_V = fliplr(G_V);
            discr.friction.data.g_Psi = fliplr(G_Psi);
            discr.friction.data.tau_Psi = fliplr(F_Psi);

            discr.setFaultTraction();
            discr.setStateEvolution();
        end
    end

    function updateSecondOrderForwardDiscr(obj)
        % TODO: Based on obj.inversionPars, select which parameters to update
        discr = obj.secondOrderForwardDiscr;

        % Update parameter values
        a = obj.pars.a;
        b = obj.pars.b;
        discr.friction.params.a = a;
        discr.friction.params.b = b;
        
        %% 
        % Update adjoint friction functions with data from forward solve
        %


        % Forward solve data
        V = obj.forwardFaultVariables.V;
        Psi = obj.forwardFaultVariables.Psi;

        % Compute coefficients needed for adjoint friction functions
        F_V = discr.friction.funs.tau_V(V, Psi, a);
        F_Psi = discr.friction.funs.tau_Psi(V, Psi, a);
        G_V = discr.friction.funs.g_V(V, Psi, a, b);
        G_Psi = discr.friction.funs.g_Psi(V, Psi, a, b);
        
        idx = strcmp(obj.supportedPars, obj.inversionPars{1}); 
        F_p = obj.F_p{idx}(V, Psi, a);
        G_p = obj.G_p{idx}(V, Psi, a, b);

        % For now, we skip interp_data
        % if interp_data
        %     error("Interpalation not yet implemented for hessian vector calculations.")
        % else
        % Update fault data and functions
        discr.friction.data.tau_V = F_V;
        discr.friction.data.g_V = G_V;
        discr.friction.data.g_Psi = G_Psi;
        discr.friction.data.tau_Psi = F_Psi;
        discr.friction.data.tau_p = F_p;
        discr.friction.data.g_p = G_p;

        discr.setFaultTraction();
        discr.setStateEvolution();
    end

    function updateSecondOrderAdjointDiscr(obj)
        % TODO: Based on obj.inversionPars, select which parameters to update
        discr = obj.secondOrderAdjointDiscr;

        % Update parameter values
        a = obj.pars.a;
        b = obj.pars.b;
        discr.friction.params.a = a;
        discr.friction.params.b = b;
        
        %% 
        % Update adjoint sources with misfit residual
        %

        % Make sure that the adjoint sources are the receivers with the correct data
        % TODO: Double check that this is the valid approach,
        % i.e. using same data but changing approx to second order
        data = obj.receiverData;
        % TODO: Figure this out. Code runs with approx = obj.forwardReceiverRecording, but not with current line
        % approx = obj.forwardReceiverRecordings;
        % disp(approx)
        approx = obj.secondOrderForwardReceiverRecordings;
        % disp(approx)
        nReceivers = numel(discr.dirac_deltas);
        interp_data = discr.interpolate_data;
        
        % Solve for negative adjoint velocity potential by reversing sign of forcing
        misfit_residual = cell(nReceivers, obj.dim);
                
        T_fwd = obj.forwardTimeIntegrationData.T;
        nStages = obj.forwardTimeIntegrationData.nStages;

        % if interp_data
        %     % Remove stages from data for interpolation
        %     t_fwd = obj.removeStagedData(T_fwd,nStages);
        %     for i = 1:nReceivers
        %         residual = approx{i} - data{i}(T_fwd); % Residual in receiver measurements including substages
        %         switch obj.misfitType
        %         case 'displacement'
        %             R_data = obj.removeStagedData(obj.integrateResidual(residual),nStages);
        %             R_pp = spline(t_fwd, R_data);
        %             R = @(t) ppval(R_pp, t);
        %             res_data_i = @(t) -R(obj.T-t);
        %         case 'velocity'
        %             r_data = obj.removeStagedData(residual,nStages);
        %             r_pp = spline(t_fwd, r_data);
        %             r = @(t) ppval(r_pp, t);
        %             res_data_i = @(t) r(obj.T-t); % Reverse in time
        %         otherwise
        %             error('misfit %s not implemented', obj.misfitType);
        %         end
        %         misfit_residual{i} = res_data_i;
        %     end
        % else
        for i = 1:nReceivers
            % residual = approx{i} - data{i}(T_fwd); % Residual in receiver measurements including substages
            residual = approx{i}; % Delta Q dagger only depends on Delta Q
            switch obj.misfitType
            case 'displacement'
                R = obj.integrateResidual(residual);
                res_data_i = -R; % Minus sign due to how the forcing enters the PDE
            case 'velocity'
                res_data_i = residual;
            otherwise
                error('misfit %s not implemented', obj.misfitType);
            end
            misfit_residual{i} = fliplr(res_data_i); % Reverse in time
        end
        % end

        %% 
        % Update adjoint friction functions with data from forward solve
        %


        % Forward solve data, these come into the functions
        V = obj.forwardFaultVariables.V;
        Psi = obj.forwardFaultVariables.Psi;
        % Adjoint solve data, these are past "raw"
        V_dagger = obj.adjointFaultVariables.V;
        Psi_dagger = obj.adjointFaultVariables.Psi;
        % Second order forward solve data, these are past "raw"
        delta_V = obj.secondOrderForwardFaultVariables.V;
        delta_Psi = obj.secondOrderForwardFaultVariables.Psi;

        % Compute coefficients needed for adjoint friction functions
        F_V = discr.friction.funs.tau_V(V, Psi, a);
        F_Psi = discr.friction.funs.tau_Psi(V, Psi, a);
        G_V = discr.friction.funs.g_V(V, Psi, a, b);
        G_Psi = discr.friction.funs.g_Psi(V, Psi, a, b);

        F_V_V = discr.friction.funs.tau_V_V(V, Psi, a);
        F_V_Psi = discr.friction.funs.tau_V_Psi(V, Psi, a);
        F_V_a = discr.friction.funs.tau_V_a(V, Psi, a);
        F_Psi_Psi = discr.friction.funs.tau_Psi_Psi(V, Psi, a);
        F_Psi_a = discr.friction.funs.tau_Psi_a(V, Psi, a);
        G_V_Psi = discr.friction.funs.g_V_Psi(V, Psi, a, b);
        G_V_V = discr.friction.funs.g_V_V(V, Psi, a, b);
        G_V_a = discr.friction.funs.g_V_a(V, Psi, a, b);
        G_Psi_Psi = discr.friction.funs.g_Psi_Psi(V, Psi, a, b);
        G_Psi_a = discr.friction.funs.g_Psi_a(V, Psi, a, b);

        % if interp_data
        %     % Update source data and functions
        %     discr.sources.funs = misfit_residual; % Source data is now a function handle
        %     discr.setInterpPointSources();
            
        %     % Update fault data and functions
        %     % Need to remove stages for the interpolation.
        %     discr.friction.data.tau_V = obj.removeStagedData(F_V, nStages);
        %     discr.friction.data.tau_Psi = obj.removeStagedData(F_Psi, nStages);
        %     discr.friction.data.g_V = obj.removeStagedData(G_V, nStages);
        %     discr.friction.data.g_Psi = obj.removeStagedData(G_Psi, nStages);
        %     discr.friction.data.t = obj.T-obj.removeStagedData(T_fwd,nStages); % Reverse in time

        %     discr.setInterpFaultTraction();
        %     discr.setInterpStateEvolution();
        % else
        % Update source data and functions
        discr.sources.data = misfit_residual;
        discr.setPointSources();

        % Update fault data and functions
        discr.friction.data.tau_V = fliplr(F_V);
        discr.friction.data.g_V = fliplr(G_V);
        discr.friction.data.g_Psi = fliplr(G_Psi);
        discr.friction.data.tau_Psi = fliplr(F_Psi);
        discr.friction.data.tau_V_V = fliplr(F_V_V);
        discr.friction.data.tau_V_Psi = fliplr(F_V_Psi);
        discr.friction.data.tau_V_a = fliplr(F_V_a);
        discr.friction.data.tau_Psi_Psi = fliplr(F_Psi_Psi);
        discr.friction.data.tau_Psi_a = fliplr(F_Psi_a);
        discr.friction.data.g_V_Psi = fliplr(G_V_Psi);
        discr.friction.data.g_V_V = fliplr(G_V_V);
        discr.friction.data.g_V_a = fliplr(G_V_a);
        discr.friction.data.g_Psi_Psi = fliplr(G_Psi_Psi);
        discr.friction.data.g_Psi_a = fliplr(G_Psi_a);

        discr.friction.data.V_dagger = V_dagger; % Already flipped (from adjoint)
        discr.friction.data.Psi_dagger = Psi_dagger; % Already flipped (from adjoint)
        discr.friction.data.delta_V = fliplr(delta_V);
        discr.friction.data.delta_Psi = fliplr(delta_Psi);

        discr.setFaultTraction();
        discr.setStateEvolution();
        % end
    end
    
    function R = integrateResidual(obj, r)
        % Utilize the same time stepper settings as used by the scheme
        % in order to integrate the residual.
        order = obj.tsOpts.forwardMethod.order;
        ks = obj.forwardTimeIntegrationData.k;
        t0 = obj.forwardTimeIntegrationData.T(1);
        
        % Solve R_t(t,u) = r(t,u), where r is discrete, i.e. one column for each time 
        % level (including stages)
        S_discr = @(idx_t,u) r(idx_t);
        ts = time.ExplicitRungeKuttaDiscreteData(0, [], S_discr, [], t0 , 0, order);
        N = length(ks);
        [~, ~, ~, nStages] = ts.getTableau();
        R = zeros(1,nStages*N);
        for i = 1:N
            ts.k = ks(i);
            ts.step();
            [~, ~, U] = ts.getV;
            ii = (i-1)*nStages+1 : i*nStages;
            R(ii) = U;
        end
        R = R-R(end);
    end
    
    function grad = computeGradient(obj)
        obj.runForward(false);
        obj.updateAdjointDiscr();
        obj.runAdjoint(false);
        grad = obj.gradientFormula();
    end

    function hessianVector = computeHessianVector(obj)
        plotSolutionsFlag = true; % Determines plotting for all runs

        disp("Run forward...")
        saveOpts = struct('videoFilename', 'assets/forward_run.avi');
        obj.runForward(plotSolutionsFlag, obj.T, saveOpts);
        disp("Run forward complete. Update adjoint discr...")
        obj.updateAdjointDiscr();
        disp("Update adjoint discr complete. Run adjoint...")
        saveOpts.videoFilename = 'assets/adjoint_run.avi';
        obj.runAdjoint(plotSolutionsFlag, obj.T, saveOpts);
        disp("Run adjoint complete. Update second order forward discr...")
        obj.updateSecondOrderForwardDiscr();
        disp("Update second order forward discr complete. Run second order forward...")
        saveOpts.videoFilename = 'assets/second_order_forward_run.avi';
        obj.runSecondOrderForward(plotSolutionsFlag, obj.T, saveOpts);
        disp("Run second order forward complete. Update second order adjoint discr...")
        obj.updateSecondOrderAdjointDiscr();
        disp("Update second order adjoint discr complete. Run second order adjoint...")
        saveOpts.videoFilename = 'assets/second_order_adjoint_run.avi';
        obj.runSecondOrderAdjoint(plotSolutionsFlag, obj.T, saveOpts);
        disp("Run second order adjoint complete. Compute hessian vector...")
        hessianVector = obj.hessianVectorFormula();
        disp("Hessian vector computation complete.")
    end
    
    function grad = gradientFormula(obj)
        % TODO: Construct the appropriate gradient formula
        % based in obj.inversionPars

        % Friction parameters
        a = obj.pars.a;
        b = obj.pars.b;
        % Ensure that parameters across optimization and discretizations are 
        % the same.
        % %assert((a == obj.forwardDiscr.friction.params.a) && (a == obj.adjointDiscr.friction.params.a));
        % %assert((b == obj.forwardDiscr.friction.params.b) && (b == obj.adjointDiscr.friction.params.b));

        % Forward variables
        V = obj.forwardFaultVariables.V;
        Psi = obj.forwardFaultVariables.Psi;
        % F_a = obj.F_p{1}(V, Psi, a, obj.forwardDiscr.friction.params.sigma0, obj.forwardDiscr.friction.params.V0, obj.forwardDiscr.friction.params.tau0);
        % G_a = obj.G_p{1}(V, Psi, a, b, obj.forwardDiscr.friction.params.f0, obj.forwardDiscr.friction.params.V0, obj.forwardDiscr.friction.params.D_c);
        idx = strcmp(obj.supportedPars, obj.inversionPars{1}); 
        F_p = obj.F_p{idx}(V, Psi, a);
        G_p = obj.G_p{idx}(V, Psi, a, b);
        
        if obj.adjointDiscr.interpolate_data
            % Remove stage data 
            nStages = obj.adjointTimeIntegrationData.nStages;
            t_adj = obj.removeStagedData(obj.adjointTimeIntegrationData.T, nStages); 
            V_adj = obj.removeStagedData(obj.adjointFaultVariables.V, nStages);
            Psi_adj = obj.removeStagedData(obj.adjointFaultVariables.Psi, nStages);
            
            % Construct interpolants and evaluate in forward stage times
            V_adj = -ppval(spline(obj.T-t_adj,V_adj), obj.forwardTimeIntegrationData.T);
            Psi_adj = ppval(spline(obj.T-t_adj,Psi_adj), obj.forwardTimeIntegrationData.T);
        else
            % Reverse adjoint variables in time for the integration  
            V_adj = fliplr(obj.adjointFaultVariables.V);
            Psi_adj = fliplr(obj.adjointFaultVariables.Psi);
        end
        grad = -(V_adj.*F_p + Psi_adj.*G_p)*obj.Ht;
    end

    function hessVec = hessianVectorFormula(obj)
        % Friction parameters
        a = obj.pars.a;
        b = obj.pars.b;
        % Ensure that parameters across optimization and discretizations are 
        % the same.
        %assert((a == obj.forwardDiscr.friction.params.a) && (a == obj.adjointDiscr.friction.params.a));
        %assert((b == obj.forwardDiscr.friction.params.b) && (b == obj.adjointDiscr.friction.params.b));

        % Forward variables
        V = obj.forwardFaultVariables.V;
        Psi = obj.forwardFaultVariables.Psi;
        % All else
        V_dagger = obj.adjointFaultVariables.V;
        Psi_dagger = obj.adjointFaultVariables.Psi;
        delta_V_dagger = obj.secondOrderAdjointFaultVariables.V;
        delta_Psi_dagger = obj.secondOrderAdjointFaultVariables.Psi;
        delta_V = obj.secondOrderForwardFaultVariables.V;
        delta_Psi = obj.secondOrderForwardFaultVariables.Psi;
        idx = strcmp(obj.supportedPars, obj.inversionPars{1}); 
        F_p = obj.F_p{idx}(V, Psi, a);
        G_p = obj.G_p{idx}(V, Psi, a, b);
        F_a_a = obj.forwardDiscr.friction.funs.tau_a_a(V, Psi, a);
        G_a_a = obj.forwardDiscr.friction.funs.g_a_a(V, Psi, a, b);
        F_V_a = obj.forwardDiscr.friction.funs.tau_V_a(V, Psi, a);
        F_Psi_a = obj.forwardDiscr.friction.funs.tau_Psi_a(V, Psi, a);
        G_V_a = obj.forwardDiscr.friction.funs.g_V_a(V, Psi, a, b);
        G_Psi_a = obj.forwardDiscr.friction.funs.g_Psi_a(V, Psi, a, b);

        parName = obj.inversionPars{1};
        if ~isfield(obj.eps_pertubations, parName)
            error('Missing epsilon perturbation for parameter %s.', parName);
        end
        eps_p = obj.eps_pertubations.(parName);
        
        % if obj.adjointDiscr.interpolate_data
        %     % Remove stage data 
        %     nStages = obj.adjointTimeIntegrationData.nStages;
        %     t_adj = obj.removeStagedData(obj.adjointTimeIntegrationData.T, nStages); 
        %     V_adj = obj.removeStagedData(obj.adjointFaultVariables.V, nStages);
        %     Psi_adj = obj.removeStagedData(obj.adjointFaultVariables.Psi, nStages);
            
        %     % Construct interpolants and evaluate in forward stage times
        %     V_adj = -ppval(spline(obj.T-t_adj,V_adj), obj.forwardTimeIntegrationData.T);
        %     Psi_adj = ppval(spline(obj.T-t_adj,Psi_adj), obj.forwardTimeIntegrationData.T);
        % else
        % Reverse adjoint variables in time for the integration  
        % V_adj = fliplr(obj.adjointFaultVariables.V);
        % Psi_adj = fliplr(obj.adjointFaultVariables.Psi);
        V_dagger = fliplr(V_dagger);
        Psi_dagger = fliplr(Psi_dagger);
        delta_V_dagger = fliplr(delta_V_dagger);
        delta_Psi_dagger = fliplr(delta_Psi_dagger);

        % end
        T1 = -(V_dagger .* F_a_a .* eps_p) * obj.Ht ./ eps_p;
        T2 = -(Psi_dagger .* G_a_a .* eps_p) * obj.Ht ./ eps_p;
        T3 = -(delta_V_dagger .* F_p) * obj.Ht ./ eps_p;
        T4 = -(delta_Psi_dagger .* G_p) * obj.Ht ./ eps_p;
        T5 = -(V_dagger .* (F_V_a .* delta_V + F_Psi_a .* delta_Psi)) * obj.Ht ./ eps_p;
        T6 = -(Psi_dagger .* (G_V_a .* delta_V + G_Psi_a .* delta_Psi)) * obj.Ht ./ eps_p;
        
        fprintf('T1: %e\nT2: %e\nT3: %e\nT4: %e\nT5: %e\nT6: %e\n', T1, T2, T3, T4, T5, T6);
        
        % NOTE: Divison by eps_a to free eps_a scaling dependence
        hessVec = (T1 + T2 + T3 + T4 + T5 + T6);
    end
    
    function grad = computeGradientFD(obj, deltaG)
        % TODO: Construct the appropriate gradient formula
        % based in obj.inversionPars

        % Compute misfit with current parameter
        pars = obj.pars;
        obj.updateForwardDiscr(pars);
        obj.runForward();
        M0 = obj.computeMisfit();
        
        % Update paramter by a step deltaG
        pars.a = pars.a + deltaG;
        obj.updateForwardDiscr(pars);
        obj.runForward();
        M = obj.computeMisfit();
        dM_dg = (M-M0)/deltaG;
        grad = dM_dg;
    end

    function hessVec = computeHessianVectorFD(obj, deltaG)
        % Compute gradient with current parameter
        pars = obj.pars;
        obj.updateForwardDiscr(pars);
        obj.runForward();
        grad0 = obj.computeGradient();
        
        % Update paramter by a step deltaG
        pars.a = pars.a + deltaG;
        obj.pars.a = pars.a;
        obj.updateForwardDiscr(pars);
        obj.runForward();
        % Have to comment out 'assert' bit from gradient function for this to run
        grad = obj.computeGradient();
        % Reset obj.pars.a
        obj.pars.a = pars.a - deltaG;
        dGrad_dg = (grad-grad0)/deltaG;
        hessVec = dGrad_dg;
    end
    
    function M = computeMisfit(obj)
        data = obj.receiverData;
        approx = obj.forwardReceiverRecordings;
        nReceiver = numel(obj.forwardReceiverRecordings);
        M = 0;
        T = obj.forwardTimeIntegrationData.T;
        for i = 1:nReceiver
            err2 = 1/2*(approx{i} - data{i}(T)).^2;
            errInt = err2*obj.Ht;
            M = M + errInt;
        end
    end
    
    % Useful for finding an approximate line search interval
    function relnorm = gradientNorm(obj, grad)
        relnorm = norm(grad/obj.pars.a, inf);
    end
    
    function [data_handle, model_handle] = plotReceiverSignals(obj, receiver, fontsize)
        default_arg('receiver',1);
        data_handle = obj.receiverData{receiver};
        nStages = obj.forwardTimeIntegrationData.nStages;
        approx = obj.removeStagedData(obj.forwardReceiverRecordings{receiver},nStages);
        t_fwd = obj.removeStagedData(obj.forwardTimeIntegrationData.T,nStages);

        data_handle = plot(t_fwd, data_handle(t_fwd), 'linewidth', 2);
        hold on
        model_handle = plot(t_fwd, approx, 'linewidth', 2);
        xlabel('t', 'fontsize', fontsize)
        xlim([0, obj.T])
        ylabel(sprintf('Rec %d',receiver), 'fontsize', fontsize);
        legend({'Data', 'Model'}, 'fontsize', fontsize, 'Location', 'northwest')
        hold off;
    end
        
    function [upPar, upRec, upMis, figure_handle] = setupPlot(obj, trueParset, initialParset, M, M0, iter)
        
        %-- Plot settings ---
        fontsize = 16;
        scrsz = get(0,'ScreenSize');
        figure_handle = figure('Position',[0.1*scrsz(3) 0.05*scrsz(4) 0.65*scrsz(3) 0.8*scrsz(4)]);
        %--------------------
        parfig = subplot(3,1,1);
        par_true = trueParset.a;
        par = initialParset.a;
        dpar = abs(par-par_true)./par_true;
        par_handle = plot(0, dpar, '-x', 'linewidth', 2, 'markersize', 8);
        xlabel('Iterations', 'fontsize', fontsize)
        ylabel('Relative error', 'fontsize', fontsize)
        xlim([0, iter-1]);

        
        parfig = subplot(3,1,2);
        [~, model_handle] = obj.plotReceiverSignals([],fontsize);
        ylimrec = 2*ylim;
        L = ylimrec(2)-ylimrec(1);
        ylimrec = [ylimrec(1)-0.1*L, ylimrec(2)+0.1*L];
        ylim(ylimrec);
        
        parfig = subplot(3,1,3);
        misfit = semilogy(0:length(M)-1, M/M0, 'linewidth', 2);
        xlabel('Iterations', 'fontsize', fontsize)
        ylabel('Misfit', 'fontsize', fontsize)
        xlim([0, iter-1]);
        ylim([1e-3, 1]);
        
        function updateParameterComparison(object)

        end
        
        function updateReceiver(object, receiver)
            default_arg('receiver', 1);
            approx = object.forwardReceiverRecordings{receiver};
            model_handle.YData = approx;
        end
        
        function updateMisfit(M)
            misfit.YData = M/M0;
        end
        
        upPar = @updateParameterComparison;
        upRec = @updateReceiver;
        upMis = @updateMisfit;
        
    end
    
    function plotGradient(obj, grad)
        warning('Not implemented');
    end

    function obj = setSyntheticReceiverData(obj, data, timeIntegrationData)
        nStages = timeIntegrationData.nStages;
        tRec = obj.removeStagedData(timeIntegrationData.T,nStages);
        for i = length(data)
            recording = obj.removeStagedData(data{i},nStages);
            data_pp = spline(tRec,recording);
            receiverData{i} = @(t) ppval(data_pp,t);
        end
        obj.receiverData = receiverData;
    end
    
end

methods (Static)

    function v = removeStagedData(V,nStages)
        if iscolumn(V)
            V = transpose(V);
        end
        v = [V(1) V(nStages:nStages:end)];
    end
    
    function s = receiverDataStruct()
        s = struct;
        s.m = []; % measured value for all time levels
        s.M = []; % measured value for all time levels including substages 
    end
    
    function s = faultVariablesStruct()
        s = struct;
        s.V = [];
        s.Psi = [];
    end
    
    function s = timeIntegrationStruct()
        s = struct;
        s.nStages = []; % Number of stages
        s.T = []; % Time levels including substages
        s.k = []; % Timesteps
        s.Ht = []; % Quadrature weights
    end
    
    function s = saveStruct()
        s = struct;
        s.grid = [];
        s.t = [];
        s.u = [];
        s.v = [];
        s.Psi = [];
    end
    
end

end