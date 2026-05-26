% Adjoint optimization for anti-plane shear in 1D with a frictional
% fault interface

classdef AntiplaneShear2DRSFrictionOptHessian < adjopt.AdjointOptimization

properties
    % Sources
    adjointReceiverRecordings 	%Recorded during LATEST adjoint simulation (at forward source coords).
    %Stored with FORWARD time convention, like everything else.
    
    % Receivers
    receiverData 				%True data, recorded by "seismometers"
    forwardReceiverRecordings	%Recorded during LATEST forward simulation
    misfitType
    RecMat

    
    % Discretizations
    forwardDiscr
    adjointDiscr
    
    % Time stepping options
    tsOpts
    
    %Data recorded during LATEST forward simulation
    forwardFaultVariables	   
    forwardTimeIntegrationData 
    
    %Data recorded during LATEST adjoint simulation
    adjointFaultVariables	    
    adjointTimeIntegrationData

    % Directional second-order (incremental) data, computed by finite differences
    secondOrderForwardReceiverRecordings
    secondOrderForwardFaultVariables
    secondOrderForwardTimeIntegrationData
    secondOrderAdjointFaultVariables
    secondOrderAdjointTimeIntegrationData
    secondOrderDirection
    secondOrderEpsilon
    
    pars			%Current parameters.
    T 				%Final time
    H               % Space quadrature
    H_b
    Ht              % Time quadrature for LATEST forward simulation
    
    % Misc
    m				%Grid points
    order
    dim
    k
    domain
    
    inversionPars
    F_p
    G_p

    % Interpolation opts for lower dimensional representation
    If2c
    Ic2f

    % Filtering opts
    filterOpts
    filter
    filterTransp
end

methods
    
    function obj = AntiplaneShear2DRSFrictionOptHessian(inputPars, inversionPars, order, k)
        default_arg('inversionPars', {'a'});
        default_arg('order', 4);
        default_arg('k', []);
        
        % Unpack parameter struct
        domain = inputPars.domain;
        opset = inputPars.opset;
        m = inputPars.m;
        material = inputPars.material;
        bc = inputPars.bc;
        friction = inputPars.friction;
        sources = inputPars.sources;
        receivers = inputPars.receivers;
        ic = inputPars.initialconditions;
        tsOpts = inputPars.tsOpts;
        T = inputPars.T;
        misfitType = inputPars.misfitType;
        filterOpts = inputPars.filterOpts;
        
        % Create forward discretization object
        forwardDiscr = elastic.AntiplaneShear2DRSFrictionFwdDiscr(opset, domain, m, order, material, bc, friction, sources, ic);
        % Initialize adjoint discretization object. Does not set any friction/receiver data.
        adjointDiscr = elastic.AntiplaneShear2DRSFrictionAdjDiscr(opset, domain, m, order, material, bc, friction, receivers);

        bidm = domain.boundaryGroups.fault_minus;

        if (isfield(inputPars,'m_p'))
            m_p = inputPars.m_p;
        else
            m_p = m;
        end
        assert(m_p <= m);
        if m_p < m % Fewer grid points used for inversion parameters. Use SBP-preserving interpolation
            order_p = 2;
            opset_p = @sbp.D2Variable;
            coarseGridDiscr = elastic.AntiplaneShear2DRSFrictionFwdDiscr(opset_p, domain, {[m_p, m_p],[m_p, m_p]}, order_p, material, bc, friction, sources, ic);
            pars = coarseGridDiscr.friction.rsParams;
            H_b = coarseGridDiscr.mbDiffOp.getBoundaryQuadrature(bidm);
            glueOps = sbp.GlueOpsOP(m,m_p,order,order_p,opset,opset_p,false);
            If2c = glueOps.Iu2v.bad;
            Ic2f = glueOps.Iv2u.good;
            % TODO: Should add methods for extracting the surface volume factors. 
            % Below is hardcoded.
            blockID = bidm{1}{1};
            boundary_name = sprintf('s_%s',bidm{1}{2});
            Jc = diag(coarseGridDiscr.mbDiffOp.diffOps{blockID}.(boundary_name));
            Jf = diag(forwardDiscr.mbDiffOp.diffOps{blockID}.(boundary_name));
            
            If2c = spdiag(1./sqrt(Jc))*If2c*spdiag(sqrt(Jf));
            Ic2f = spdiag(1./sqrt(Jf))*Ic2f*spdiag(sqrt(Jc));
            
            % Interpolate coarsegrid parameters
            interpParams = obj.interpParameters(pars, Ic2f);
            % Update discrs with interpolated values
            forwardDiscr.friction.rsParams = interpParams;
            forwardDiscr.setFaultTraction();
            forwardDiscr.setStateEvolution();
            adjointDiscr.friction.rsParams = interpParams;
        else % No interpolation needed
            pars = forwardDiscr.friction.rsParams;
            
            H_b = forwardDiscr.mbDiffOp.getBoundaryQuadrature(bidm);

            If2c = speye(m);
            Ic2f = speye(m);
        end
        
        F_p = cell(numel(inversionPars), 1);
        G_p = cell(numel(inversionPars), 1);
        for i = 1:numel(inversionPars)
            parName = inversionPars{i};
            [F_p{i}, G_p{i}] = elastic.friction.inversion.getPartialDerivativeFun(parName);
        end
        
        % Set instance variables
        obj.dim = 2;
        obj.forwardDiscr = forwardDiscr;
        obj.adjointDiscr = adjointDiscr; % Note: Must be updated with data prior to calling runAdjoint
        
        obj.tsOpts = tsOpts;
        
        obj.forwardReceiverRecordings = [];
        obj.adjointReceiverRecordings = [];
        obj.receiverData = [];
        obj.misfitType = misfitType;

        H = forwardDiscr.H;
        switch misfitType
        case 'displacement'
            % Record u
            Erec = forwardDiscr.E.u;
        case 'velocity'
            % Record v
            Erec = forwardDiscr.E.v;
        end
        obj.RecMat = adjointDiscr.dirac_deltas'*H*Erec;

        filt = [];
        filtTransp = [];
        if ~isempty(filterOpts)
            [b,a] = butter(filterOpts.order,filterOpts.cutoffFreq);
            Nfilt = inputPars.T/(0.5*tsOpts.k);
            It = eye(Nfilt);
            if filterOpts.zero_delay            
                filt = filtfilt(b,a,It)';
            else
                filt = filter(b,a,It,[],2);
            end
            if filterOpts.filterResidual
                filtTransp = filt';
            end
        end
        
        obj.forwardFaultVariables = obj.faultVariablesStruct();
        obj.adjointFaultVariables = obj.faultVariablesStruct();
        
        obj.forwardTimeIntegrationData = obj.timeIntegrationStruct();
        obj.adjointTimeIntegrationData = obj.timeIntegrationStruct();
        obj.secondOrderForwardReceiverRecordings = [];
        obj.secondOrderForwardFaultVariables = obj.faultVariablesStruct();
        obj.secondOrderForwardTimeIntegrationData = obj.timeIntegrationStruct();
        obj.secondOrderAdjointFaultVariables = obj.faultVariablesStruct();
        obj.secondOrderAdjointTimeIntegrationData = obj.timeIntegrationStruct();
        obj.secondOrderDirection = struct;
        obj.secondOrderEpsilon = [];

        obj.domain = domain;
        obj.m = m;
        obj.order = order;
        obj.T = T;
        obj.pars = pars;
        obj.inversionPars = inversionPars;
        obj.H = H;
        obj.H_b = H_b;
        obj.If2c = If2c;
        obj.Ic2f = Ic2f; 
        obj.Ht = [];
        obj.k = k;
        obj.F_p = F_p;
        obj.G_p = G_p;
        obj.filterOpts = filterOpts;
        obj.filter = filt;
        obj.filterTransp = filtTransp;
    end
    
    function runForward(obj, plotFlag, T, saveOpts, progressBar)
        default_arg('plotFlag', false);
        default_arg('T', obj.T);
        default_arg('saveOpts', []);
        default_arg('progressBar', []);
        
        discr = obj.forwardDiscr;
        % Matrix for measuring misfit
        Rec = obj.RecMat;
        % Time stepping options
        opts.method = obj.tsOpts.forwardMethod;
        opts.T = T;
        opts.k = obj.tsOpts.k;
        opts.cont_time = true;

        if opts.method.adaptive
            [recData, faultData, timeData] = ...
            obj.runSimulationAdaptive(discr, opts, Rec, plotFlag, saveOpts, progressBar);
        else
            [recData, faultData, timeData] = ...
            obj.runSimulation(discr, opts, Rec, plotFlag, saveOpts, progressBar);
        end
        obj.forwardReceiverRecordings = recData;
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
        % Time stepping options
        opts.method = obj.tsOpts.adjointMethod;
        opts.T = T;
        opts.cont_time = false;
        
        % Time integrate using a standard RK method, but with time steps taken from forward problem
        opts.k = fliplr(obj.forwardTimeIntegrationData.k); % Reverse timestep vector
        [~, faultData, timeData] = obj.runSimulation(discr, opts, [], plotFlag, saveOpts, progressBar);
        obj.adjointFaultVariables = faultData;
        obj.adjointTimeIntegrationData = timeData;
    end

    function runSecondOrderForward(obj, plotFlag, T, saveOpts, progressBar)
        default_arg('plotFlag', false);
        default_arg('T', obj.T);
        default_arg('saveOpts', []);
        default_arg('progressBar', false);

        if isempty(fieldnames(obj.secondOrderDirection)) || isempty(obj.secondOrderEpsilon)
            error('Second-order direction is not set. Call setSecondOrderDirection first.');
        end

        eps2 = obj.secondOrderEpsilon;
        parsBase = obj.pars;
        disp("eps2 = " + num2str(eps2));

        % Base trajectory
        obj.updateForwardDiscr(parsBase);
        obj.runForward(plotFlag, T, saveOpts, progressBar);
        recBase = obj.forwardReceiverRecordings;
        faultBase = obj.forwardFaultVariables;
        timeBase = obj.forwardTimeIntegrationData;
        HtBase = obj.Ht;

        % Perturbed trajectory p + eps * dp
        parsPlus = obj.addDirectionToPars(parsBase, obj.secondOrderDirection, eps2);
        obj.setParameters(parsPlus);
        obj.updateForwardDiscr(parsPlus);
        obj.runForward(plotFlag, T, saveOpts, progressBar);
        recPlus = obj.forwardReceiverRecordings;
        faultPlus = obj.forwardFaultVariables;

        % Directional derivative approximation
        obj.secondOrderForwardReceiverRecordings = (recPlus - recBase)/eps2;
        obj.secondOrderForwardFaultVariables.V = (faultPlus.V - faultBase.V)/eps2;
        obj.secondOrderForwardFaultVariables.Psi = (faultPlus.Psi - faultBase.Psi)/eps2;
        obj.secondOrderForwardTimeIntegrationData = timeBase;

        % Restore base state
        obj.setParameters(parsBase);
        obj.updateForwardDiscr(parsBase);
        obj.forwardReceiverRecordings = recBase;
        obj.forwardFaultVariables = faultBase;
        obj.forwardTimeIntegrationData = timeBase;
        obj.Ht = HtBase;
    end

    function runSecondOrderAdjoint(obj, plotFlag, T, saveOpts, progressBar)
        default_arg('plotFlag', false);
        default_arg('T', obj.T);
        default_arg('saveOpts', []);
        default_arg('progressBar', false);

        if isempty(fieldnames(obj.secondOrderDirection)) || isempty(obj.secondOrderEpsilon)
            error('Second-order direction is not set. Call setSecondOrderDirection first.');
        end

        eps2 = obj.secondOrderEpsilon;
        parsBase = obj.pars;

        % Base forward+adjoint trajectory
        obj.updateForwardDiscr(parsBase);
        obj.runForward(plotFlag, T, saveOpts, progressBar);
        recBase = obj.forwardReceiverRecordings;
        faultFwdBase = obj.forwardFaultVariables;
        timeFwdBase = obj.forwardTimeIntegrationData;
        HtBase = obj.Ht;
        obj.updateAdjointDiscr();
        obj.runAdjoint(plotFlag, T, saveOpts, progressBar);
        faultAdjBase = obj.adjointFaultVariables;
        timeAdjBase = obj.adjointTimeIntegrationData;

        % Perturbed forward+adjoint trajectory p + eps * dp
        parsPlus = obj.addDirectionToPars(parsBase, obj.secondOrderDirection, eps2);
        obj.setParameters(parsPlus);
        obj.updateForwardDiscr(parsPlus);
        obj.runForward(plotFlag, T, saveOpts, progressBar);
        obj.updateAdjointDiscr();
        obj.runAdjoint(plotFlag, T, saveOpts, progressBar);
        faultAdjPlus = obj.adjointFaultVariables;

        % Directional derivative approximation
        obj.secondOrderAdjointFaultVariables.V = (faultAdjPlus.V - faultAdjBase.V)/eps2;
        obj.secondOrderAdjointFaultVariables.Psi = (faultAdjPlus.Psi - faultAdjBase.Psi)/eps2;
        obj.secondOrderAdjointTimeIntegrationData = timeAdjBase;

        % Restore base state and consistent discretizations
        obj.setParameters(parsBase);
        obj.updateForwardDiscr(parsBase);
        obj.forwardReceiverRecordings = recBase;
        obj.forwardFaultVariables = faultFwdBase;
        obj.forwardTimeIntegrationData = timeFwdBase;
        obj.Ht = HtBase;
        obj.updateAdjointDiscr();
        obj.adjointFaultVariables = faultAdjBase;
        obj.adjointTimeIntegrationData = timeAdjBase;
    end

    function setSecondOrderDirection(obj, direction, epsilon)
        default_arg('epsilon', 1e-6);
        obj.secondOrderDirection = obj.normalizeDirectionStruct(direction);
        obj.secondOrderEpsilon = epsilon;
    end

    function hv = computeHessianVector(obj, direction, epsilon, plotFlag, progressBar)
        default_arg('epsilon', 1e-6);
        default_arg('plotFlag', false);
        default_arg('progressBar', false);

        directionNative = obj.normalizeDirectionStruct(direction);
        obj.setSecondOrderDirection(directionNative, epsilon);

        parsBase = obj.pars;

        % Base gradient g(p)
        gBase = obj.computeGradient(plotFlag, progressBar);

        % Perturbed gradient g(p + eps*dp)
        parsPlus = obj.addDirectionToPars(parsBase, directionNative, epsilon);
        obj.setParameters(parsPlus);
        obj.updateForwardDiscr(parsPlus);
        gPlus = obj.computeGradient(plotFlag, progressBar);

        % Directional derivative of the gradient
        hv = struct;
        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            hv.(parName) = (gPlus.(parName) - gBase.(parName))/epsilon;
        end

        % Cache second-order trajectories for diagnostics
        obj.setParameters(parsBase);
        obj.updateForwardDiscr(parsBase);
        obj.runSecondOrderForward(false, obj.T, [], progressBar);
        obj.runSecondOrderAdjoint(false, obj.T, [], progressBar);

        % Restore optimization state at base point
        obj.setParameters(parsBase);
        obj.updateForwardDiscr(parsBase);
    end

    function [epsList, relErr] = validateHessianVectorFD(obj, direction, epsList, progressBar)
        default_arg('epsList', 10.^(-3:-1:-7));
        default_arg('progressBar', false);

        directionNative = obj.normalizeDirectionStruct(direction);
        hvRef = obj.computeHessianVector(directionNative, min(epsList)/10, false, progressBar);
        hvRefVec = obj.gradientNativeToVector(hvRef);
        relErr = zeros(size(epsList));

        for i = 1:numel(epsList)
            hvFD = obj.computeHessianVector(directionNative, epsList(i), false, progressBar);
            hvFDVec = obj.gradientNativeToVector(hvFD);
            relErr(i) = norm(hvFDVec - hvRefVec)/(norm(hvRefVec) + eps);
        end
    end

    % =====================================================================
    % ANALYTIC SECOND-ORDER METHODS (BOILERPLATE)
    % =====================================================================

    function hv = computeHessianVectorAnalytic(obj, direction)
        directionNative = obj.normalizeDirectionStruct(direction);
        
        % 1. Run Base Trajectories (Forward & Adjoint)
        obj.runForward(false, obj.T, [], false);
        obj.updateAdjointDiscr();
        obj.runAdjoint(false, obj.T, [], false);
        
        pars = obj.forwardDiscr.friction.rsParams;
        V = obj.forwardFaultVariables.V;
        Psi = obj.forwardFaultVariables.Psi;
        V_adj = fliplr(obj.adjointFaultVariables.V);
        Psi_adj = fliplr(obj.adjointFaultVariables.Psi);
        
        % 2. Solve Linearized Forward Problem for Delta V and Delta Psi
        [dV, dPsi] = obj.runLinearizedForward(directionNative);
        
        % 3. Assemble the analytic Hessian-Vector Product
        hv = struct;
        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            eps_p = directionNative.(parName); % perturbation delta p
            
            % TODO: Define these function handles in your friction properties
            % F_pp = obj.F_pp{i}(V, Psi, pars...); 
            % G_pp = obj.G_pp{i}(V, Psi, pars...);
            % F_pV = obj.F_pV{i}(V, Psi, pars...);
            % F_pPsi = obj.F_pPsi{i}(V, Psi, pars...);
            % G_pV = obj.G_pV{i}(V, Psi, pars...);
            % G_pPsi = obj.G_pPsi{i}(V, Psi, pars...);
            % F_p = obj.F_p{i}(V, Psi, pars...);
            % G_p = obj.G_p{i}(V, Psi, pars...);
            
            % Placeholder arrays (zeros) to represent evaluated dependencies
            F_pp = zeros(size(V)); G_pp = zeros(size(V));
            F_pV = zeros(size(V)); F_pPsi = zeros(size(V));
            G_pV = zeros(size(V)); G_pPsi = zeros(size(V));
            F_p = zeros(size(V));  G_p = zeros(size(V));
            
            % Linearized adjoint variables (if you fully implemented runLinearizedAdjoint)
            % dV_adj_rev = fliplr(dV_adj);
            % dPsi_adj_rev = fliplr(dPsi_adj);
            dV_adj_rev = zeros(size(V));    % placeholder
            dPsi_adj_rev = zeros(size(V));  % placeholder
            
            % Term 1: Second derivative w.r.t parameter squared
            term1 = V_adj .* F_pp .* eps_p + Psi_adj .* G_pp .* eps_p;
            
            % Term 2 & 3: Coupling with Adjoint Variations
            term2 = dV_adj_rev .* F_p + dPsi_adj_rev .* G_p;
            
            % Term 4 & 5: Mixed second derivatives with State Variations
            term3 = V_adj .* (F_pV .* dV + F_pPsi .* dPsi);
            term4 = Psi_adj .* (G_pV .* dV + G_pPsi .* dPsi);
            
            % Assemble the exact spatial-temporal integral
            total_integrand = term1 + term2 + term3 + term4;
            
            hv.(parName) = -(obj.H_b * obj.If2c * total_integrand) * obj.Ht;
        end
    end

    function [dV, dPsi] = runLinearizedForward(obj, directionNative)
        % BOILERPLATE: Integrate the Second-Order Forward Problem
        % Delta F_j = (dF/dV)*dV + (dF/dPsi)*dPsi + (dF/dp)*eps_p
        % Delta G_j = (dG/dV)*dV + (dG/dPsi)*dPsi + (dG/dp)*eps_p
        
        % In practice, you must construct a specialized timestepper 
        % that applies Delta F as the boundary traction and Delta G as the 
        % state evolution, utilizing the base trajectory matrices.
        
        [nf, nSteps] = size(obj.forwardFaultVariables.V);
        dV = zeros(nf, nSteps);
        dPsi = zeros(nf, nSteps);
        
        % TODO: Provide linearized time-stepping loop here
    end

    function [dV_adj, dPsi_adj] = runLinearizedAdjoint(obj, directionNative, dV, dPsi)
        % BOILERPLATE: Integrate the Second-Order Adjoint Problem
        % Driven by Delta Q_adj = sum Delta \dot{u} delta(x - x_r)
        
        [nf, nSteps] = size(obj.adjointFaultVariables.V);
        dV_adj = zeros(nf, nSteps);
        dPsi_adj = zeros(nf, nSteps);
        
        % TODO: Provide linearized adjoint time-stepping loop here
    end
    
    
    function [receiverRecordings, faultData, timeData] = runSimulation(obj, discr, tsOpts,...
        Rec, plotFlag, saveOpts, progressBar)
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
        
        [nf, ~] = size(discr.E.usim);
        
        % Initialize data structs
        if ~isempty(Rec)
            [nReceivers,~] = size(Rec);
            receiverRecordings = zeros(nReceivers, N*nStages);
        else
            receiverRecordings = [];
        end
        
        faultData = obj.faultVariablesStruct();
        faultData.V = zeros(nf, N*nStages); 
        faultData.Psi = zeros(nf, N*nStages); 
        
        
        timeData = obj.timeIntegrationStruct();
        timeData.nStages = nStages;
        timeData.T = zeros(1, N*nStages);
        timeData.Ht = zeros(1, N*nStages);
        timeData.k = ks;
        
        % Setup plot
        if plotFlag
            plotOpts.plot_variables = 'all';
            [update, fh] = discr.setupPlot(plotOpts);
            fh();
            r = discr.getTimeSnapshot(0);
            update(r);
        end

        % Save solution if requested
        % Initial save. Include grid and fault in addition to solution fields
        if ~isempty(saveOpts)
            saveData.grid = discr.grid;
            saveData.fault = discr.grid.getBoundary(obj.domain.boundaryGroups.fault_minus);

            [w, t] = ts.getV;
            saveData.t = t;
            saveData.u = discr.E.u*w;
            saveData.v = discr.E.v*w;
            saveData.Psi = discr.E.Psi*w;
            saveData.D = discr.fault_jump(saveData.u);
            if tsOpts.cont_time
                saveData.V = discr.V_star(t, w);
                saveData.tau = discr.penalty_fault\discr.fault_traction(t, w, saveData.V);
            else
                saveData.V = discr.V_star(1, w);
                saveData.tau = discr.penalty_fault\discr.fault_traction(1, w, saveData.V);
            end

            save_counter = 1;
            filename = sprintf('%s/%d.mat',saveOpts.saveDir, save_counter);
            save(filename, 'saveData', '-v7.3','-nocompression');
        end
        
        % Initialize progress bar
        if progressBar
            s = util.replace_string('','   %d %%',0);
        end
        
        for i = 1:N;
            ts.k = ks(i); % Set timestep and step forward.
            ts.step();
            
            [~, ~, W, T] = ts.getV;
            i_stages = (i-1)*nStages+1 : i*nStages;
            
            % Update progress bars
            if progressBar
                s = util.replace_string(s,'   %.2f %%',i/N*100);
            end
            
            % Store at receivers
            if ~isempty(Rec)	
                receiverRecordings(:,i_stages) = Rec*W;
            end
            
            % Store fault data and time integration data for the iteration
            if tsOpts.cont_time
                for is = 1:nStages
                    faultData.V(:, i_stages(is)) = discr.V_star(T(is), W(:,is));
                    faultData.Psi(:, i_stages(is)) = discr.E.Psi*W(:,is);
                end
            else
                for is = 1:nStages
                    faultData.V(:, i_stages(is)) = discr.V_star(i_stages(is), W(:,is));
                    faultData.Psi(:, i_stages(is)) = discr.E.Psi*W(:,is);
                end
            end
            %faultData.V(:, i_stages) = discr.E.Vs*W;
            
            timeData.T(i_stages) = T;
            timeData.Ht(i_stages) = ts.getTimeStepQuadrature();
            
            % Plot
            if( plotFlag && (mod(i,100) == 0) )
                fh();
                update(discr.getTimeSnapshot(ts));
                drawnow;
            end
            
            % Save solution if requested
            if ~isempty(saveOpts) && (mod(i, saveOpts.frameSpacing) == 1)
                [w, t] = ts.getV;
                saveData.t = t;
                saveData.u = discr.E.u*w;
                saveData.v = discr.E.v*w;
                saveData.Psi = faultData.Psi(:,i_stages(end));
                saveData.D = discr.fault_jump(saveData.u);
                saveData.V = faultData.V(:,i_stages(end));
                if tsOpts.cont_time
                    saveData.tau = discr.penalty_fault\discr.fault_traction(t, w, saveData.V);
                else
                    saveData.tau = discr.penalty_fault\discr.fault_traction(i_stages(end), w, saveData.V);
                end

                save_counter = save_counter+1;
                filename = sprintf('%s/%d.mat',saveOpts.saveDir, save_counter);
                save(filename, 'saveData', '-v7.3','-nocompression');
            end
        end
        % Plot at final time
        if plotFlag
            fh();
            update(discr.getTimeSnapshot(ts));
            drawnow;
        end
        
        % End progress bar
        if progressBar
            s = util.replace_string(s,'');
        end
        
        % Save last time level if requested
        if ~isempty(saveOpts)
            [w, t] = ts.getV;
            saveData.t = t;
            saveData.u = discr.E.u*w;
            saveData.v = discr.E.v*w;
            saveData.Psi = faultData.Psi(:,i_stages(end));
            saveData.D = discr.fault_jump(saveData.u);
            saveData.V = faultData.V(:,i_stages(end));
            if tsOpts.cont_time
                saveData.tau = discr.penalty_fault\discr.fault_traction(t, w, saveData.V);
            else
                saveData.tau = discr.penalty_fault\discr.fault_traction(i_stages(end), w, saveData.V);
            end
            
            save_counter = save_counter+1;
            filename = sprintf('%s/%d.mat',saveOpts.saveDir, save_counter);
            save(filename, 'saveData', '-v7.3','-nocompression');
        end
    end
    
    function [receiverRecordings, faultData, timeData] = runSimulationAdaptive(obj, discr, tsOpts,...
        Rec, plotFlag, saveOpts, progressBar)
        default_arg('saveOpts', []);
        default_arg('progressBar', false)
        saveData = obj.saveStruct();
        
        [ts, N] = discr.getTimestepper(tsOpts.method, tsOpts.T, tsOpts.k);
        [~, b, ~, nStages] = ts.getTableau();
        
        [nf, ~] = size(discr.E.usim);
        
        % Initialize data structs
        buffsz = 1e5; % Initial guess of number of timelevels.
        if ~isempty(Rec)
            [nReceivers,~] = size(Rec);
            receiverRecordings = zeros(nReceivers, buffsz*nStages);
        else
            receiverRecordings = [];
        end

        
        faultData = obj.faultVariablesStruct();
        faultData.V = zeros(nf, buffsz*nStages); 
        faultData.Psi = zeros(nf, buffsz*nStages); 
        
        timeData = obj.timeIntegrationStruct();
        timeData.T = zeros(1, buffsz*nStages);
        timeData.Ht = zeros(1, buffsz*nStages);
        timeData.nStages = nStages;

        % Setup plot
        if plotFlag
            plotOpts.plot_variables = 'all';
            [update, fh] = discr.setupPlot(plotOpts);
            fh();
            r = discr.getTimeSnapshot(0);
            update(r);
        end
        
        % Initialize progress bar
        if progressBar
            s = util.replace_string('','   %d %%',0);
        end

        % Save solution if requested
        % Initial save. Include grid and fault in addition to solution fields
        if ~isempty(saveOpts)
            saveData.grid = discr.grid;
            saveData.fault = discr.grid.getBoundary(obj.domain.boundaryGroups.fault_minus);

            [w, t] = ts.getV;
            saveData.t = t;
            saveData.u = discr.E.u*w;
            saveData.v = discr.E.v*w;
            saveData.Psi = discr.E.Psi*w;
            saveData.D = discr.fault_jump(saveData.u);
            saveData.V = discr.V_star(t, w);

            save_counter = 1;
            filename = sprintf('%s/%d.mat',saveOpts.saveDir, save_counter);
            save(filename, 'saveData', '-v7.3','-nocompression');
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
                if ~isempty(Rec)
                    receiverRecordings = [receiverRecordings, zeros(nReceivers, buffsz*nStages)];
                end
                 
                timeData.k = [timeData.k, zeros(1, buffsz)];
                timeData.T = [timeData.T, zeros(1, buffsz*nStages)];
                timeData.Ht = [timeData.Ht, zeros(1, buffsz*nStages)];

                faultData.V = [faultData.V, zeros(nf, buffsz*nStages)];
                faultData.Psi = [faultData.Psi, zeros(nf, buffsz*nStages)];
                
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
            if ~isempty(Rec)
                receiverRecordings(:,i_stages) = R*W;
            end
            
            % Store fault data and time integration data for the iteration
            k = T(end) - T(1); % Extract the timestep used for the stages;
            for is = 1:nStages
                faultData.V(:, i_stages(is)) = discr.V_star(T(is), W(:,is));
                faultData.Psi(:, i_stages(is)) = discr.E.Psi*W(:,is);
            end
            
            timeData.T(i_stages) = T;
            timeData.Ht(i_stages) = k*b; % Compute timestep quadrature for current timestep
            timeData.k(i) = k;
            
            % Plot
            if( plotFlag && (mod(i,100) == 0) )
                fh();
                update(discr.getTimeSnapshot(ts));
                drawnow;
            end
            
            % Save solution if requested
            if ~isempty(saveOpts) && (mod(i, saveOpts.frameSpacing) == 1)
                [w, t] = ts.getV;
                saveData.t = t;
                saveData.u = discr.E.u*w;
                saveData.v = discr.E.v*w;
                saveData.Psi = faultData.Psi(:,i_stages(end));
                saveData.D = discr.fault_jump(saveData.u);
                saveData.V = faultData.V(:,i_stages(end));
                saveData.V = faultData.V(:,i_stages(end));

                save_counter = save_counter+1;
                filename = sprintf('%s/%d.mat',saveOpts.saveDir, save_counter);
                save(filename, 'saveData', '-v7.3','-nocompression');
            end
            % Update iteration counter
            i = i+1; 
        end
        % Plot at final time
        if plotFlag
            fh();
            update(discr.getTimeSnapshot(ts));
            drawnow;
        end

        % Prune vectors
        if ~isempty(Rec)
            receiverRecordings = receiverRecordings(:,1:i_stages(end));
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
        
        % Save last time level if requested
        if ~isempty(saveOpts)
            [w, t] = ts.getV;
            saveData.t = t;
            saveData.u = discr.E.u*w;
            saveData.v = discr.E.v*w;
            saveData.Psi = faultData.Psi(:,i_stages(end));
            saveData.D = discr.fault_jump(saveData.u);
            saveData.V = faultData.V(:,i_stages(end));
            
            save_counter = save_counter+1;
            filename = sprintf('%s/%d.mat',saveOpts.saveDir, save_counter);
            save(filename, 'saveData', '-v7.3','-nocompression');
        end
    end

    function checkParameters(obj)
        pars_fine = obj.interpParameters(obj.pars, obj.Ic2f);
        fwd_pars = obj.forwardDiscr.friction.rsParams;
        adj_pars = obj.adjointDiscr.friction.rsParams;
        
        names = fieldnames(pars_fine);
        for i = 1:numel(names)
            name = names{i};
            p = pars_fine.(name);
            p_fwd = fwd_pars.(name);
            p_adj = adj_pars.(name);
            assert(all(p == p_fwd) && all(p == p_adj));
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
        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            obj.pars.(parName) = pars.(parName);
        end
    end
    
    % Converts inversion parameters from native format to one column vector
    % parsNative -- parameters in native fomat
    function parsVec = parametersNativeToVector(obj, parsNative)
        pars = cell(numel(obj.inversionPars), 1);
        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            pars{i} = parsNative.(parName);
        end
        parsVec = cell2mat(pars);
    end
    
    % Converts vector of inversion parameters to the native format
    % parsVec: column vector of inversion parameters
    function parsNative = parametersVectorToNative(obj, parsVec)
        parsNative = struct;
        index = 1;
        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            m = length(obj.pars.(parName));
            parsNative.(parName) = parsVec(index : index+m-1);
            index = index + m;
        end
    end
    
    % Converts gradient from native format to one column vector
    % gradNative -- gradient in native fomat
    function gradVec = gradientNativeToVector(obj, gradNative)
        gradVec = obj.parametersNativeToVector(gradNative);
    end
    
    % Converts gradient vector to native format
    % gradVec -- gradient vector
    function gradNative = gradientVectorToNative(obj, gradVec)
        gradNative = obj.parametersVectorToNative(gradVec);
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
        pars_fine = obj.interpParameters(pars, obj.Ic2f);
        obj.forwardDiscr.friction.rsParams = pars_fine;
        % Update friction functions
        obj.forwardDiscr.setFaultTraction();
        obj.forwardDiscr.setStateEvolution();
    end
    
    % TODO: Consider splitting into 2 functions dep on if data is interpolated.
    function updateAdjointDiscr(obj)
        discr = obj.adjointDiscr;
        pars = obj.forwardDiscr.friction.rsParams;

        % Update parameter values
        
        discr.friction.rsParams = pars;
        %% 
        % Update adjoint sources with misfit residual
        %

        % Make sure that the adjoint sources are the receivers with the correct data
        data = obj.receiverData;
        approx = obj.forwardReceiverRecordings;
        [nReceivers, ~] = size(approx);
        
        % Solve for negative adjoint velocity potential by reversing sign of forcing
        misfit_residual = cell(nReceivers, 1);       
        T_fwd = obj.forwardTimeIntegrationData.T;
        filteredSignal = zeros(2,numel(T_fwd)/2);
        % Compute receiver misfit from forward solve
        for i = 1:nReceivers
            approx_i = approx(i,:);
            data_i = data{i}(T_fwd);

            % Perform filtering if applicable
            if ~isempty(obj.filterOpts)
                if obj.filterOpts.filterResidual
                    residual = approx_i - data_i;
                    % First filtering
                    filteredSignal(1,:) = (obj.filter*(residual(1:2:end-1)'))';
                    filteredSignal(2,:) = (obj.filter*(residual(2:2:end)'))';
                    % Adjoint filtering
                    residual = obj.Ht'.*filteredSignal(:)';
                    filteredSignal(1,:) = (obj.filterTransp*(residual(1:2:end-1)'))';
                    filteredSignal(2,:) = (obj.filterTransp*(residual(2:2:end)'))';
                    residual = (1./obj.Ht').*filteredSignal(:)';
                else % Filter data only
                    filteredSignal(1,:) = (obj.filter*(data_i(1:2:end-1)'))';
                    filteredSignal(2,:) = (obj.filter*(data_i(2:2:end)'))';
                    residual = approx_i - filteredSignal(:)';
                end
            else
                residual = approx_i - data_i;
            end
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
        %% 
        % Update adjoint friction functions with data from forward solve
        %
        % Forward solve data
        V = obj.forwardFaultVariables.V;
        Psi = obj.forwardFaultVariables.Psi;

        % Compute coefficients needed for adjoint friction functions
        
        F_V = discr.friction.funs.F_V(V, Psi, pars.a, pars.sigma0, pars.V0, pars.tau0);
        F_Psi = discr.friction.funs.F_Psi(V, Psi, pars.a, pars.sigma0, pars.V0, pars.tau0);
        G_V = discr.friction.funs.G_V(V, Psi, pars.a, pars.b, pars.f0, pars.V0, pars.D_c);
        G_Psi = discr.friction.funs.G_Psi(V, Psi, pars.a, pars.b, pars.f0, pars.V0, pars.D_c);

        % Update source data and functions
        discr.sources.data = misfit_residual;
        discr.setPointSources();

        % Update fault data and functions
        discr.friction.data.F_V = fliplr(F_V);
        discr.friction.data.G_V = fliplr(G_V);
        discr.friction.data.G_Psi = fliplr(G_Psi);
        discr.friction.data.F_Psi = fliplr(F_Psi);

        discr.setFaultTraction();
        discr.setStateEvolution();
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
    
    function grad = computeGradient(obj, plotFlag, progressBar)
        default_arg('plotFlag', false)
        default_arg('progressBar', false)
        obj.runForward(plotFlag, obj.T, [], progressBar);
        obj.updateAdjointDiscr();
        obj.runAdjoint(plotFlag, obj.T, [], progressBar);
        grad = obj.gradientFormula();
    end
    
    function grad = gradientFormula(obj)
        obj.checkParameters();
        % Friction parameters
        pars = obj.forwardDiscr.friction.rsParams;

        % Forward variables
        V = obj.forwardFaultVariables.V;
        Psi = obj.forwardFaultVariables.Psi;

        % Reverse adjoint variables in time for the integration  
        V_adj = fliplr(obj.adjointFaultVariables.V);
        Psi_adj = fliplr(obj.adjointFaultVariables.Psi);
        
        grad = struct;
        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            F_p = obj.F_p{i}(V, Psi, pars.a, pars.sigma0, pars.V0, pars.tau0);
            G_p = obj.G_p{i}(V, Psi, pars.a, pars.b, pars.f0, pars.V0, pars.D_c);
            grad.(parName) = -(obj.H_b*obj.If2c*(V_adj.*F_p + Psi_adj.*G_p))*obj.Ht;
        end
    end
    
    function grad = computeGradientFD(obj, deltaG)
        % TODO: Construct the appropriate gradient formula
        % based in obj.inversionPars

        % Compute misfit with current parameter vector
        pars = obj.pars;
        obj.updateForwardDiscr(pars);
        obj.runForward();
        M0 = obj.computeMisfit();
        
        % Update each element in paramter vector by a step deltaG
        for j = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{j};
            grad.(parName) = zeros(size(pars.(parName)));
            for i = 1:length(pars.(parName))
                pars_tmp = pars;
                pars_tmp.(parName)(i) = pars_tmp.(parName)(i) + deltaG;
                %pars = pars_tmp; %
                obj.updateForwardDiscr(pars_tmp);
                obj.runForward();
                M = obj.computeMisfit();
                dM_dg = (M-M0)/deltaG;
                grad.(parName)(i) = dM_dg;
            end
        end
    end
    
    function M = computeMisfit(obj)
        data = obj.receiverData;
        approx = obj.forwardReceiverRecordings;
        [nReceiver,~] = size(approx);
        M = 0;
        T = obj.forwardTimeIntegrationData.T;
        filteredSignal = zeros(2,numel(T)/2);
        for i = 1:nReceiver
            approx_i = approx(i,:);
            data_i = data{i}(T);
            if ~isempty(obj.filterOpts)
                if obj.filterOpts.filterResidual
                    residual = approx_i - data_i;
                    filteredSignal(1,:) = (obj.filter*(residual(1:2:end-1)'))';
                    filteredSignal(2,:) = (obj.filter*(residual(2:2:end)'))';
                    residual = filteredSignal(:)';        
                else % Only filter data
                    filteredSignal(1,:) = (obj.filter*(data_i(1:2:end-1)'))';
                    filteredSignal(2,:) = (obj.filter*(data_i(2:2:end)'))';
                    residual = approx_i - filteredSignal(:)';
                end
            else 
                residual = approx_i - data_i;
            end
            err2 = residual.^2;
            errInt = 1/2*err2*obj.Ht;
            M = M + errInt;
        end
    end
    
    % Useful for finding an approximate line search interval
    function relnorm = gradientNorm(obj, grad)
        relnorm = zeros(1, numel(obj.inversionPars) );
        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            p = obj.pars.(parName);
            grad_p = grad.(parName);
            relnorm(i) = norm(grad_p./p, inf);
        end
        relnorm = max(relnorm);
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
            approx = object.forwardReceiverRecordings(receiver,:);
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

    function directionNative = normalizeDirectionStruct(obj, direction)
        if isnumeric(direction)
            directionNative = obj.parametersVectorToNative(direction(:));
        elseif isstruct(direction)
            directionNative = direction;
        else
            error('Direction must be a vector or a struct in native parameter format.');
        end

        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            if ~isfield(directionNative, parName)
                error('Direction is missing inversion parameter %s.', parName);
            end
            p = obj.pars.(parName);
            dp = directionNative.(parName);
            if ~isequal(size(p), size(dp))
                error('Direction size mismatch for parameter %s.', parName);
            end
        end
    end

    function parsOut = addDirectionToPars(obj, parsIn, direction, alpha)
        parsOut = parsIn;
        for i = 1:numel(obj.inversionPars)
            parName = obj.inversionPars{i};
            parsOut.(parName) = parsIn.(parName) + alpha*direction.(parName);
        end
    end

    function obj = setSyntheticReceiverData(obj, data, timeIntegrationData)
        nStages = timeIntegrationData.nStages;
        tRec = obj.removeStagedData(timeIntegrationData.T,nStages);
        [nRec, ~] = size(data);
        receiverData = cell(nRec,1);
        for i = 1:nRec
            recording = obj.removeStagedData(data(i,:),nStages);
            data_pp = spline(tRec,recording);
            receiverData{i} = @(t) ppval(data_pp,t);
        end
        obj.receiverData = receiverData;
    end
    
end

methods (Static)

    function pars_fine = interpParameters(pars, Ic2f)
        pars_fine = struct;
        names = fieldnames(pars);
        for i = 1:numel(names)
            name = names{i};
            pars_fine.(name) = Ic2f*pars.(name);
        end
    end

    function v = removeStagedData(V,nStages)
        if isvector(V)
            if iscolumn(V)
                V = transpose(V);
            end
            v = [V(1) V(nStages:nStages:end)];
        else
            assert(ismatrix(V));
            v = [V(:,1) V(:,nStages:nStages:end)];
        end
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
    end
    
end

end