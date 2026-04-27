% Adjoint optimization for 2D elastic, Cartesian

classdef Elastic2DMultiblockOpt < adjopt.AdjointOptimization

	properties

		% Sources
		sourcePoints				%Coordinates
		sourceDeltas				%Delta functions
		adjointReceiverRecordings 	%Recorded during LATEST adjoint simulation (at forward source coords).
									%Stored with FORWARD time convention, like everything else.

		% Receivers
		receiverPoints				%Coordinates
		receiverData 				%True data, recorded by "seismometers"
		receiverDeltas				%Delta functions
		forwardReceiverRecordings	%Recorded during LATEST forward simulation

		forwardDiscr
		adjointDiscr
		pars			%Current medium parameters.

		forwardField	%Recorded during LATEST forward simulation
		adjointField	%Recorded during LATEST adjoint simulation

		T 				%Final time
		k 				%Time step
		tvec 			%Includes substage points

		m
		order
		H           	%Total space quadrature, for both components
		Ht 				%Time quadrature, a column vector of weights.
		dim

		inversionPars

		% Tensors and matrices for inverting for lambda, mu
		B
		D1
		optFlag

		% Dummy diffOps with (lambda = 1, mu = 0) and (lambda = 0, mu = 1)
		diffOp_mu_1
		diffOp_lambda_1

	end

	methods

		function obj = Elastic2DMultiblockOpt(pars, inversionPars, order, k)
			default_arg('inversionPars', {'rho'});
			default_arg('order', 2);
			default_arg('k', []);

			% Unpack parameter struct
			xlim = pars.xlim;
			ylim = pars.ylim;
			lambda = pars.lambda;
			mu = pars.mu;
			rho = pars.rho;
			bc = pars.bc;
			sources = pars.sources;
			receivers = pars.receivers;
			T = pars.T;
			m = pars.m;
			ip = inversionPars;

			% If inverting for lambda or mu, extra optimization preparations are required
			if ismember('lambda', ip) || ismember('mu',ip)
				optFlag = true;
			else
				optFlag = false;
			end
			obj.optFlag = optFlag;

			% Create forward discretization object
			forwardDiscr = elastic.elasticDiscrInterface(m, order, xlim, ylim, lambda, mu, rho,...
                                               [], bc, sources, [], optFlag);
			obj.dim = forwardDiscr.dim;
			nBlocks = forwardDiscr.grid.nBlocks;

			if isa(lambda, 'function_handle')
                lambda = grid.evalOn(forwardDiscr.grid, lambda);
            end
            if isa(mu, 'function_handle')
                mu = grid.evalOn(forwardDiscr.grid, mu);
            end
            if isa(rho, 'function_handle')
                rho = grid.evalOn(forwardDiscr.grid, rho);
            end

            if iscell(lambda)
            	if isa(lambda{1}, 'function_handle')
            		lambda = multiblock.evalOn(forwardDiscr.grid, lambda);
            	end
            end
            if iscell(mu)
            	if isa(mu{1}, 'function_handle')
            		mu = multiblock.evalOn(forwardDiscr.grid, mu);
            	end
            end
            if iscell(rho)
            	if isa(rho{1}, 'function_handle')
            		rho = multiblock.evalOn(forwardDiscr.grid, rho);
            	end
            end
            pars.lambda = lambda;
            pars.mu = mu;
            pars.rho = rho;

		 	% Set zero initial data
			u0 = @(x,y) [0*x; 0*x];
			ut0 = @(x,y) [0*x; 0*x];
			forwardDiscr.v0 = grid.evalOn(forwardDiscr.grid, u0);
			forwardDiscr.v0t = grid.evalOn(forwardDiscr.grid, ut0);

			% Create discrete delta functions corresponding to sources and receivers.
			obj.sourceDeltas = forwardDiscr.generateDeltaFunctions(sources.x);
			obj.receiverDeltas = forwardDiscr.generateDeltaFunctions(receivers.x);

			% Choose time-step
			if isempty(k)
				CFL = 0.25;
				ip = inversionPars;
				% If inverting for strcture, use smaller time step
				if ismember('rho', ip) || ismember('lambda', ip) || ismember('mu', ip)
					CFL_margin = 2.5;
					CFL = CFL/CFL_margin;
				end
				k = forwardDiscr.getTimestep([], CFL);
				k = alignedTimestep(k, T);
			else
				k_limit = forwardDiscr.getTimestep([], 0.5);
				if(k > k_limit)
					k_1 = forwardDiscr.getTimestep([], 1);
					CFL = k / k_limit;
					warn_str = ['Chosen time-step is likely unstable with current parameters. '...
								'It corresponds to a CFL of ' num2str(CFL) '.' ];
					warning('ON','all');
					warning(warn_str);
					warning('OFF','all');
				end
			end
			obj.k = k;

			% Get time quadrature
			[ts, Nt] = forwardDiscr.getTimestepper([], T, obj.k);
			Ht_local = ts.getTimeStepQuadrature;
			obj.Ht = kron(ones(Nt,1), Ht_local);

			% Build time vector that includes stage points
			[~, ~, c, nStages] = ts.getTableau();
			tvec = zeros(1,Nt*nStages);
			dt = ts.k;
			for i = 1:Nt
				ii = (i-1)*nStages+1 : i*nStages;
				t_local = (i-1)*dt + c'*dt;
				tvec(ii) = t_local;
			end

			% Make sure that sources are discrete and not function handles
			for s = 1:numel(sources.x)
				for d = 1:obj.dim
					if isa(sources.g{s}{d}, 'function_handle')
						sources.g{s}{d} = sources.g{s}{d}(tvec);
					end
				end
			end
			pars.sources = sources;

			obj.B = [];
			obj.D1 = [];
			obj.diffOp_mu_1 = [];
			obj.diffOp_lambda_1 = [];

			% If optimizing for lambda or mu
			if ismember('lambda', ip) || ismember('mu',ip)

				nPar = length(lambda);

				% ============= Second derivatives ==================
	   			% B, used for adjoint optimization
	            % M_ij(a) + B_ijk*a_k

	            % Create empty cells
	            B = cell(obj.dim, 1);
	            for i = 1:obj.dim
	            	B{i} = cell(nPar, 1);
	            	for p = 1:nPar
	            		B{i}{p} = cell(nBlocks, nBlocks);
	            		for b = 1:nBlocks
	            			mb = prod(m{b});
	            			B{i}{p}{b, b} = sparse(mb, mb);
	            		end
	            	end
	            end

	            % Fill with single-block diffOps
	            for i = 1:obj.dim
	            	ind = 0;
	            	for j = 1:nBlocks
	            		m_local = prod(m{j});
	            		for k = 1:m_local
	            			B{i}{ind + k}{j, j} = forwardDiscr.diffOp.diffOps{j}.B{i}{k};
	            		end
	            		ind = ind + m_local;
	            	end
	            end

	            % Convert cell matrices to double
	            for i = 1:obj.dim
	            	for j = 1:nPar
	            		B{i}{j} = blockmatrix.toMatrix(B{i}{j});
	            	end
	            end
	            obj.B = B;
            	% ============================================= %

	            % ============ First derivatives ============= %
	            % D1 = forwardDiscr.diffOp.D1;
	            D1 = cell(obj.dim, 1);
	            D1{1} = multiblockOperator(forwardDiscr.grid, forwardDiscr.diffOp, 'D1{1}');
	            D1{2} = multiblockOperator(forwardDiscr.grid, forwardDiscr.diffOp, 'D1{2}');
				obj.D1 = D1;
				% ============================================= %

				% Construct derivatives of traction operators etc
				% by creating diffOp with lambda = 1, mu = 1;
				lambda_dummy = 0*lambda + 1;
				mu_dummy = 0*mu;
				toyDiscr = elastic.elasticDiscrInterface(m, order, xlim, ylim, lambda_dummy, mu_dummy, rho,...
                                               [], bc, sources);
                obj.diffOp_lambda_1 = toyDiscr.diffOp;
                mu_dummy = 0*mu + 1;
				lambda_dummy = 0*lambda;
				toyDiscr = elastic.elasticDiscrInterface(m, order, xlim, ylim, lambda_dummy, mu_dummy, rho,...
                                               [], bc, sources);
                obj.diffOp_mu_1 = toyDiscr.diffOp;
	        end

			% Set instance variables
			obj.forwardDiscr = forwardDiscr;
			obj.adjointDiscr = [];
			obj.sourcePoints = sources.x;
			obj.receiverPoints = receivers.x;
			obj.forwardReceiverRecordings = [];
			obj.adjointReceiverRecordings = [];
			obj.receiverData = [];

			obj.forwardField = struct;
			obj.forwardField.u = [];
			obj.forwardField.ut = [];
			obj.adjointField = struct;
			obj.adjointField.u = [];
			obj.adjointField.ut = [];

			obj.m = m;
			obj.order = order;
			obj.H = forwardDiscr.H;
			obj.T = T;
			obj.tvec = tvec;
			obj.pars = pars;
			obj.inversionPars = inversionPars;

		end


		function runForward(obj, plotFlag, T, saveOpts)
			default_arg('plotFlag', false);
			default_arg('T', obj.T);
			default_arg('saveOpts', []);

			discr = obj.forwardDiscr;
			receiverDeltas = obj.receiverDeltas;
			receiverPoints = obj.receiverPoints;
			[receiverData, waveField] = obj.runSimulation(discr, T, receiverPoints, receiverDeltas, plotFlag, saveOpts);
			obj.forwardReceiverRecordings = receiverData;
			obj.forwardField = waveField;
		end

		function runAdjoint(obj, plotFlag, T, saveOpts)
			default_arg('plotFlag', false);
			default_arg('T', obj.T);
			default_arg('saveOpts', []);

			discr = obj.adjointDiscr;
			receiverDeltas = obj.sourceDeltas;
			receiverPoints = obj.sourcePoints;

			[receiverData, waveField] = obj.runSimulation(discr, T, receiverPoints, receiverDeltas, plotFlag, saveOpts);

			nSources = numel(obj.sourcePoints);
			for i = 1:nSources
				for d = 1:obj.dim
					receiverData{i,d} = rot90(receiverData{i,d},2);
				end
			end
			obj.adjointReceiverRecordings = receiverData;
			waveField.u = fliplr(waveField.u);
			waveField.ut = -fliplr(waveField.ut);
			obj.adjointField = waveField;
		end

		function [receiverRecordings, waveField] = runSimulation(obj, discr, T, ...
												receiverPoints, receiverDeltas, plotFlag, saveOpts)
			default_arg('saveOpts', []);
			saveData = struct;
			saveData.t = [];
			saveData.u = [];

			[ts, N] = discr.getTimestepper([], T, obj.k);
			[~, ~, ~, nStages] = ts.getTableau();
			H = obj.H;

			% For storing data at receivers.
			nReceivers = numel(receiverPoints);
			receiverRecordings = cell(nReceivers, obj.dim);
			for i = 1:nReceivers
				for d = 1:obj.dim
					receiverRecordings{i,d} = sparse(1, length(obj.Ht));
				end
			end

			% For storing full field
			[m, ~] = size(discr.diffOp.D);
			waveField = struct;
			waveField.u = sparse(m, length(obj.Ht));
			waveField.ut = sparse(m, length(obj.Ht));

			% Setup plot
			if plotFlag
				[update, fig] = discr.setupPlot();
				r = struct;
			end

			for i = 1:N
				ts.step();
				[~, ~, V, ~, Vt] = ts.getV;
				ii = (i-1)*nStages+1 : i*nStages;

				% Store at receivers
				for j = 1:nReceivers
					for d = 1:obj.dim
						receiverRecordings{j,d}(ii) = (H*receiverDeltas{j,d})'*V;
					end
				end

				% Store time-derivative if inverting for rho
				ip = obj.inversionPars;
				if ismember('rho', ip)
					waveField.ut(:, ii) = Vt;
				end

				% Store solution if inverting for lambda or mu
				ip = obj.inversionPars;
				if ismember('lambda', ip) || ismember('mu', ip)
					waveField.u(:, ii) = V;
				end

				% Plot
				if( plotFlag && (mod(i,10) == 0) )
					[v, t] = ts.getV;
					r.v = v;
					r.t = t;
					update(r);
					drawnow;
				end

				% Save solution if requested
				if ~isempty(saveOpts) && (mod(i, saveOpts.frameSpacing) == 0)
					[v, t] = ts.getV;
					saveData.u = [saveData.u, v];
					saveData.t = [saveData.t, t];
				end
			end

			if ~isempty(saveOpts)
				saveData.g = obj.forwardDiscr.grid;
				saveData.Div = obj.forwardDiscr.Div;
				saveData.LAMBDA = obj.forwardDiscr.LAMBDA;
				save(saveOpts.filename, 'saveData','-v7.3','-nocompression');
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
				switch parName
				case 'source'
					obj.pars.sources.g = pars.sources.g;
				case {'rho', 'lambda', 'mu'}
					obj.pars.(parName) = pars.(parName);
				end
			end
		end

		% Converts inversion parameters from native format to one column vector
		% parsNative -- parameters in native fomat
		function parsVec = parametersNativeToVector(obj, parsNative)
			pars = cell(numel(obj.inversionPars), 1);
			for i = 1:numel(obj.inversionPars)
				parName = obj.inversionPars{i};
				switch parName
				case 'source'
					nSources = numel(obj.sourcePoints);
					vec = [];
					for s = 1:nSources
						for d = 1:obj.dim
							vec = [vec; parsNative.sources.g{s}{d}];
						end
					end
					pars{i} = vec;
				case {'rho','lambda','mu'}
					pars{i} = parsNative.(parName);
				end
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
				switch parName
				case 'source'
					nSources = numel(obj.sourcePoints);
					m = length(obj.pars.sources.g{1}{1});
					parsNative.sources = struct;
					parsNative.sources.g = cell(nSources, obj.dim);
					for s = 1:nSources
						for d = 1:obj.dim
							parsNative.sources.g{s}{d} = parsVec(index : index+m-1);
							index = index + m;
						end
					end
				case {'rho','lambda','mu'}
					m = length(obj.pars.(parName));
					parsNative.(parName) = parsVec(index : index+m-1);
					index = index + m;
				end
			end
		end

		% Converts gradient from native format to one column vector
		% gradNative -- gradient in native fomat
		function gradVec = gradientNativeToVector(obj, gradNative)
			grad = cell(numel(obj.inversionPars), 1);
			for i = 1:numel(obj.inversionPars)
				parName = obj.inversionPars{i};
				switch parName
				case 'source'
					nSources = numel(obj.sourcePoints);
					vec = [];
					for s = 1:nSources
						for d = 1:obj.dim
							vec = [vec; gradNative{i}{s,d}];
						end
					end
					grad{i} = vec;
				case {'rho','lambda','mu'}
					grad{i} = gradNative{i};
				end
			end
			gradVec = cell2mat(grad);
		end

		% Converts gradient vector to native format
		% gradVec -- gradient vector
		function gradNative = gradientVectorToNative(obj, gradVec)
			gradNative = cell(numel(obj.inversionPars), 1);
			index = 1;
			for i = 1:numel(obj.inversionPars)
				parName = obj.inversionPars{i};
				switch parName
				case 'source'
					nSources = numel(obj.sourcePoints);
					m = length(obj.pars.sources.g{1}{1});
					gradNative{i} = cell(nSources, obj.dim);
					for s = 1:nSources
						for d = 1:obj.dim
							gradNative{i}{s,d} = gradVec(index : index+m-1);
							index = index + m;
						end
					end
				case {'rho','lambda','mu'}
					m = length(obj.pars.(parName));
					gradNative{i} = gradVec(index : index+m-1);
					index = index + m;
				end
			end
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
			default_arg('pars', obj.pars);

			% Create discr
			obj.forwardDiscr = obj.createElasticDiscr(pars);

		end

		function updateAdjointDiscr(obj)

			pars = obj.pars;

			% Make sure that the adjoint sources are the receivers with the correct data
			sources = pars.receivers;
			data = obj.receiverData;
			approx = obj.forwardReceiverRecordings;
			nReceivers = numel(obj.receiverPoints);
			for rec = 1:nReceivers
				for comp = 1:obj.dim
					err = approx{rec,comp} - data{rec,comp};
					sources.g{rec}{comp} = rot90(err,2);
				end
			end
			pars.sources = sources;

			% Create discr
			obj.adjointDiscr = obj.createElasticDiscr(pars);
		end

		% Helper function that is used to create both forward and adjoint discr.
		function discr = createElasticDiscr(obj, pars)
			order = obj.order;

			% Unpack parameter struct
			xlim = pars.xlim;
			ylim = pars.ylim;
			lambda = pars.lambda;
			mu = pars.mu;
			rho = pars.rho;
			bc = pars.bc;
			sources = pars.sources;
			T = pars.T;
			m = pars.m;

			% Create discretization object
			discr = elastic.elasticDiscrInterface(m, order, xlim, ylim, lambda, mu, rho,...
                                               [], bc, sources, [], obj.optFlag);

			% Check if time step is stable
			k_limit = discr.getTimestep([], 0.5);
			if(obj.k > k_limit)
				k_1 = discr.getTimestep([], 1);
				CFL = obj.k / k_limit;
				warn_str = ['Chosen time-step is likely unstable with current parameters. '...
							'It corresponds to a CFL of ' num2str(CFL) '.' ];
				warning('ON','all');
				warning(warn_str);
				warning('OFF','all');
			end

		 	% Set zero initial data
			u0 = @(x,y) [0*x; 0*x];
			ut0 = @(x,y) [0*x; 0*x];
			discr.v0 = grid.evalOn(discr.grid, u0);
			discr.v0t = grid.evalOn(discr.grid, ut0);

		end

		function grad = computeGradient(obj)
			obj.runForward(false);
			obj.updateAdjointDiscr();
			obj.runAdjoint(false);
			grad = obj.gradientFormula();
		end

		function grad = gradientFormula(obj)

			grad = cell(numel(obj.inversionPars));

			for i = 1:numel(obj.inversionPars)
				switch obj.inversionPars{i}

				case 'source'
					nSources = numel(obj.sourcePoints);
					grad{i} = cell(nSources, obj.dim);
					for s = 1:nSources
						for comp = 1:obj.dim
							grad{i}{s, comp} = obj.adjointReceiverRecordings{s, comp}.*(obj.Ht');
						end
					end

				case 'rho'
					E = obj.forwardDiscr.Ecomp;
					grad{i} = 0;
					% Sum over components
					for j = 1:obj.dim
						grad{i} = grad{i} + obj.timeIntegration((E{j}'*obj.forwardField.ut) .* (E{j}'*obj.adjointField.ut));
					end
					% Multiply by space quadrature weights
					grad{i} = (E{1}'*diag(obj.H)).*grad{i};

				case 'lambda'
					D1 = obj.D1;

					E = obj.forwardDiscr.Ecomp;
					for k = 1:obj.dim
						E{k} = transpose(E{k});
					end
					grad{i} = 0;

            		H_kron = obj.H;
            		H = E{1}*H_kron*transpose(E{1});
            		H_row_vec = diag(H)';

            		psi = obj.adjointField.u;
            		u = obj.forwardField.u;

					% Mixed derivatives
					for j = 1:obj.dim
						% For k =/= j
						for k = 1:j-1
							grad{i} = grad{i} - (D1{j}*(E{j}*psi)) .* (D1{k}*(E{k}*u));
						end
						for k = j+1:obj.dim
							grad{i} = grad{i} - (D1{j}*(E{j}*psi)) .* (D1{k}*(E{k}*u));
						end
					end
					% Integrate in time
					grad{i} = obj.timeIntegration(grad{i});

					% Multiply by space quadrature weights
					grad{i} = H*grad{i};

					% Non-mixed derivatives
					nPar = length(obj.pars.lambda);
					for j = 1:obj.dim

						Uj = E{j}*u;
						PSIj = E{j}*psi;

						for q = 1:nPar
							Bjq = obj.B{j}{q};
							grad_vec = H_row_vec*(PSIj.*(Bjq*Uj));

							% Integrate in time
							grad{i}(q) = grad{i}(q) - obj.timeIntegration(grad_vec);
						end

					end

				case 'mu'
					D1 = obj.D1;

					E = obj.forwardDiscr.Ecomp;
					for k = 1:obj.dim
						E{k} = transpose(E{k});
					end
					grad{i} = 0;

            		H_kron = obj.H;
            		H = E{1}*H_kron*transpose(E{1});
            		H_row_vec = diag(H)';

            		psi = obj.adjointField.u;
            		u = obj.forwardField.u;

					% Mixed derivatives
					for j = 1:obj.dim
						% For k =/= j
						for k = 1:j-1
							grad{i} = grad{i} - (D1{k}*(E{j}*psi)) .* (D1{j}*(E{k}*u));
						end
						for k = j+1:obj.dim
							grad{i} = grad{i} - (D1{k}*(E{j}*psi)) .* (D1{j}*(E{k}*u));
						end
					end
					% Integrate in time
					grad{i} = obj.timeIntegration(grad{i});

					% Multiply by space quadrature weights
					grad{i} = H*grad{i};

					% Non-mixed derivatives
					nPar = length(obj.pars.mu);
					for j = 1:obj.dim
						Uj = E{j}*u;
						PSIj = E{j}*psi;

						for q = 1:nPar
							%----- Terms from d_j mu d_i u_j
							Bjq = obj.B{j}{q};
							grad_vec = H_row_vec*(PSIj.*(Bjq*Uj));

							% Integrate in time
							grad{i}(q) = grad{i}(q) - obj.timeIntegration(grad_vec);
							%------------------------------

							%------ Terms from d_j mu d_j u_i --------
							grad_vec = 0;
							for l = 1:obj.dim
								Blq = obj.B{l}{q};
								grad_vec = grad_vec + H_row_vec*(PSIj.*(Blq*Uj));
							end
							% Integrate in time
							grad{i}(q) = grad{i}(q) - obj.timeIntegration(grad_vec);
							%-------------------------------------
						end
					end
				end

				if ismember(obj.inversionPars{i}, {'lambda','mu'})

					switch obj.inversionPars{i}
					case 'mu'
						diffOp = obj.diffOp_mu_1;
					case 'lambda'
						diffOp = obj.diffOp_lambda_1;
					end

					% Add correction terms for displacement BC, for lambda and mu
					for b = 1:numel(obj.pars.bc)
						boundary = obj.pars.bc{b}.boundary;
						comp = obj.pars.bc{b}.type{1};
						type = obj.pars.bc{b}.type{2};

						if ismember(type, {'D','d','dirichlet','Dirichlet'})

							% Get boundary operators
							e 		= diffOp.getBoundaryOperator('e', boundary);
							e_c 	= diffOp.getBoundaryOperator(['e' num2str(comp)], boundary);
							H_b 	= diffOp.getBoundaryQuadrature(boundary);
							tau_c 	= diffOp.getBoundaryOperator(['tau' num2str(comp)], boundary);
							alpha_c = diffOp.getBoundaryOperator(['alpha' num2str(comp)], boundary);
							tuning = 1.2;

							% Create e for scalar field
							e_scalar = (e'*E{1}')';
                            e_scalar = e_scalar(:, 1:obj.dim:end);

							% Boundary quadrature for scalar field
							H_b = e_scalar'*E{1}*e*H_b*e'*transpose(E{1})*e_scalar;

							% Add gradient contributions
							grad_temp = (e_c'*psi) .* (H_b*tau_c'*u);
							grad_temp = grad_temp + (tau_c'*psi) .* (H_b*e_c'*u);
							grad_temp = grad_temp - tuning*(e_c'*psi) .* (H_b*alpha_c'*u);

							% Integrate in time
							grad_temp = obj.timeIntegration(grad_temp);

							% Project and add to gradient
							grad{i} = grad{i} + e_scalar * grad_temp;
						end
					end

					% Add interface correction terms
					conn = obj.forwardDiscr.grid.connections;
					nBlocks = obj.forwardDiscr.grid.nBlocks;
					for b1 = 1:nBlocks
						for b2 = 1:nBlocks
							if ~isempty(conn{b1, b2})
								boundary{1} = {b1, conn{b1, b2}{1}};
								boundary{2} = {b2, conn{b1, b2}{2}};

								% Build boundary operators for both blocks
								for blockId = 1:2
									H_b_cell{blockId} = diffOp.getBoundaryQuadrature(boundary{blockId});
									alpha_cell{blockId} = diffOp.getBoundaryOperator('alpha', boundary{blockId});
									e_cell{blockId} = diffOp.getBoundaryOperator('e', boundary{blockId});
									tau_cell{blockId} = diffOp.getBoundaryOperator('tau', boundary{blockId});
	                                e_scalar_cell{blockId} = (e_cell{blockId}'*E{1}')';
	                                e_scalar_cell{blockId} = e_scalar_cell{blockId}(:, 1:obj.dim:end);
								end

								% Loop over both sides of interface
								for local_block = 1:2

									% Assign operators corresponding to own block and neighbor block
									switch local_block
									case 1
										e = e_cell{1};
										e_nei = e_cell{2};
                                        e_scalar = e_scalar_cell{1};

										H_b = H_b_cell{1};
										H_b_nei = H_b_cell{2};

										tau = tau_cell{1};
										alpha = alpha_cell{1};
									case 2
										e = e_cell{2};
										e_nei = e_cell{1};
										e_scalar = e_scalar_cell{2};

										H_b = H_b_cell{2};
										H_b_nei = H_b_cell{1};

										tau = tau_cell{2};
										alpha = alpha_cell{2};
									end
									tuning = 1.2;

									% Boundary terms from operator
									grad_local_tot = (e'*psi) .* (H_b*tau'*u);

									% Interface SATs on own side
									grad_local_tot = grad_local_tot + 1/2*(tau'*psi) .* (H_b * (e'-e_nei') * u);
									grad_local_tot = grad_local_tot - 1/2*(e'*psi) .* (H_b * tau' * u);
									grad_local_tot = grad_local_tot - tuning*1/4*(e'*psi) .* (H_b*alpha'*e*(e' - e_nei')*u);

									% Interface SATs on opposite side
									grad_local_tot = grad_local_tot - tuning*1/4*(e_nei'*psi) .* (H_b_nei*alpha'*e*(e_nei' - e')*u);
									grad_local_tot = grad_local_tot - 1/2*(e_nei'*psi) .* (H_b_nei * tau' * u);

									grad_local = 0;
									for k = 1:obj.dim
										grad_local = grad_local + e_scalar'*E{k}*e*grad_local_tot;
									end

                                    % Integrate in time
                                    grad_local = obj.timeIntegration(grad_local);

                                    % Project and add to gradient
                                    grad{i} = grad{i} + e_scalar*grad_local;
								end
							end
						end
					end
				end
			end
		end

		% Should deltaG be different for different kinds of parameters?'
		function grad = computeGradientFD(obj, deltaG)


			% Compute misfit with current parameters
			pars = obj.pars;
			obj.updateForwardDiscr(pars);
			obj.runForward();
			M0 = obj.computeMisfit();

			grad = cell(numel(obj.inversionPars), 1);

			for i = 1:numel(obj.inversionPars)

				switch obj.inversionPars{i}

				case {'rho', 'lambda', 'mu' }
					parName = obj.inversionPars{i};
					par = getfield(pars, parName);
					grad{i} = zeros(size(par));

					for j = 1:length(par);
						par_temp = par;
						par_temp(j) = par_temp(j) + deltaG;
						pars = setfield(pars, parName, par_temp);
						obj.updateForwardDiscr(pars);
						obj.runForward(false);
						M = obj.computeMisfit();
						dM_dg = (M-M0)/deltaG;
						grad{i}(j) = dM_dg;
					end

				case 'source'
					tvec = obj.tvec;
					Nt = length(tvec);
					pars = obj.pars;
					sources = pars.sources;
					g = sources.g;

					nSources = numel(obj.sourcePoints);
					grad{i} = cell(nSources, obj.dim);
					for s = 1:nSources
						for d = 1:obj.dim

							% Evaluate source function handles if needed
							if isa(g{s}{d}, 'function_handle')
								g{s}{d} = g{s}{d}(tvec);
							end

							% Allocate gradient
							grad{i}{s, d} = zeros(1, Nt);
						end
					end

					for j = 1:Nt
						for s = 1:nSources
							for d = 1:obj.dim
								g_temp = g;
								g_temp{s}{d}(j) = g_temp{s}{d}(j) + deltaG;
								pars.sources.g = g_temp;
								obj.updateForwardDiscr(pars);
								obj.runForward(false);
								M = obj.computeMisfit();
								dM_dg = (M-M0)/deltaG;
								grad{i}{s,d}(j) = dM_dg;
							end
						end
					end
				end
			end
		end

		% v is assumed to be a matrix where each column
		% corresponds to a time point
		function vInt = timeIntegration(obj, v)
			vInt = v*obj.Ht;
		end

		function M = computeMisfit(obj)
			data = obj.receiverData;
			approx = obj.forwardReceiverRecordings;
			nReceivers = numel(obj.receiverPoints);

			M = 0;
			for i = 1:nReceivers
				for d = 1:obj.dim
					err2 = 1/2*(approx{i,d} - data{i,d}).^2;
					errInt = obj.timeIntegration(err2);
					M = M + errInt;
				end
			end
		end

		% Useful for finding an approximate line search interval
		function relnorm = gradientNorm(obj, grad)

			gradnorm = zeros(1, numel(obj.inversionPars) );
			relnorm = zeros(1, numel(obj.inversionPars) );
			parnorm = zeros(1, numel(obj.inversionPars) );

			for i = 1:numel(obj.inversionPars)

				switch obj.inversionPars{i}
				case 'source'
					nSources = numel(obj.sourcePoints);
					pars = obj.pars.sources.g;
					for s = 1:nSources
						for d = 1:obj.dim
							parnorm(i) = parnorm(i) + norm(pars{s}{d});
							gradnorm(i) = gradnorm(i) + norm(grad{i}{s,d});
						end
					end
					relnorm(i) = gradnorm(i)/parnorm(i);
				case 'rho'
					relnorm(i) = norm(grad{i}./obj.pars.rho, inf);
				case 'lambda'
					relnorm(i) = norm(grad{i}./obj.pars.lambda, inf);
				case 'mu'
					relnorm(i) = norm(grad{i}./obj.pars.mu, inf);
				end
			end
			relnorm = max(relnorm);

		end

		function [data_handle, model_handle] = plotReceiverSignals(obj, receiver, component, fontsize)
			default_arg('receiver', 1);
			default_arg('component', 1);

			data = obj.receiverData;
			approx = obj.forwardReceiverRecordings;

			data_handle = plot(obj.tvec, data{receiver,component}, 'linewidth', 2);
			hold on
			model_handle = plot(obj.tvec, approx{receiver,component}, 'linewidth', 2);
			xlabel('t', 'fontsize', fontsize)
			xlim([0, obj.T])
			ylabel(sprintf('Rec %d, comp %d',receiver,component), 'fontsize', fontsize);
			legend({'Data', 'Model'}, 'fontsize', fontsize, 'Location', 'northwest')
			hold off;
		end

		% Plots the parameter parName for the sets in parsets
		% Example:
		% parName = 'source' (or 'rho')
		% parsets = {parset1, parset2, ...}
		function sur = plotParameterComparison(obj, parName, parsets, sur, fontsize)
			default_arg('sur',[]);

			switch parName
			case 'source'
				source = 1;
				comp = 1;

				% Plot
				for i = 1:numel(parsets)
					g = parsets{i}.sources.g{source}{comp};
					if isa(g, 'function_handle')
						g = g(obj.tvec);
					end
					plot(obj.tvec(1:4:end), g(1:4:end), 'linewidth', 2);
					hold on;
				end
				hold off;
				ylimpar = ylim;
				L = ylimpar(2)-ylimpar(1);
				ylimpar = [ylimpar(1)-0.1*L, ylimpar(2)+0.1*L];
				ylim(ylimpar);
				xlabel('t', 'fontsize', fontsize)
				ylabel('Source time function', 'fontsize', fontsize)
				legend('True', 'Current guess', 'Initial guess', 'fontsize', fontsize)

			case {'rho','lambda','mu'}
				g = obj.forwardDiscr.grid;

				% Assume that there are two parameter sets and plot the relative difference
				par = cell(2,1);
				for i = 1:2
					par{i} = getfield(parsets{i}, parName);
					if isa(par{i}, 'function_handle')
						par{i} = grid.evalOn(g, par{i});
					end
				end
				dpar = abs(par{2}-par{1})./par{1};

				if isempty(sur)
					sur = multiblock.Surface(g, dpar);
					view(0,90);
					shading interp;
					colorbar
					title(['Relative error in ' parName], 'fontsize', fontsize)
					xlabel('x', 'fontsize', fontsize)
					ylabel('y', 'fontsize', fontsize)
					caxis([0, 1.2])
				else
					sur.ZData = dpar;
					sur.CData = dpar;
					caxis([0, 1.2])
				end


			end
		end

		function [upPar, upRec, upMis, figure_handle] = setupPlot(obj, trueParset, initialParset, M, M0, iter)

			%-- Plot settings ---
			fontsize = 16;
			scrsz = get(0,'ScreenSize');
			figure_handle = figure('Position',[0.1*scrsz(3) 0.05*scrsz(4) 0.65*scrsz(3) 0.8*scrsz(4)]);
			%--------------------

			parfig = subplot(3,1,1);
			sur = obj.plotParameterComparison(obj.inversionPars{1} ,{trueParset, obj.pars, initialParset}, [], fontsize);

			subplot(3,1,2);
			[~, model_handle] = obj.plotReceiverSignals([],[],fontsize);
			ylimrec = 1.2*ylim;
			L = ylimrec(2)-ylimrec(1);
			ylimrec = [ylimrec(1)-0.1*L, ylimrec(2)+0.1*L];
			ylim(ylimrec);

			subplot(3,1,3);
			misfit = semilogy(M/M0, 'linewidth', 2);
			xlabel('Iterations', 'fontsize', fontsize)
			ylabel('Misfit', 'fontsize', fontsize)
			xlim([1, iter]);

			function updateParameter(object, sur)
				axes(parfig);
				object.plotParameterComparison(obj.inversionPars{1} ,{trueParset, object.pars, initialParset}, sur, fontsize);
			end

			function updateReceiver(object, receiver, component)
				default_arg('receiver', 1);
				default_arg('component', 1);

				approx = object.forwardReceiverRecordings{receiver, component};
				model_handle.YData = approx;
			end

			function updateMisfit(M)
				misfit.YData = M/M0;
			end

			upPar = @(object)updateParameter(object, sur);
			upRec = @updateReceiver;
			upMis = @updateMisfit;

		end

		function plotGradient(obj, grad)

			% Plot first gradient component
			parName = obj.inversionPars{1};
			switch parName
			case 'source'

			case {'rho','lambda','mu'}
				g = obj.forwardDiscr.grid;

				multiblock.Surface(g, grad{1});
				view(0,90);
				shading interp;
				colorbar
				ylabel(['Gradient with respect to ' parName])
			end
		end

	end

end