classdef EmbeddedRungeKutta < time.Timestepper
    properties
        F                       % RHS function F(t,v)
        k                       % time step
        t                       % current time
        v                       % current solution vector
        m                       % size of solution vector
        n                       % current iteration                  
        order                   % order of the timestepper
        a, b, bStar, c, s       % Butcher tableau
        K                       % Stage rates
        U                       % Stage approximations
        T                       % Stage times
        relTol                  % relative tolerance used
        errorCheckHandle        % Handle for computing the error used for accepting updates
        constraintCheckHandle   % Handle for checking constraint violations in the time-stepping, e.g positivity.
        reportRetry             % Flag determining whether retrys should be reported.
    end


    methods
        function obj = EmbeddedRungeKutta(F, k, t0, v0, order, relTol, errorCheckCallback, constraintCheckCallback, reportRetry)
            default_arg('order', 3);
            default_arg('relTol', 1e-3);
            default_arg('errorCheckCallback', []);
            default_arg('constraintCheckCallback', []);
            default_arg('reportRetry', false);
            assert(relTol > 0, 'Error tolerance must be positive.');

            obj.F = F;
            obj.k = k;
            obj.t = t0;
            obj.v = v0;
            obj.m = length(v0);
            obj.n = 0;
            obj.order = order;
            obj.relTol = relTol;
            obj.reportRetry = reportRetry;

            switch order
            case 3
                [obj.a, obj.b, obj.bStar, obj.c, obj.s] = time.rkparameters.rk32();
            case 5
                [obj.a, obj.b, obj.bStar, obj.c, obj.s] = time.rkparameters.rk54();
            otherwise
                error('That embedded RK method is not available');
            end

            obj.K = zeros(obj.m, obj.s);
            obj.U = zeros(obj.m, obj.s);

            if isempty(errorCheckCallback)
                % Default: Use maximum norm error
                obj.errorCheckHandle = @(vNew, vStar) max(abs(vNew - vStar)./(abs(vNew)+1e-16));
            else
                obj.errorCheckHandle = errorCheckCallback;
            end

            if isempty(constraintCheckCallback)
                % Default: No user constraint
                obj.constraintCheckHandle = [];
            else
                obj.constraintCheckHandle = constraintCheckCallback;
            end

        end

        function [v,t,U,T,K] = getV(obj)
            v = obj.v;
            t = obj.t;
            U = obj.U; % Stage approximations in previous time step.
            T = obj.T; % Stage times in previous time step.
            K = obj.K; % Stage rates in previous time step.
        end

        function [a,b,c,s,bStar] = getTableau(obj)
            a = obj.a;
            b = obj.b;
            c = obj.c;
            s = obj.s;
            bStar = obj.bStar;
        end

        % Returns quadrature weights for stages in one time step
        % with stepsize k
        function quadWeights = getTimeStepQuadrature(obj)
            quadWeights = obj.k*obj.b;
        end

        function obj = step(obj)
            v = obj.v;
            a = obj.a;
            b = obj.b;
            bStar = obj.bStar;
            c = obj.c;
            s = obj.s;
            K = obj.K;
            U = obj.U;
            F = obj.F;

            % Adaptive step. Retries until error is sufficiently small
            tol = obj.relTol;
            err = 2*obj.relTol;
            count = 1;

            while err > tol
                dt = obj.k;
                if isnan(dt)
                    error('dt = NaN');
                end
                if dt < 1e-30
                    error('dt too small');
                end

                for i = 1:s
                    U(:,i) = v;
                    for j = 1:i-1
                        U(:,i) = U(:,i) + dt*a(i,j)*K(:,j);
                    end

                    [violation, violationstr] = obj.checkViolation(U(:,i));
                    if violation
                        break;
                    end

                    obj.T(i) = obj.t + c(i)*dt;
                    K(:,i) = F(obj.T(i), U(:,i));

                    [violation, violationstr] = obj.checkViolation(K(:,i));
                    if violation
                        break;
                    end
                end

                vNew = v + dt*K*b;
                vStar = v + dt*K*bStar;

                err = obj.errorCheckHandle(vNew, vStar);
                
                % Update timestep and count try if retaking the same step
                if violation
                    % Violation encountered. Notify user and update timestep
                    if obj.reportRetry
                        fprintf('Failed try %d, due to %s, dt = %0.1e, t = %0.8e\n',...
                                count, violationstr, dt, obj.t);
                    end
                    err = 2*tol;
                    obj.k = dt/2;
                    count = count + 1;
                elseif err > tol 
                    % Tolerance exceeded. Notify user and update timestep
                    if obj.reportRetry
                        fprintf('Failed try %d, err = %0.1e, dt = %0.1e, t = %0.8e\n',...
                                count, err, dt, obj.t);
                    end
                    obj.k = 0.9*dt*(tol/err)^(1/obj.order);
                    count = count + 1;
                elseif err == 0
                    obj.k = dt*2;
                else
                    obj.k = 0.9*dt*(tol/err)^(1/obj.order);
                end
            end

            obj.v = vNew;
            obj.t = obj.t + dt;
            obj.n = obj.n + 1;
            obj.U = U;
            obj.K = K;
        end

        function [isviolated, violationstr] = checkViolation(obj, U)
            % Check if positivity constraints are violated
            isviolated = false;
            violationstr = '';
            
            % Check user constraints
            if ~isempty(obj.constraintCheckHandle)
                [isviolated, violationstr] = obj.constraintCheckHandle(U);
                return
            end

            % Check for NaNs
            if any(isnan(U))
                isviolated = true;
                violationstr = 'NaN violation';
                return;
            end

            % Check for infs
            if any(isinf(U))
                isviolated = true;
                violationstr = 'inf violation';
                return;
            end
        end
    end

end