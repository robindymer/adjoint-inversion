classdef ExplicitRungeKuttaDiscreteData < time.Timestepper
    properties
        D
        S_cont           % Function handle @(t, v) for continuous time-dependent data
        S_discr          % Function handle @(idx_t, v) for discrete time-depedent data
                         % One column per stage
        k
        t
        v
        m
        n
        order
        a, b, c, s  % Butcher tableau
        K           % Stage rates
        U           % Stage approximations
        T           % Stage times
    end


    methods
        function obj = ExplicitRungeKuttaDiscreteData(D, S_cont, S_discr, k, t0, v0, order)
            default_arg('order', 4);
            default_arg('S_cont', []);
            default_arg('S_discr', []);

            obj.D = D;
            obj.S_cont = S_cont;
            obj.S_discr = S_discr;
            obj.k = k;
            obj.t = t0;
            obj.v = v0;
            obj.m = length(v0);
            obj.n = 0;
            obj.order = order;

            switch order
            case 3
                [obj.a, obj.b, ~, obj.c, obj.s] = time.rkparameters.rk32();
            case 4
                [obj.a, obj.b, obj.c, obj.s] = time.rkparameters.rk4();
            case 5
                [obj.a, obj.b, ~, obj.c, obj.s] = time.rkparameters.rk54();
            case 6
                [obj.a, obj.b, obj.c, obj.s] = time.rkparameters.rk6();
            otherwise
                error('That RK method is not available');
            end

            obj.K = zeros(obj.m, obj.s);
            obj.U = zeros(obj.m, obj.s);

        end

        function [v,t,U,T,K] = getV(obj)
            v = obj.v;
            t = obj.t;
            U = obj.U; % Stage approximations in previous time step.
            T = obj.T; % Stage times in previous time step.
            K = obj.K; % Stage rates in previous time step.
        end

        function [a,b,c,s] = getTableau(obj)
            a = obj.a;
            b = obj.b;
            c = obj.c;
            s = obj.s;
        end

        % Returns quadrature weights for stages in one time step
        function quadWeights = getTimeStepQuadrature(obj)
            [~, b] = obj.getTableau();
            quadWeights = obj.k*b;
        end

        function obj = step(obj)
            v = obj.v;
            a = obj.a;
            b = obj.b;
            c = obj.c;
            s = obj.s;
            S_cont = obj.S_cont;
            S_discr = obj.S_discr;
            dt = obj.k;
            K = obj.K;
            U = obj.U;
            D = obj.D;

            for i = 1:s
                U(:,i) = v;
                for j = 1:i-1
                    U(:,i) = U(:,i) + dt*a(i,j)*K(:,j);
                end

                K(:,i) = D*U(:,i);
                obj.T(i) = obj.t + c(i)*dt;

                % Data from continuous function and discrete time-points.
                if ~isempty(S_cont)
                    K(:,i) = K(:,i) + S_cont(obj.T(i),U(:,i));
                end
                if ~isempty(S_discr)
                    idx_t = obj.n*s + i;
                    K(:,i) = K(:,i) + S_discr(idx_t,U(:,i));
                end

            end

            obj.v = v + dt*K*b;
            obj.t = obj.t + dt;
            obj.n = obj.n + 1;
            obj.U = U;
            obj.K = K;
        end
    end


    methods (Static)
        function k = getTimeStep(lambda, order)
            default_arg('order', 4);
            switch order
            case 4
                k = time.rk4.get_rk4_time_step(lambda);
            otherwise
                error('Time-step function not available for this order');
            end
        end
    end

end