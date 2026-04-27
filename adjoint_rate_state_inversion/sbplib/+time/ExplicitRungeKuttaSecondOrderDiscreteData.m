classdef ExplicitRungeKuttaSecondOrderDiscreteData < time.Timestepper
    properties
        k
        t
        w
        m
        D
        E
        M
        C_cont % Continuous part (function handle) of forcing on first order form.
        C_discr% Discrete part (matrix) of forcing on first order form.
        n
        order
        tsImplementation % Time stepper object, RK first order form,
                         % which this wraps around.
    end


    methods
        % Solves u_tt = Du + Eu_t + S by
        % Rewriting on first order form:
        %   w_t = M*w + C(t,w)
        % where
        %   M = [
        %      0, I;
        %      D, E;
        %   ]
        % and
        %   C(t,w) = [
        %      0;
        %      S(t,w)
        %   ]
        % D, E, should be matrices (or empty for zero)
        % They can also be omitted by setting them equal to the empty matrix.
        % S = S_cont + S_discr, where S_cont is a function handle
        % S_discr a function handle or a matrix of data vectors, with one column per stage.
        function obj = ExplicitRungeKuttaSecondOrderDiscreteData(D, E, S_cont, S_discr, k, t0, v0, v0t, order)
            default_arg('order', 4);
            default_arg('S_cont', []);
            default_arg('S_discr', []);
            obj.D = D;
            obj.E = E;
            obj.m = length(v0);
            obj.n = 0;

            default_arg('D', sparse(obj.m, obj.m) );
            default_arg('E', sparse(obj.m, obj.m) );

            obj.k = k;
            obj.t = t0;
            obj.w = [v0; v0t];

            I = speye(obj.m);
            O = sparse(obj.m,obj.m);

            obj.M = [
                O, I;
                D, E;
            ];

            % Build C_cont
            if ~isempty(S_cont)
                % Ensure that S_cont is of the correct form for
                % ExplicitRungeKuttaDiscreteData
                switch nargin(S_cont)
                case 1 % S_cont(t)
                    S_cont = @(t, v ,vt) S_cont(t);
                case 2 % S_cont(t,v)
                    S_cont = @(t, v, vt) S_cont(t, v);
                otherwise
                    % S_cont(t,v,vt) or error
                    assert(nargin(S_cont) == 3,'Incorrect number or arguments in S_cont');
                end
                obj.C_cont = @(t,w)[
                    sparse(obj.m, 1);
                    S_cont(t, w(1:end/2), w(end/2+1:end))
                            ];
            else
                obj.C_cont = [];
            end

            % Build C_discr
            if ~isempty(S_discr)
                % Ensure that S_discr is of the correct form for
                % ExplicitRungeKuttaDiscreteData
                if isa(S_discr,'function_handle')
                    switch nargin(S_discr)
                    case 1 % S_discr(t)
                        S_discr = @(idx_t, v, vt) S_discr(idx_t);
                    case 2 % S_discr(t,v)
                        S_discr = @(idx_t, v, vt) S_discr(idx_t, v);
                    otherwise
                        % S_discr(t,v,vt) or error
                        assert(nargin(S_discr) == 3,'Incorrect number or arguments in S_discr');
                    end
                else % S_discr is a matrix
                    S_discr = @(idx_t, v, vt) S_discr(:, idx_t);
                end
                obj.C_discr = @(idx_t, w) [
                    sparse(obj.m, 1);
                    S_discr(idx_t, w(1:end/2), w(end/2+1:end))
                            ];
            else
                obj.C_discr = [];
            end
            obj.tsImplementation = time.ExplicitRungeKuttaDiscreteData(obj.M, obj.C_cont, obj.C_discr,...
                                                                        k, obj.t, obj.w, order);
        end

        function [v,t,U,T,K] = getV(obj)
            [w,t,U,T,K] = obj.tsImplementation.getV();

            v = w(1:end/2);
            U = U(1:end/2, :); % Stage approximations in previous time step.
            K = K(1:end/2, :); % Stage rates in previous time step.
            % T: Stage times in previous time step.
        end

        function [vt,t,U,T,K] = getVt(obj)
            [w,t,U,T,K] = obj.tsImplementation.getV();

            vt = w(end/2+1:end);
            U = U(end/2+1:end, :); % Stage approximations in previous time step.
            K = K(end/2+1:end, :); % Stage rates in previous time step.
            % T: Stage times in previous time step.
        end

        function [a,b,c,s] = getTableau(obj)
            [a,b,c,s] = obj.tsImplementation.getTableau();
        end

        % Returns quadrature weights for stages in one time step
        function quadWeights = getTimeStepQuadrature(obj)
            [~, b] = obj.getTableau();
            quadWeights = obj.k*b;
        end

        % Use RK for first order form to step
        function obj = step(obj)
            obj.tsImplementation.step();
            [v, t] = obj.tsImplementation.getV();
            obj.w = v;
            obj.t = t;
            obj.n = obj.n + 1;
        end
    end

    methods (Static)
        function k = getTimeStep(lambda, order)
            default_arg('order', 4);
            k = obj.tsImplementation.getTimeStep(lambda, order);
        end
    end

end