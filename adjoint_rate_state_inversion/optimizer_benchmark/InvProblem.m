% Counted, cached wrapper around AntiplaneShear2DRSFrictionOptHessian for optimizer benchmarking.
% f(x): forward only. fg(x): forward (cached) + adjoint. Hv(x,V): reuses first-order, 2nd-order fwd+adj per column.
classdef InvProblem < handle
    properties
        obj
        xtrue
        nFwd = 0
        nAdj = 0
        nHv = 0
        t0
        maxTime = inf
        maxCost = inf
        lastFwdX = []
        lastF
        lastGradX = []
        lastG
        hist = struct('t',{[]},'f',{[]},'err',{[]},'nFwd',{[]},'nAdj',{[]},'nHv',{[]},'gnorm',{[]},'x',{[]},'info',{{}})
    end
    methods
        function P = InvProblem(obj, xtrue)
            P.obj = obj;
            P.xtrue = xtrue;
            P.t0 = tic;
        end
        function setX(P, x)
            o = P.obj;
            o.setParameters(o.parametersVectorToNative(x));
            o.updateForwardDiscr();
        end
        function f = f(P, x)
            if any(x <= 0)
                f = inf; return
            end
            if isequal(x, P.lastFwdX)
                f = P.lastF; return
            end
            P.setX(x);
            P.lastGradX = []; % adjoint state no longer matches the object's forward state
            P.obj.runForward(false);
            P.nFwd = P.nFwd + 1;
            f = P.obj.computeMisfit();
            P.lastFwdX = x; P.lastF = f;
        end
        function [f, g] = fg(P, x)
            if any(x <= 0)
                f = inf; g = nan(size(x)); return
            end
            f = P.f(x);
            if isequal(x, P.lastGradX)
                g = P.lastG; return
            end
            o = P.obj;
            o.updateAdjointDiscr();
            o.runAdjoint(false);
            P.nAdj = P.nAdj + 1;
            o.firstOrderParsVec = o.parametersNativeToVector(o.forwardDiscr.friction.rsParams);
            g = o.gradientNativeToVector(o.gradientFormula());
            P.lastGradX = x; P.lastG = g;
        end
        function HV = Hv(P, x, V)
            P.fg(x); % make sure first-order solutions are valid at x
            o = P.obj;
            HV = zeros(size(V));
            for k = 1:size(V, 2)
                o.setPerturbationDirection(o.parametersVectorToNative(V(:,k)));
                HV(:,k) = o.gradientNativeToVector(o.computeHessianVector());
                P.nHv = P.nHv + 1;
            end
        end
        % Equivalent gradient evaluations, using measured solve costs (fwd 4.3s, adj 1.3s, Hv 4.1s).
        function c = cost(P)
            c = (4.3*P.nFwd + 1.3*P.nAdj + 4.1*P.nHv)/5.6;
        end
        function log(P, x, f, g, info)
            default_arg('info', '');
            h = P.hist;
            h.t(end+1) = toc(P.t0);
            h.f(end+1) = f;
            h.err(end+1) = norm(x - P.xtrue)/norm(P.xtrue);
            h.nFwd(end+1) = P.nFwd; h.nAdj(end+1) = P.nAdj; h.nHv(end+1) = P.nHv;
            h.gnorm(end+1) = norm(g);
            h.x(:,end+1) = x;
            h.info{end+1} = info;
            P.hist = h;
            fprintf('[%6.0fs] it %3d  f %10.4e  |g| %9.3e  err %7.4f  fwd %4d adj %4d Hv %4d  cost %6.1f  %s\n', ...
                h.t(end), numel(h.f)-1, f, norm(g), h.err(end), P.nFwd, P.nAdj, P.nHv, P.cost(), info);
        end
        function out = timeUp(P)
            out = toc(P.t0) > P.maxTime || P.cost() > P.maxCost;
        end
    end
end
