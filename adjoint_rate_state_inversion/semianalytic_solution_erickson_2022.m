function [u_m, u_p, V] = semianalytic_solution_erickson_2022(F,x_s,sigma)
    syms x t
    u0(x) = exp(-((x-x_s)/sigma).^2);
    u0_t(x) = subs(diff(u0(x-t),'t'),t,0);
    u0_x(x) = diff(u0,'x');

    u0 = matlabFunction(u0);
    u0_t = matlabFunction(u0_t);
    u0_x = matlabFunction(u0_x);
    clear x t;
   
    w0_m = @(t) u0_t(-t) - u0_x(-t);
    w0_p = @(t) u0_t(t) + u0_x(t);
    eta = @(t) w0_p(t) - w0_m(t);
    g = @(t,V) 2*F(V) + V - eta(t);
    V = @(t) arrayfun(@(t) (t>0)*fzero(@(V) g(t,V), [0, eta(t)]), t);

    Q_m = @(t) w0_m(t) + 2*F(V(t));
    Q_p = @(t) w0_p(t) - 2*F(V(t));

    Psi_m = @(tau) 1/2*(tau>0)*integral(Q_m, 0, tau,'AbsTol',1e-10,'RelTol',1e-10);
    Psi_p = @(tau) 1/2*(tau>0)*integral(Q_p, 0, tau,'AbsTol',1e-10,'RelTol',1e-10);

    u_m = @(x,t) u0(x - t) + arrayfun(@(x) Psi_m(t + x), x);
    u_p = @(x,t) u0(x + t) + arrayfun(@(x) Psi_p(t - x), x);
end