classdef D1StaggeredUpwind < sbp.OpSet
    % Compatible staggered and upwind operators by Ken Mattsson and Ossian O'reilly
    properties
        % x_primal: "primal" grid with m points. Equidistant. Called Plus grid in Ossian's paper.
        % x_dual: "dual" grid with m+1 points. Called Minus grid in Ossian's paper.

        % D1_primal takes FROM dual grid TO primal grid
        % D1_dual takes FROM primal grid TO dual grid

        D1_primal % SBP operator approximating first derivative
        D1_dual % SBP operator approximating first derivative

        Dplus_primal  % Upwind operator on primal grid
        Dminus_primal % Upwind operator on primal grid
        Dplus_dual % Upwind operator on dual grid
        Dminus_dual % Upwind operator on dual grid

        H_primal % Norm matrix
        H_dual % Norm matrix
        H_primalI % H^-1
        H_dualI % H^-1
        e_primal_l % Left boundary operator
        e_dual_l % Left boundary operator
        e_primal_r % Right boundary operator
        e_dual_r % Right boundary operator
        m % Number of grid points.
        m_primal % Number of grid points.
        m_dual % Number of grid points.
        h % Step size
        x_primal % grid
        x_dual % grid
        x
        borrowing % Struct with borrowing limits for different norm matrices
    end

    methods
        function obj = D1StaggeredUpwind(m,lim,order)

          xl = lim{1};
          xr = lim{2};
          L = xr-xl;
          h = L/(m-1);

          m_primal = m;
          m_dual = m+1;

          switch order
          case 2
            [~, ~, obj.H_primal, obj.H_dual,...
            obj.H_primalI, obj.H_dualI,...
            obj.D1_primal, obj.D1_dual, obj.Dplus_primal, obj.Dminus_primal,...
            obj.Dplus_dual, obj.Dminus_dual] = sbp.implementations.d1_staggered_upwind_2(m, L);
          case 4
            [~, ~, obj.H_primal, obj.H_dual,...
            obj.H_primalI, obj.H_dualI,...
            obj.D1_primal, obj.D1_dual, obj.Dplus_primal, obj.Dminus_primal,...
            obj.Dplus_dual, obj.Dminus_dual] = sbp.implementations.d1_staggered_upwind_4(m, L);
          case 6
            [~, ~, obj.H_primal, obj.H_dual,...
            obj.H_primalI, obj.H_dualI,...
            obj.D1_primal, obj.D1_dual, obj.Dplus_primal, obj.Dminus_primal,...
            obj.Dplus_dual, obj.Dminus_dual] = sbp.implementations.d1_staggered_upwind_6(m, L);
          case 8
            [~, ~, obj.H_primal, obj.H_dual,...
            obj.H_primalI, obj.H_dualI,...
            obj.D1_primal, obj.D1_dual, obj.Dplus_primal, obj.Dminus_primal,...
            obj.Dplus_dual, obj.Dminus_dual] = sbp.implementations.d1_staggered_upwind_8(m, L);
          otherwise
           error('Invalid operator order %d.',order);
          end

          obj.m = m;
          obj.m_primal = m_primal;
          obj.m_dual = m_dual;
          obj.h = h;

          obj.x_primal = linspace(xl, xr, m)';
          obj.x_dual = [xl, linspace(xl+h/2, xr-h/2, m-1), xr]';

          obj.e_primal_l = sparse(m_primal,1);
          obj.e_primal_r = sparse(m_primal,1);
          obj.e_primal_l(1) = 1;
          obj.e_primal_r(m_primal) = 1;

          obj.e_dual_l = sparse(m_dual,1);
          obj.e_dual_r = sparse(m_dual,1);
          obj.e_dual_l(1) = 1;
          obj.e_dual_r(m_dual) = 1;

          obj.borrowing = [];
          obj.x = [];

        end

        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end
    end
end
