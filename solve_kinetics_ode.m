function [t, exp_of_x4] = solve_kinetics_ode(c1, c2, c3, x1_0, x2_0, interval, solver)

% x(t) = [E(t); S(t); C(t); P(t)]
x0 = [x1_0; x2_0; 0; 0];
[t, cons] = solver(@mma_deriva, interval, x0);

exp_of_x4 = cons(:, 4);

    function xdot = mma_deriva(t, x)

    xdot = zeros(4, 1);

    xdot(1) = -c1 * x(1) * x(2) + (c2 + c3)*x(3);
    xdot(2) = -c1 * x(1) * x(2) + c2 * x(3);
    xdot(3) = -xdot(1);
    xdot(4) = c3*x(3);

    end

end