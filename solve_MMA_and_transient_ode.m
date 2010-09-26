function [time_samples, transient_prob, exp_of_x4] = solve_MMA_and_transient_ode(c1, c2, c3, x1_0, x2_0, interval, solver)

km = (c2 + c3) / c1;

% x(t) = [pi(0); ... p(x2_0); E(t); S(t); C(t); P(t)]
x0 = zeros(x2_0 + 5, 1);
% pi(1) = 1
x0(1) = 1;
% initial number of mulecules
x0(x2_0 + 2:x2_0 + 5) = [x1_0; x2_0; 0; 0];
[time_samples, solution] = solver(@deriva, interval, x0);
transient_prob = solution(:, 1:x2_0 + 1);
exp_of_x4 = solution(:, x2_0 + 5);

    function xdot = deriva(t, x)

    xdot = zeros(x2_0 + 5, 1);

    % Michaelis-Menten-Approximation
    xdot(x2_0 + 3) = -c1 * x1_0 * x(x2_0 + 3) + x1_0 * x(x2_0 + 3) * (c1 * x(x2_0 + 3) + c2) / ( x(x2_0 + 3) + km ) ;
    xdot(x2_0 + 5) = c3 * x1_0 * x(x2_0 + 3) / ( x(x2_0 + 3) + km );


    % Markov chain ode
%     rate = xdot(x2_0 + 5) * (1 - x(x2_0 + 1)) + xdot(x2_0 + 5) * x(x2_0 + 1) / (1.02 - x(x2_0 + 1));

    rate = xdot(x2_0 + 5);
    xdot(1) = - rate * x(1);
    for i = 2:x2_0,
        xdot(i) = rate * x(i - 1) - rate * x(i);
    end
    xdot(x2_0 + 1) = rate * x(x2_0);

    end

end