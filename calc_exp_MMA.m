function [t, exp_of_x4] = calc_exp_MMA(c1, c2, c3, x1_0, x2_0, interval)

km = (c2 + c3) / c1;

x0 = [x2_0; 0];
[t, cons] = ode45(@mma_deriva, interval, x0);

exp_of_x4 = cons(:, 2);

    function xdot = mma_deriva(t, x)

    xdot = zeros(2, 1);

    xdot(1) = -c1 * x1_0 * x(1) + x1_0 * x(1) * (c1 * x(1) + c2) / ( x(1) + km ) ;
    xdot(2) = c3 * x1_0 * x(1) / ( x(1) + km );

    end

end
