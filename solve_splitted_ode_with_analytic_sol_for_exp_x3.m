function [time_samples, transient_prob]= solve_splitted_ode_with_analytic_sol_for_exp_x3(c1, c2, c3, x1_0, x2_0, interval, solver)

pi0 = zeros(x2_0 + 1, 1);
pi0(1) = 1;
[time_samples, transient_prob] = solver(@deriva_main, interval, pi0);


    function pdot = deriva_main(t, p)
        pdot = zeros(x2_0 + 1, 1);
        
        pdot(1) = - calc_lambda(t, 1) * p(1);
        for j = 2:x2_0,
            pdot(j) = calc_lambda(t, j - 1) * p(j - 1) - calc_lambda(t, j) * p(j);
        end
        pdot(x2_0 + 1) = calc_lambda(t, x2_0) * p(x2_0);
        
    end

    function lambda_i_t = calc_lambda(t, i)
        % calculate outgoing lambda from state i at time t
        lambda_i_t = c3 * calc_x3_inside_macro_state(c1, c2, x1_0, x2_0 - i + 1, t);
    end

end
