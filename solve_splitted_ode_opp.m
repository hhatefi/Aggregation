function [time_samples, transient_prob]= solve_splitted_ode_opp(c1, c2, c3, x1_0, x2_0, interval, solver)

% Q_exact is the matrix of orignal Markov chain
% Q_splitted_T is the collection of matrices that each of them is the 
%   transpose of generator matrix of a splitted markov chain

Q_splitted_T = cell(x2_0, 1);
time_samples = cell(x2_0, 1);
exp_of_x3s = cell(x2_0, 1);

for i = 0:x2_0 - 1,
	x = [x1_0 x2_0 - i 0 i];
    nState4ithStage = min([x1_0 x2_0 - i]) + 1;
    nTransitions = (nState4ithStage - 1) * 2;

    transition_index = 1;
    state_index = 1;
    I = zeros(1, nState4ithStage + nTransitions);
    J = zeros(1, nState4ithStage + nTransitions);
    NZ = zeros(1, nState4ithStage + nTransitions);
    while x(1) >= 0 && x(2) >= 0,
		exit_rate = 0;
		% R1
		rate = c1 * x(1) * x(2);
		if rate ~= 0,
			I(transition_index) = state_index;
			J(transition_index) = state_index + 1;
			NZ(transition_index) = rate;
			exit_rate = exit_rate + rate;
			transition_index = transition_index + 1;
		end
		% R2
		rate = c2 * x(3);
		if rate ~= 0,
			I(transition_index) = state_index;
			J(transition_index) = state_index - 1;
			NZ(transition_index) = rate;
			exit_rate = exit_rate + rate;
			transition_index = transition_index + 1;
		end
		% Add diagonal element
		if exit_rate ~= 0,
			I(transition_index) = state_index;
			J(transition_index) = state_index; 
			NZ(transition_index) = - exit_rate;
			transition_index = transition_index + 1;
		end
		state_index = state_index + 1;
		x(1) = x(1) - 1;
		x(2) = x(2) - 1;
		x(3) = x(3) + 1;
    end
    Q_splitted_T{i + 1} = sparse(I, J, NZ, nState4ithStage, nState4ithStage);    
    
    % calculate the approximated expectation of x3 
    [time_sample_i, exp_of_x3_i] = app_exp_of_x3_in_birth_death_process(Q_splitted_T{i + 1}, interval);
    time_samples{i + 1} = time_sample_i;
    exp_of_x3s{i + 1} = exp_of_x3_i;

end
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
        lambda_i_t = -1;
        if i < 1 || i > x2_0,
            lambda_i_t = 0;
        else
            len = length(time_samples{i});
            for k = 1:len - 1,
                if t >= time_samples{i}(k) && t <= time_samples{i}(k + 1),
                    exp_of_x3_t = (exp_of_x3s{i}(k + 1) - exp_of_x3s{i}(k)) / (time_samples{i}(k + 1) - time_samples{i}(k)) * (t - time_samples{i}(k)) + exp_of_x3s{i}(k);
                    lambda_i_t = c3 * exp_of_x3_t;
                    break;
                end
            end
            if lambda_i_t == -1,
                fprintf(2, 'ERROR: No interpolation interval found.\n');
            end
        end
    end

end
