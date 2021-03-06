function [time_samples, exp_of_x3] = app_exp_of_x3_in_birth_death_process(Q, interval)

% Precision for detection of steady state
PRECISION = 1E-4;
cur_exp_of_x3 = 0;
cur_time = 0;

i = 1;

while cur_time < interval(2),
    
    time_samples(i) = cur_time;
    exp_of_x3(i) = cur_exp_of_x3;
    
    i = i + 1;

    % compute approximated probabilities and states
    
    % Lower State and its Probability
    lower_state = floor(cur_exp_of_x3);
    lower_state_prob = lower_state + 1 - cur_exp_of_x3;

    % Upper State and its Probability
    upper_state = lower_state + 1;
    upper_state_prob = 1 - lower_state_prob;

    % expected sojourn time in lower state
    ls_exp_sojourn_time = - 1 / Q(lower_state + 1, lower_state + 1);
    
    % expected sojourn time in upper state
    if upper_state <= size(Q, 1) - 1,
        us_exp_sojourn_time = - 1 / Q(upper_state + 1, upper_state + 1);
    else
        if upper_state_prob > 1e-5,
            fprintf(2, 'Error: Probability is %f.\n', upper_state_prob);
        end
        us_exp_sojourn_time = 0;
    end
    
    % Expected Sojourn Time
    exp_sojourn_time = lower_state_prob * ls_exp_sojourn_time + upper_state_prob * us_exp_sojourn_time;
    
    % Expection of x3 for lower state
    if lower_state == 0,
        ls_exp_of_x3 = 1;
    elseif lower_state == size(Q, 1) - 1,
        ls_exp_of_x3 = lower_state - 1;
    else
        prob = 1;% - exp(Q(lower_state + 1, lower_state + 1) * exp_sojourn_time);
        ls_exp_of_x3 = (-(lower_state - 1) * Q(lower_state + 1, lower_state)/Q(lower_state + 1, lower_state + 1) - ...
            (lower_state + 1) * Q(lower_state + 1, lower_state + 2) / Q(lower_state + 1, lower_state + 1)) * prob  + lower_state * (1 - prob);
    end

    % Expection of x3 for upper state
    if upper_state == 0,
        us_exp_of_x3 = 1;
    elseif upper_state == size(Q, 1) - 1,
        us_exp_of_x3 = upper_state - 1;
    elseif upper_state > size(Q, 1) - 1,
        if upper_state_prob > 1e-5,
            fprintf(2, 'Error: Probability is %f.\n', upper_state_prob);
        end
        us_exp_of_x3 = 0;
    else
        prob = 1;% - exp(Q(upper_state + 1, upper_state + 1) * exp_sojourn_time);
        us_exp_of_x3 = (-(upper_state - 1) * Q(upper_state + 1, upper_state)/Q(upper_state + 1, upper_state + 1) - ...
            (upper_state + 1) * Q(upper_state + 1, upper_state + 2) / Q(upper_state + 1, upper_state + 1)) * prob + upper_state * (1 - prob);
    end
            
    cur_time = cur_time + exp_sojourn_time;
    cur_exp_of_x3 = lower_state_prob * ls_exp_of_x3 + upper_state_prob * us_exp_of_x3;
    
    if abs(cur_exp_of_x3 - exp_of_x3(i - 1)) < PRECISION || cur_exp_of_x3 - exp_of_x3(i - 1) < 0,
        time_samples(i) = cur_time;
        exp_of_x3(i) = calc_exp_of_x3_ss();
        
        cur_time = interval(2);
        cur_exp_of_x3 = exp_of_x3(i);
        
        i = i + 1;

    end

end
time_samples(i) = cur_time;
exp_of_x3(i) = cur_exp_of_x3;


    function exp_of_x3_ss = calc_exp_of_x3_ss()
        nStates = size(Q, 1);
        p = zeros(1, nStates);
        p(nStates) = 1;

        for j = (nStates - 1):(-1):1,
            p(j) = Q(j + 1, j) / Q(j, j + 1) * p(j + 1);
        end

        sum_p = 0;
        exp_of_x3_ss = 0;
        for j = 1:nStates,
            sum_p = sum_p + p(j);
            exp_of_x3_ss = exp_of_x3_ss + (j - 1) * p(j);
        end
        exp_of_x3_ss = exp_of_x3_ss / sum_p;

    end
end