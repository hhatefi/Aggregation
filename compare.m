function compare(c1, c2, c3, x1_0, x2_0, interval)
% Exact method
[Q_exact, nStates] = gen_mat_exact(c1, c2, c3, x1_0, x2_0);
pi0 = zeros(1, nStates);
pi0(1) = 1;
[t_exact, tr_prob] = tr_analysis(Q_exact, pi0, interval,@ode15s);
[exp_of_x4_exact, trans_prob_exact] = calc_exp_exact(x1_0, x2_0, tr_prob);

% Aggregation method
[Q_aggregated, nStates] = gen_mat_aggregated(c1, c2, c3, x1_0, x2_0);
pi0 = zeros(1, nStates);
pi0(1) = 1;
[t_aggregated, tr_prob] = tr_analysis(Q_aggregated, pi0, interval, @ode45);
exp_of_x4_aggregated = calc_exp_from_prob(x2_0, tr_prob);
[t_mse, mse] = calc_error_between_dists(t_aggregated, tr_prob, t_exact, trans_prob_exact);
figure;
plot(t_mse, sqrt(mse));
title('MSE of distribution for aggregated CTMC');
mean = sum(mse.^0.5) / size(t_mse,1);
sd = sqrt(sum(mse) / size(t_mse,1) - mean^2);
xlabel(sprintf('Mean = %0.3f,\tSD = %0.3f', mean, sd));

% NAA
exp_of_x4_NAA = calc_exp_aggregated_app(Q_aggregated, t_exact);

% Michaelis-Menten Approximation (Exact solution of Differetial Eq.)
[t_MMA, exp_of_x4_MMA] = calc_exp_MMA(c1, c2, c3, x1_0, x2_0, interval);

% Michaelis-Menten Approximation for: 
% 1) Rate = c3.x1_0 (in case x2(t) >> Km),
% 2) Rate = c3.x1_0.x2(t)/km (in case x2(t) << Km),
% 3) Rate = 0.5.c3.x1_0 (in case x2(t) = Km),
% 4) Rate = c3.x1_0.x2(t)/(x2(t) + km)
% 5) Rate is combination of (1), (2) and (3)

[Q_MMA, nStates] = gen_mat_for_MMA(c1, c2, c3, x1_0, x2_0);
pi0 = zeros(1, nStates);
pi0(1) = 1;
% MMA (1)
[t_MMA_1, tr_prob] = tr_analysis(Q_MMA{1}, pi0, interval, @ode45);
exp_of_x4_MMA_1 = calc_exp_from_prob(x2_0, tr_prob);
[t_mse, mse] = calc_error_between_dists(t_MMA_1, tr_prob, t_exact, trans_prob_exact);
figure;
plot(t_mse, sqrt(mse));
title('MSE of distribution for MMA, Rate = c3.x1(0)');
mean = sum(mse.^0.5) / size(t_mse,1);
sd = sqrt(sum(mse) / size(t_mse,1) - mean^2);
xlabel(sprintf('Mean = %0.3f,\tSD = %0.3f', mean, sd));

% MMA (2)
[t_MMA_2, tr_prob] = tr_analysis(Q_MMA{2}, pi0, interval, @ode45);
exp_of_x4_MMA_2 = calc_exp_from_prob(x2_0, tr_prob);
[t_mse, mse] = calc_error_between_dists(t_MMA_2, tr_prob, t_exact, trans_prob_exact);
figure;
plot(t_mse, sqrt(mse));
title('MSE of distribution for MMA, Rate = c3.x1(0).x2(t)/km');
mean = sum(mse.^0.5) / size(t_mse,1);
sd = sqrt(sum(mse) / size(t_mse,1) - mean^2);
xlabel(sprintf('Mean = %0.3f,\tSD = %0.3f', mean, sd));

% MMA (3)
[t_MMA_3, tr_prob] = tr_analysis(Q_MMA{3}, pi0, interval, @ode45);
exp_of_x4_MMA_3 = calc_exp_from_prob(x2_0, tr_prob);
[t_mse, mse] = calc_error_between_dists(t_MMA_3, tr_prob, t_exact, trans_prob_exact);
figure;
plot(t_mse, sqrt(mse));
title('MSE of distribution for MMA, Rate = 0.5.c3.x1(0)');
mean = sum(mse.^0.5) / size(t_mse,1);
sd = sqrt(sum(mse) / size(t_mse,1) - mean^2);
xlabel(sprintf('Mean = %0.3f,\tSD = %0.3f', mean, sd));

% MMA (4)
[t_MMA_4, tr_prob] = tr_analysis(Q_MMA{4}, pi0, interval, @ode45);
exp_of_x4_MMA_4 = calc_exp_from_prob(x2_0, tr_prob);
[t_mse, mse] = calc_error_between_dists(t_MMA_3, tr_prob, t_exact, trans_prob_exact);
figure;
plot(t_mse, sqrt(mse));
title('MSE of distribution for MMA with combined rate');
mean = sum(mse.^0.5) / size(t_mse,1);
sd = sqrt(sum(mse) / size(t_mse,1) - mean^2);
xlabel(sprintf('Mean = %0.3f,\tSD = %0.3f', mean, sd));

% MMA (5)
[t_MMA_5, tr_prob] = tr_analysis(Q_MMA{5}, pi0, interval, @ode45);
exp_of_x4_MMA_5 = calc_exp_from_prob(x2_0, tr_prob);
[t_mse, mse] = calc_error_between_dists(t_MMA_3, tr_prob, t_exact, trans_prob_exact);
figure;
plot(t_mse, sqrt(mse));
title('MSE of distribution for MMA, Rate = c3.x1(0).x2(t)/(x2(t) + km)');
mean = sum(mse.^0.5) / size(t_mse,1);
sd = sqrt(sum(mse) / size(t_mse,1) - mean^2);
xlabel(sprintf('Mean = %0.3f,\tSD = %0.3f', mean, sd));

figure;
hold on 
plot(t_exact, exp_of_x4_exact, 'LineWidth', 2);
plot(t_aggregated, exp_of_x4_aggregated, t_MMA, exp_of_x4_MMA, t_exact, exp_of_x4_NAA,...
    t_MMA_1, exp_of_x4_MMA_1,t_MMA_2, exp_of_x4_MMA_2,t_MMA_3, exp_of_x4_MMA_3,t_MMA_4, exp_of_x4_MMA_4,t_MMA_5, exp_of_x4_MMA_5, 'LineWidth', 1);
legend('Exact', 'Aggregated', 'MMA', 'NAA', 'MMA (S(t)>>Km)', 'MMA (S(t)<<Km)', 'MMA (S(t)=Km)', 'MMA CTMC', 'MMA Combination')
title('Expectation of X4(t)');
hold off

fprintf( 'MSE for: \nAggregated CTMC (Exact Solution): %f\nNAA: %f\nMMA (Exact Solution): %f\nMMA (S(t)>>Km): %f\nMMA (S(t)<<Km): %f\nMMA (S(t)=Km): %f\nMMA CTMC: %f\nMMA Combination: %f\n',...
calc_MSE(t_exact, exp_of_x4_exact, t_aggregated, exp_of_x4_aggregated),...
calc_MSE(t_exact, exp_of_x4_exact, t_exact, exp_of_x4_NAA),...
calc_MSE(t_exact, exp_of_x4_exact, t_MMA, exp_of_x4_MMA),...
calc_MSE(t_exact, exp_of_x4_exact, t_MMA_1, exp_of_x4_MMA_1),...
calc_MSE(t_exact, exp_of_x4_exact, t_MMA_2, exp_of_x4_MMA_2),...
calc_MSE(t_exact, exp_of_x4_exact, t_MMA_3, exp_of_x4_MMA_3),...
calc_MSE(t_exact, exp_of_x4_exact, t_MMA_4, exp_of_x4_MMA_4),...
calc_MSE(t_exact, exp_of_x4_exact, t_MMA_5, exp_of_x4_MMA_5)...
)

return;

function [exp_of_x4, trans_prob] = calc_exp_exact(x1_0, x2_0, tr_prob)

nSamples = size(tr_prob, 1);
nStates = size(tr_prob, 2);
exp_of_x4 = zeros(1, nSamples);
trans_prob = zeros(nSamples, x2_0 + 1);
for time_index = 1:nSamples,
	state_index = 1;
	% for each column: there are x2_0 + 1 columns and each state in column i has x(4) = i
	for i = 0:x2_0,
		x = [x1_0 x2_0 - i 0 i];
		while x(1) >= 0 && x(2) >= 0,
			exp_of_x4(time_index) = exp_of_x4(time_index) + x(4) * tr_prob(time_index, state_index);
            trans_prob(time_index, x(4) + 1) = trans_prob(time_index, x(4) + 1) + tr_prob(time_index, state_index);
			state_index = state_index + 1;
			x(1) = x(1) - 1;
			x(2) = x(2) - 1;
			x(3) = x(3) + 1;
        end
    end
end

return;

function exp_of_x4 = calc_exp_from_prob(x2_0, tr_prob)

nSamples = size(tr_prob, 1);
nStates = size(tr_prob, 2);
exp_of_x4 = zeros(1, nSamples);
for time_index = 1:nSamples,
	% for each column: there are x2_0 + 1 columns and each state in column i has x(4) = i
	for i = 0:x2_0,
		exp_of_x4(time_index) = exp_of_x4(time_index) + i * tr_prob(time_index, i + 1);
	end
end

return;

function exp_of_x4 = calc_exp_aggregated_app(Q_aggregated, t)
% t is column vector 
nSamples = size(t, 1);
nStates = size(Q_aggregated, 1);
exp_of_x4 = zeros(1, nSamples);

left = 0;
right = 1 / Q_aggregated(1, 2);
j = 1;
for i = 2:nStates - 1,
	while j <= nSamples && t(j) <= right,
		if t(j) > left && t(j) <= right,
			exp_of_x4(j) = i - 1;
		end
		j = j + 1;
	end
	left = right;
	right = left + 1 / Q_aggregated(i, i + 1);
end

while j <= nSamples,
	exp_of_x4(j) = nStates - 1;
	j = j + 1;
end


function [Q_MMA, nStates] = gen_mat_for_MMA(c1, c2, c3, x1_0, x2_0)

    N = 5;
    % Q_MMA is a cell array 
    Q_MMA = cell(1, N);
    km = (c2 + c3) / c1; 
    nStates =  x2_0 + 1;
    Rows = zeros(1, 2 * x2_0 );
    Cols = zeros(N, 2 * x2_0);
    NZ = zeros(N, 2 * x2_0);
    
    for i = 1: x2_0,
        Rows(2*i - 1) = i;
        Rows(2*i) = i;

        for j = 1:N,
            Cols(j, 2*i - 1) = i;
            Cols(j, 2*i) = i + 1;
        end
        
        % Adding non-zero elements
        NZ(1, 2*i - 1) = - c3 * x1_0;
        NZ(2, 2*i - 1) = - c3 * x1_0 * (x2_0 - i + 1) / km;
        NZ(3, 2*i - 1) = - 0.5 * c3 * x1_0;
        NZ(4, 2*i - 1) = - c3 * x1_0 * (x2_0 - i + 1) / (km + x2_0 - i + 1);

        ratio = (x2_0 - i + 1) / (km + x2_0 - i + 1);
        [min_error, min_error_index] = min(abs([ratio - 1, ratio - (x2_0 - i + 1) / km, ratio - 0.5]));
        NZ(5, 2*i - 1) = NZ(min_error_index, 2*i - 1);
        
        for j = 1:N,
            NZ(j, 2*i) = - NZ(j, 2*i - 1);
        end
    end
    
    for i = 1:N,
        Q_MMA{i} = sparse(Rows, Cols(i,:), NZ(i,:), nStates, nStates);
    end


    



