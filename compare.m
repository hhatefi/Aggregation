function compare(c1, c2, c3, x1_0, x2_0, interval, PLOT, FILEID, THRESHOLD)

% Validate input arguments
error(nargchk(6, 9, nargin));

% Default Parameters
if nargin < 7 || isempty(PLOT)
    PLOT = 1;
end

if nargin < 8 || isempty(FILEID)
    FILEID = 0;
end

if nargin < 9 || isempty(THRESHOLD),
    THRESHOLD = 0.05;
end

% % %Test
% [time_samples, transient_prob]= solve_splitted_ode(c1, c2, x1_0, x2_0, interval, @ode15s);
% return
% % %Test

% Exact method
[Q_exact, nStates] = gen_mat_exact(c1, c2, c3, x1_0, x2_0);
pi0 = zeros(1, nStates);
pi0(1) = 1;
[t_exact, tr_prob] = tr_analysis(Q_exact, pi0, interval,@ode15s);
[exp_of_x4_exact, trans_prob_exact] = calc_exp_exact(x1_0, x2_0, tr_prob);

% Kinetics ode solution
[t_kinetics_ode, exp_of_x4_kinetics_ode] = solve_kinetics_ode(c1, c2, c3, x1_0, x2_0, interval, @ode45);

% Aggregation method
[Q_aggregated, nStates] = gen_mat_aggregated(c1, c2, c3, x1_0, x2_0);
pi0 = zeros(1, nStates);
pi0(1) = 1;
[t_aggregated, tr_prob] = tr_analysis(Q_aggregated, pi0, interval, @ode45);
exp_of_x4_aggregated = calc_exp_from_prob(x2_0, tr_prob);
[t_error, mse, mre] = calc_error_between_dists(t_aggregated, tr_prob, t_exact, trans_prob_exact, THRESHOLD);
mean_rmse = sum(mse.^0.5) / size(t_error,1);
sd_rmse = sqrt(sum(mse) / size(t_error,1) - mean_rmse^2);
mean_mre = sum(mre) / size(t_error,1);
sd_mre = sqrt(sum(mre .* mre) / size(t_error,1) - mean_mre^2);
if PLOT > 0,
    figure;
    plot(t_error, sqrt(mse), t_error, mre);
    title('Error of distribution for aggregated CTMC');
    legend('RMSE', 'MRE');
    xlabel(sprintf('RMSE: Mean = %0.5f,\tSD = %0.5f\nMRE: Mean = %0.5f,\tSD = %0.5f', mean_rmse, sd_rmse, mean_mre, sd_mre));
end
if FILEID > 0,
    fprintf(FILEID, '%f,%f,%f,%f\n', mean_rmse, sd_rmse, mean_mre, sd_mre);
end

% NAA
exp_of_x4_NAA = calc_exp_NAA(Q_aggregated, t_exact);

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
[t_error, mse, mre] = calc_error_between_dists(t_aggregated, tr_prob, t_exact, trans_prob_exact, THRESHOLD);
mean_rmse = sum(mse.^0.5) / size(t_error,1);
sd_rmse = sqrt(sum(mse) / size(t_error,1) - mean_rmse^2);
mean_mre = sum(mre) / size(t_error,1);
sd_mre = sqrt(sum(mre .* mre) / size(t_error,1) - mean_mre^2);
if PLOT > 0,
    figure;
    plot(t_error, sqrt(mse), t_error, mre);
    title('Error of distribution for MMA, Rate = c3.x1(0)');
    legend('RMSE', 'MRE');
    xlabel(sprintf('RMSE: Mean = %0.5f,\tSD = %0.5f\nMRE: Mean = %0.5f,\tSD = %0.5f', mean_rmse, sd_rmse, mean_mre, sd_mre));
end
if FILEID > 0,
    fprintf(FILEID, '%f,%f,%f,%f\n', mean_rmse, sd_rmse, mean_mre, sd_mre);
end



% MMA (2)
[t_MMA_2, tr_prob] = tr_analysis(Q_MMA{2}, pi0, interval, @ode45);
exp_of_x4_MMA_2 = calc_exp_from_prob(x2_0, tr_prob);
[t_error, mse, mre] = calc_error_between_dists(t_MMA_2, tr_prob, t_exact, trans_prob_exact, THRESHOLD);
mean_rmse = sum(mse.^0.5) / size(t_error,1);
sd_rmse = sqrt(sum(mse) / size(t_error,1) - mean_rmse^2);
mean_mre = sum(mre) / size(t_error,1);
sd_mre = sqrt(sum(mre .* mre) / size(t_error,1) - mean_mre^2);
if PLOT > 0,
    figure;
    plot(t_error, sqrt(mse), t_error, mre);
    title('Error of distribution for MMA, Rate = c3.x1(0).x2(t)/km');
    legend('RMSE', 'MRE');
    xlabel(sprintf('RMSE: Mean = %0.5f,\tSD = %0.5f\nMRE: Mean = %0.5f,\tSD = %0.5f', mean_rmse, sd_rmse, mean_mre, sd_mre));
end
if FILEID > 0,
    fprintf(FILEID, '%f,%f,%f,%f\n', mean_rmse, sd_rmse, mean_mre, sd_mre);
end


% MMA (3)
[t_MMA_3, tr_prob] = tr_analysis(Q_MMA{3}, pi0, interval, @ode45);
exp_of_x4_MMA_3 = calc_exp_from_prob(x2_0, tr_prob);
[t_error, mse, mre] = calc_error_between_dists(t_MMA_3, tr_prob, t_exact, trans_prob_exact, THRESHOLD);
mean_rmse = sum(mse.^0.5) / size(t_error,1);
sd_rmse = sqrt(sum(mse) / size(t_error,1) - mean_rmse^2);
mean_mre = sum(mre) / size(t_error,1);
sd_mre = sqrt(sum(mre .* mre) / size(t_error,1) - mean_mre^2);
if PLOT > 0,
    figure;
    plot(t_error, sqrt(mse), t_error, mre);
    title('Error of distribution for MMA, Rate = 0.5.c3.x1(0)');
    legend('RMSE', 'MRE');
    xlabel(sprintf('RMSE: Mean = %0.5f,\tSD = %0.5f\nMRE: Mean = %0.5f,\tSD = %0.5f', mean_rmse, sd_rmse, mean_mre, sd_mre));
end
if FILEID > 0,
    fprintf(FILEID, '%f,%f,%f,%f\n', mean_rmse, sd_rmse, mean_mre, sd_mre);
end


% MMA (4)
[t_MMA_4, tr_prob] = tr_analysis(Q_MMA{4}, pi0, interval, @ode45);
exp_of_x4_MMA_4 = calc_exp_from_prob(x2_0, tr_prob);
[t_error, mse, mre] = calc_error_between_dists(t_MMA_4, tr_prob, t_exact, trans_prob_exact, THRESHOLD);
mean_rmse = sum(mse.^0.5) / size(t_error,1);
sd_rmse = sqrt(sum(mse) / size(t_error,1) - mean_rmse^2);
mean_mre = sum(mre) / size(t_error,1);
sd_mre = sqrt(sum(mre .* mre) / size(t_error,1) - mean_mre^2);
if PLOT > 0,
    figure;
    plot(t_error, sqrt(mse), t_error, mre);
    title('Error of distribution for MMA with combined rate');
    legend('RMSE', 'MRE');
    xlabel(sprintf('RMSE: Mean = %0.5f,\tSD = %0.5f\nMRE: Mean = %0.5f,\tSD = %0.5f', mean_rmse, sd_rmse, mean_mre, sd_mre));
end
if FILEID > 0,
    fprintf(FILEID, '%f,%f,%f,%f\n', mean_rmse, sd_rmse, mean_mre, sd_mre);
end


% MMA (5)
[t_MMA_5, tr_prob] = tr_analysis(Q_MMA{5}, pi0, interval, @ode45);
exp_of_x4_MMA_5 = calc_exp_from_prob(x2_0, tr_prob);
[t_error, mse, mre] = calc_error_between_dists(t_MMA_5, tr_prob, t_exact, trans_prob_exact, THRESHOLD);
mean_rmse = sum(mse.^0.5) / size(t_error,1);
sd_rmse = sqrt(sum(mse) / size(t_error,1) - mean_rmse^2);
mean_mre = sum(mre) / size(t_error,1);
sd_mre = sqrt(sum(mre .* mre) / size(t_error,1) - mean_mre^2);
if PLOT > 0,
    figure;
    plot(t_error, sqrt(mse), t_error, mre);
    title('Error of distribution for MMA, Rate = c3.x1(0).x2(t)/(x2(t) + km)');
    legend('RMSE', 'MRE');
    xlabel(sprintf('RMSE: Mean = %0.5f,\tSD = %0.5f\nMRE: Mean = %0.5f,\tSD = %0.5f', mean_rmse, sd_rmse, mean_mre, sd_mre));
end
if FILEID > 0,
    fprintf(FILEID, '%f,%f,%f,%f\n', mean_rmse, sd_rmse, mean_mre, sd_mre);
end


% Combining kinetics ode with transient analysis
[t_ckt, tr_prob] = solve_kinetics_and_transient_ode(c1, c2, c3, x1_0, x2_0, interval, @ode45);
exp_of_x4_ckt = calc_exp_from_prob(x2_0, tr_prob);
[t_error, mse, mre] = calc_error_between_dists(t_ckt, tr_prob, t_exact, trans_prob_exact, THRESHOLD);
mean_rmse = sum(mse.^0.5) / size(t_error,1);
sd_rmse = sqrt(sum(mse) / size(t_error,1) - mean_rmse^2);
mean_mre = sum(mre) / size(t_error,1);
sd_mre = sqrt(sum(mre .* mre) / size(t_error,1) - mean_mre^2);
if PLOT > 0,
    figure;
    plot(t_error, sqrt(mse), t_error, mre);
    title('Error of Combining kinetics ode with transient analysis, Rate = c3.x3(t)');
    legend('RMSE', 'MRE');
    xlabel(sprintf('RMSE: Mean = %0.5f,\tSD = %0.5f\nMRE: Mean = %0.5f,\tSD = %0.5f', mean_rmse, sd_rmse, mean_mre, sd_mre));
end
if FILEID > 0,
    fprintf(FILEID, '%f,%f,%f,%f\n', mean_rmse, sd_rmse, mean_mre, sd_mre);
end



% analysing splitted Markov chain
[t_smc, tr_prob]= solve_splitted_ode(c1, c2, c3, x1_0, x2_0, interval, @ode15s);
exp_of_x4_smc = calc_exp_from_prob(x2_0, tr_prob);
[t_error, mse, mre] = calc_error_between_dists(t_smc, tr_prob, t_exact, trans_prob_exact, THRESHOLD);
mean_rmse = sum(mse.^0.5) / size(t_error,1);
sd_rmse = sqrt(sum(mse) / size(t_error,1) - mean_rmse^2);
mean_mre = sum(mre) / size(t_error,1);
sd_mre = sqrt(sum(mre .* mre) / size(t_error,1) - mean_mre^2);
if PLOT > 0,
    figure;
    plot(t_error, sqrt(mse), t_error, mre);
    title('Error of splitted CTMC');
    legend('RMSE', 'MRE');
    xlabel(sprintf('RMSE: Mean = %0.5f,\tSD = %0.5f\nMRE: Mean = %0.5f,\tSD = %0.5f', mean_rmse, sd_rmse, mean_mre, sd_mre));
end
if FILEID > 0,
    fprintf(FILEID, '%f,%f,%f,%f\n', mean_rmse, sd_rmse, mean_mre, sd_mre);
end


if PLOT > 0,
    figure;
    hold on 
    plot(t_exact, exp_of_x4_exact, t_kinetics_ode, exp_of_x4_kinetics_ode, '-.', 'LineWidth', 2);
    plot( t_aggregated, exp_of_x4_aggregated,  t_MMA, exp_of_x4_MMA, t_exact, exp_of_x4_NAA,...
        t_MMA_1, exp_of_x4_MMA_1,t_MMA_2, exp_of_x4_MMA_2,t_MMA_3, exp_of_x4_MMA_3,t_MMA_4, exp_of_x4_MMA_4,t_MMA_5, exp_of_x4_MMA_5, 'LineWidth', 1);
    legend('Exact', 'Kinetics ODE', 'Aggregated', 'MMA', 'NAA', 'MMA (S(t)>>Km)', 'MMA (S(t)<<Km)', 'MMA (S(t)=Km)', 'MMA CTMC', 'MMA Combination')
    title('Expectation of X4(t)');
    hold off
end

fprintf( 'MSE of E[X4(t)] for: \nKinetics ODE solution: %f\nAggregated CTMC (Exact Solution): %f\nNAA: %f\nMMA (Exact Solution): %f\nMMA (S(t)>>Km): %f\nMMA (S(t)<<Km): %f\nMMA (S(t)=Km): %f\nMMA CTMC: %f\nMMA Combination: %f\n',...
calc_RMSE(t_exact, exp_of_x4_exact, t_kinetics_ode, exp_of_x4_kinetics_ode),...
calc_RMSE(t_exact, exp_of_x4_exact, t_aggregated, exp_of_x4_aggregated),...
calc_RMSE(t_exact, exp_of_x4_exact, t_exact, exp_of_x4_NAA),...
calc_RMSE(t_exact, exp_of_x4_exact, t_MMA, exp_of_x4_MMA),...
calc_RMSE(t_exact, exp_of_x4_exact, t_MMA_1, exp_of_x4_MMA_1),...
calc_RMSE(t_exact, exp_of_x4_exact, t_MMA_2, exp_of_x4_MMA_2),...
calc_RMSE(t_exact, exp_of_x4_exact, t_MMA_3, exp_of_x4_MMA_3),...
calc_RMSE(t_exact, exp_of_x4_exact, t_MMA_4, exp_of_x4_MMA_4),...
calc_RMSE(t_exact, exp_of_x4_exact, t_MMA_5, exp_of_x4_MMA_5)...
)

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

function exp_of_x4 = calc_exp_NAA(Q_aggregated, t)
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
