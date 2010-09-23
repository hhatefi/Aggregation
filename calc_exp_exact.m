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