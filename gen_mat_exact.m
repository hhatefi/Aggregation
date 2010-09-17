function [Q, nStates] = gen_mat_exact(c1, c2, c3, x1_0, x2_0)

%nTransitions = 3 * (x2_0 - x1_0) * x1_0 + 3 * x1_0 * (x1_0 + 1) / 2;
%nStates = (x1_0 + 1) * (x2_0 + 1 - x1_0 / 2);

% calculating number of states and transitions
nStates = 0;
nTransitions = 0;
for i = 0:x2_0,
    nState4ithStage = min([x1_0 x2_0 - i]) + 1;
    nStates = nStates + nState4ithStage;
    nTransitions = nTransitions + (nState4ithStage - 1) * 3;
end
% Number of nonzero elements is nTransitions plus diagonal elemnts minus one (one state is absorbing)
I = zeros (1, nTransitions + nStates - 1);
J = zeros (1, nTransitions + nStates - 1);
NZ = zeros(1, nTransitions + nStates - 1);

transition_index = 1;
state_index = 1;
% for each column: there are x2_0 + 1 columns and each state in column i has x(4) = i
for i = 0:x2_0,
	x = [x1_0 x2_0 - i 0 i];
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
		% R3
		rate = c3 * x(3);
		if rate ~= 0,
			I(transition_index) = state_index;
			J(transition_index) = state_index + min(x1_0, x2_0 - i);
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
end
Q = sparse(I, J, NZ, nStates, nStates);



