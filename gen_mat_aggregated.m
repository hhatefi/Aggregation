function [Q, nStates] = gen_mat_aggregated(c1, c2, c3, x1_0, x2_0)

nStates =  x2_0 + 1;
I = zeros(1, 2 * x2_0 );
J = zeros(1, 2 * x2_0);
NZ = zeros(1, 2 * x2_0);
p = zeros(1, x1_0);
for k = 0:x2_0 - 1,
	m_k = min(x1_0, x2_0 - k);

	p(1) = 1;
	for j = 2:m_k + 1,
		p(j) = p(j - 1) * (x2_0 + 2 - k - j) * (x1_0 + 2 - j) / (c2 * (j - 1));
	end
	sum_p = 0;
	exp_of_x3 = 0;
	for i = 1:m_k + 1,
		sum_p = sum_p + p(i);
		exp_of_x3 = exp_of_x3 + (i - 1) * p(i);
	end
	exp_of_x3 = exp_of_x3 / sum_p;

	I(2 * k + 1) = k + 1;
	J(2 * k + 1) = k + 1;
	NZ(2 * k + 1) = - exp_of_x3 * c3;

	I(2 * k + 2) = k + 1;
	J(2 * k + 2) = k + 2;
	NZ(2 * k + 2) = -NZ(2 * k + 1);

end
Q = sparse(I, J, NZ, nStates, nStates);
return;
