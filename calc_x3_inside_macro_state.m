function x3 = calc_x3_inside_macro_state(c1, c2, x1_0, x2_0, t)

root_vector = roots([c1 -(c1*x1_0 + c1*x2_0 + c2) c1*x1_0*x2_0]);

a1 = max(root_vector);
a2 = min(root_vector);

alpha = -c1*(a1 - a2)*t;
x3 = a1 * (exp(alpha) - 1) ./ (exp(alpha) - a1/a2);