% monomial_test.m
% Author: David Palma, PhD, DE <david.palma@uniud.it>
r = 66;
c = 22;
B = build_matrix_B(r, c, 0);
d = zeros(size(B, 1), 1);
d(floor(c / 2)) = 8;

tspan = [0 1];

fprintf('Test #1');
tic
[tm1, v1] = monomial_ode(B, d, tspan);
toc

tspan = [0 4];

fprintf('Test #2');
tic
[tm2, v2] = monomial_ode(B, d, tspan);
toc
