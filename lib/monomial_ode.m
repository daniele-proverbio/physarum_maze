% monomial_ode.m
% Author: David Palma, PhD, DE <david.palma@uniud.it>
function [tm, v] = monomial_ode(B, d, tspan)
    v0 = double(zeros(size(B, 1), 1));
    epsilon = (1/2 - rand(size(B' * v0))) .* 0.4;
    coeff = 20;
    gamma = 0.5 + epsilon;

    whos B

    [tm, v] = ode15s(@v_dot, tspan, v0);

    function dvdt = v_dot(t, v)
        L = (B * B');
        dvdt = L \ (-B * (((B' * v) ./ gamma).^(2 * coeff + 1))) + d;
    end
end
