% dynamic_system_simulation.m
% Simulation of the time varying behaviour of the dynamic system running
% the mathematical model associated with said system.
%
% Author: David Palma, PhD, DE <david.palma@uniud.it>
%
% Arguments:
%   r:     no. of rows
%   c:     no. of columns
%   mode:  mode -> ground (0), end-to-end (1)
%   delta: range of the air dielectric rigidity
%   T:     time
%   debug: debug mode
%   th:    th -> monomial (0), piecewise-linear (1)
%
% Returns:
%   none
function obj = dynamic_system_simulation(r, c, mode, delta, T, debug, th)

    %% Add "lib" folder and subfolders to path
    addpath(genpath('../lib'));

    %% build the generalized incidence matrix of the graph
    B = build_matrix_B(r, c, mode);

    %% total current outflowing from the nodes
    d = zeros(size(B, 1), 1);
    d(floor(c / 2)) = 8;

    %% simulation
    if th == 1
        [obj, v, u, epsilon] = piecewise_linear(B, d, delta, T, debug, r, c);
    elseif th == 0
        [~, ~, obj] = monomial(B, d, delta, T, debug, r, c);
    elseif th == 2
        [obj, v, u, epsilon] = linear(B, d, delta, T, debug, r, c);
    end

    % save('wspace');
    fprintf('done.\n\n');

end
