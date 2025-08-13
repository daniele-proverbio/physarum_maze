% piecewise_linear.m
% The function implements the dynamic system used to explain why the
% lightning always chooses a specific single path, using a piecewise-linear
% threshold function.
%
% Author: David Palma, PhD, DE <david.palma@uniud.it>
%
% Usage:
%   [v,u,obj] = piecewise_linear(B,d,delta,T,debug,r,c)
%
% Arguments:
%   B:     incidence matrix
%   d:     total current outflowing from the nodes
%   delta: range of the air dielectric rigidity
%   T:     time
%   debug: debug mode
%   r:     no. of rows
%   c:     no. of columns
%
% Returns:
%   v:     vector of the node tensions (1,2,...,n)
%   u:     vector of the arc tensions (1,2,...,m)
%   obj:   sequence of images
function [obj, v, u, epsilon] = piecewise_linear(B, d, delta, T, debug, r, c)

    %% Check inputs
    narginchk(2, 7);

    if nargin < 7
        debug = 0; % debug mode
    elseif nargin < 4
        T = 4; % total time
    elseif nargin < 3
        delta = 0.3; % range of the air dielectric rigidity
    end

    %% Check dimensions/values and initialise variables
    assert(abs(delta) < 1, ...
        'Error: Incorrect value for delta (must be between 0 and 1).');

    assert(T > 0, ...
        'Error: Incorrect value for T (must be positive).');

    assert(size(B, 1) == size(d, 1), ...
        'Error: Incorrect dimensions of B and d.');

    v = double(zeros(size(B, 1), 1)); % node tensions
    u = double(zeros(size(B, 2), 1)); % current
    tau = 1.260e-3; % sampling time (k*tau < 1)
    iters = round(T / tau); % no. of iterations

    % vector used to randomise the air dielectric rigidity
    epsilon = (1/2 - rand(size(B' * v))) .* delta; % epsilon in (-delta/2,+delta/2)

    % pre-compute the inverse of B*B'
    L = B * B';
    L_inv = inv(L);

    if debug
        whos B L d v u
        fprintf('Total time: %g\n', T);
        fprintf('Sampling time: %g\n', tau);
        fprintf('Number of iterations: %d\n', iters);
        fprintf('Range of the air dielectric rigidity: (%g,%g)\n', ...
            -delta / 2, delta / 2);

        % sequence of images
        obj = zeros(r, c, iters);

        fprintf('Computation in progress... %3d%%', 0);
    end

    %% Dynamic algorithm
    for t = 1:iters
        u = sat(B' * v, epsilon);
        v = v + tau * L_inv * ((-B * u) + d); %#ok<MINV>

        % debug mode
        if debug
            cnt = 0;
            tmp = zeros(size(u));
            tmp(u > 0) = u(u > 0);

            for i = 1:r
                for j = 1:c - 1
                    % rows
                    cnt = cnt + 1;
                    obj(i, j:j + 1, t) = tmp(cnt);
                end
                for j = 1:c
                    % columns
                    cnt = cnt + 1;
                    obj(i:i + 1, j, t) = tmp(cnt);
                end
            end
            fprintf('\b\b\b\b%3d%%', round(100 * t / iters));
        end
    end

    % normalisation
    if debug
        obj = obj ./ (max(u(:)));
    end

    if debug
        fprintf('\ndone.\n\n');
    end
end

function u = sat(y, epsilon)
    % Electrostatic threshold function that relates the vector of the arc
    % tensions y and the tension law u. The expression of the tension y on the
    % arcs (on the capacitors) given the node tension v is y=B'v.

    %% Check inputs and initialize variables
    narginchk(1, 2)

    if nargin == 1
        epsilon = 0;
    end

    %% Initialise variables
    kappa = 800; % 800 : close to ideal; 300: farther from ideal
    gamma = 0.5 + epsilon; % heterogeneity
    u = zeros(size(y)); % arc tensions
    epsi = 10e-5; % coefficient

    %% Piecewise-linear threshold function
    u(y > +gamma) = kappa .* (y(y > +gamma) - gamma(y > +gamma)) + ...
        gamma(y > +gamma) .* epsi;
    u(y < -gamma) = kappa .* (y(y < -gamma) + gamma(y < -gamma)) - ...
        gamma(y < -gamma) .* epsi;
    u(y >= -gamma & y <= +gamma) = y(y >= -gamma & y <= +gamma) .* epsi;
end
