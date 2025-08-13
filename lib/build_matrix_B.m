% build_matrix_B.m
% The function builds the generalized incidence matrix of the graph with n
% rows corresponding to the nodes, and m columns corresponding to the arcs.
%
% Author: David Palma, PhD, DE <david.palma@uniud.it>
%
% Usage:
%   B = build_matrix_B(rows, cols, method)
%
% Arguments:
%   rows:   rows of the graph
%   cols:   columns of the graph
%   method: 0 (ground) or 1 (end-to-end)
%
% Returns:
%   B:      incidence matrix
function B = build_matrix_B(rows, cols, method)

    %% Check inputs and initialize variables
    narginchk(2, 3)

    if nargin < 3
        method = 0;
    end

    nodes = cols * rows;
    arcs = (cols - 1) * rows + cols * (rows - 1);
    B = zeros(nodes, arcs);
    alpha = 1;

    %% Build matrix B
    for r = 0:rows - 1
        % compute the indices
        idx = r * (cols - 1) + cols * r;
        i = r * cols + 1;
        j = i + cols - 1;

        % horizontal
        k = idx + 1; % first arc of the i-th row
        l = idx + cols - 1; % last arc of the i-th row

        % update matrix B
        B(i:j - 1, k:l) = alpha * eye(cols - 1);
        B(i + 1:j, k:l) = B(i + 1:j, k:l) - alpha * eye(cols - 1);

        % vertical
        if r < rows - 1
            k = l + 1; % first arc of the i-th row
            l = l + cols; % last arc of the i-th row

            % update matrix B
            B(i:j, k:l) = alpha * eye(cols);
            B(i + cols:j + cols, k:l) = -alpha * eye(cols);
        end
    end

    %% Add 'cols' columns (ground) or 1 column (end-to-end)
    switch method
        case 0% ground
            B(:, end + cols) = 0;
            B(end - (cols - 1):end, end - (cols - 1):end) = -eye(cols);

        otherwise % end-to-end
            B(:, end + cols) = 0;
            B(end, end) = -1;
    end
