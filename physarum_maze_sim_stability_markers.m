%% PHYSARUM CIRCUIT-NETWORK MAZE SIMULATION
%
% Implements the dynamical model of Physarum polycephalum's maze-solving
% network described by Eq. 7.
%
% The equation is integrated as a MASS-MATRIX ODE
%       A * dv/dt = -( B*M(B'*v) - dbar ),      A = B*C*B'  (constant)
% to boost MATLAB's stiff solver ode15s.
%
% The script:
%   1) builds an example maze on a regular 2D grid (4-neighbour lattice),
%   2) maps it into the graph / incidence matrix B 
%   3) integrates the equation forward in time,
%   4) saves a snapshot of the potential field v(t), AND of the
%      deviation-from-equilibrium variable x(t) = v(t) - v_bar.
%

clear; close all; clc;
rng(2);  % seed

%% ------------------------- USER PARAMETERS -------------------------

% --- grid / maze size 
% A proper labyrinth (see generate_labyrinth_maze below) is carved on a
% grid of "cells" separated by single-cell-wide walls.
% Non-odd values are rounded to the nearest valid size automatically.

Nrows = 41;                 % number of grid rows
Ncols = 27;                 % number of grid columns

% --- maze "braiding"
% braidProb randomly knocks down a small fraction of the remaining dead-
% end walls to add a few loops/shortcuts; 0 = perfect maze (no loops).
braidProb = 0.05;

% --- circuit / dynamics parameters 
Ck_value   = 1;             % capacitance per link, C_k = 1 default
alpha      = 1e-5;          % lower slope of the characteristic function M_k, eq.(4)
beta       = 800;           % upper (steep) slope of M_k -> close to ideal threshold
VT_mean    = 0.5;           % mean activation threshold V_Tk
heteroVT   = true;          % true -> heterogeneous thresholds (scenario (c)/(d))
delta_VT   = 0.5;           % half-width of the uniform distribution for V_Tk
d_amplitude = 1;            % magnitude of the input flow at the maze entry

% --- simulation horizon and snapshot times 
T_final = 400;                                  % total simulated time
snapshot_times = [0, 1, 5, 20, 60, 150, 300, T_final];  % <-- user editable (times to take snapshots)

% --- output 
outDir = fullfile(pwd,'physarum_snapshots_het');
if ~exist(outDir,'dir'); mkdir(outDir); end

%% ------------------------- BUILD THE MAZE -------------------------
% wallMask(r,c) = true  -> cell is an obstacle (no node there)
%              = false -> free cell (becomes a graph node)
[wallMask, Nrows, Ncols] = generate_labyrinth_maze(Nrows, Ncols, braidProb);

% Node map: nodeID(r,c) = index of the node at free cell (r,c), 0 if wall
nodeID = zeros(Nrows,Ncols);
freeIdx = find(~wallMask);
nodeID(freeIdx) = 1:numel(freeIdx);
n = numel(freeIdx);   % number of nodes

% Ensure the maze graph (internal 4-neighbour connectivity) is connected;
% if not, keep only the largest connected component as the usable maze.
[wallMask, nodeID, n] = enforce_connected(wallMask, nodeID);

%% ------------------------- ENTRY / EXIT NODES -------------------------
% Entry: a single free cell on the top row (Physarum's starting point,
% "node 1" of the paper). We search outward from the top-centre until a
% free cell is found.
entry_r = 1;
entry_c = pick_free_in_row(wallMask, entry_r, ceil(Ncols/2));
entry_node = nodeID(entry_r, entry_c);

% Exit / food source: spread over the whole bottom row.
bottom_r = Nrows;
exit_cols  = find(nodeID(bottom_r,:) > 0);
exit_nodes = nodeID(bottom_r, exit_cols)';

%% ------------------------- BUILD INCIDENCE MATRIX B -------------------------
% Internal links: 4-neighbour connections between free cells .
% External links: one per entry (into entry_node) and one per exit node
% (out of that node).

linkList = zeros(0,2);   % [startNode, endNode], internal links only (B_hk=1,B_rk=-1)
for r = 1:Nrows
    for c = 1:Ncols
        if wallMask(r,c); continue; end
        h = nodeID(r,c);
        % right neighbour
        if c < Ncols && ~wallMask(r,c+1)
            linkList(end+1,:) = [h, nodeID(r,c+1)]; %#ok<AGROW>
        end
        % down neighbour
        if r < Nrows && ~wallMask(r+1,c)
            linkList(end+1,:) = [h, nodeID(r+1,c)]; %#ok<AGROW>
        end
    end
end

m_internal = size(linkList,1);
m_entry    = 1;                     % single input link into the maze
m_exit     = numel(exit_nodes);     % one output link per bottom node
m = m_internal + m_entry + m_exit;

rows_ = zeros(2*m_internal + m_entry + m_exit, 1);
cols_ = zeros(size(rows_));
vals_ = zeros(size(rows_));
p = 0;

% internal links: B_hk = +1 (start), B_rk = -1 (end)
for k = 1:m_internal
    h = linkList(k,1); r_ = linkList(k,2);
    p = p+1; rows_(p)=h;  cols_(p)=k; vals_(p)=  1;
    p = p+1; rows_(p)=r_; cols_(p)=k; vals_(p)= -1;
end

% entry link: comes from the external environment into entry_node
col_entry = m_internal + 1;
p = p+1; rows_(p)=entry_node; cols_(p)=col_entry; vals_(p) = -1;

% exit links: leave the maze (bottom nodes) toward the food source
for j = 1:m_exit
    k = m_internal + m_entry + j;
    p = p+1; rows_(p)=exit_nodes(j); cols_(p)=k; vals_(p)=1;
end

rows_ = rows_(1:p); cols_ = cols_(1:p); vals_ = vals_(1:p);
B = sparse(rows_, cols_, vals_, n, m);

%% ------------------------- OTHER INGREDIENTS
% CAPACITANCE MATRIX C 
C = Ck_value * speye(m);

% THRESHOLDS V_Tk (possibly heterogeneous)
if heteroVT
    VT = VT_mean - delta_VT + 2*delta_VT*rand(m,1);   % U[VT_mean-delta, VT_mean+delta]
    VT = max(VT, 1e-3);                                % keep strictly positive
else
    VT = VT_mean * ones(m,1);
end

% INPUT FLOW VECTOR dbar 
dbar = zeros(n,1);
dbar(entry_node) = d_amplitude;

%% ------------------------- MASS MATRIX & ODE SETUP -------------------------
A = B * C * B.';           % constant, sparse, symmetric positive definite

odefun = @(t, v) -( B * charFunc_M(B.' * v, VT, alpha, beta) - dbar );

opts = odeset('Mass', A, 'MassSingular', 'no', 'MStateDependence', 'none', ...
              'RelTol', 1e-6, 'AbsTol', 1e-8, ...
              'Jacobian', @(t,v) odeJacobian(v, B, VT, alpha, beta));

v0 = zeros(n,1);

fprintf('Simulating Physarum network: %d nodes, %d links...\n', n, m);
tic;
[tOut, vOut] = ode15s(odefun, snapshot_times, v0, opts);
fprintf('Done in %.2f s.\n', toc);

%% ------------------------- STEADY STATE v_bar & DEVIATION x(t) -------------

fprintf('Continuing integration to estimate the steady state v_bar...\n');
convOpts = odeset(opts, 'Events', @(t,v) steadyStateEvent(t, v, odefun));
vLast = vOut(end,:).';
tExtra = max(10*T_final, 10);
[~, ~, ~, vEnd, ~] = ode15s(odefun, [tOut(end), tOut(end)+tExtra], vLast, convOpts);
if ~isempty(vEnd)
    vbar = vEnd(end,:).';
else
    % did not fully converge within tExtra: fall back to the last state
    [~, vTail] = ode15s(odefun, [tOut(end), tOut(end)+tExtra], vLast, opts);
    vbar = vTail(end,:).';
end

xOut = vOut - vbar.';   % deviation variable, one row per snapshot time

%% ------------------------- SAVE SNAPSHOTS -------------------------

nSnap = numel(tOut);
snapshots(nSnap) = struct('t', [], 'grid_v', [], 'grid_x', []);

for s = 1:nSnap
    v_s = vOut(s,:).';
    x_s = xOut(s,:).';

    grid_v = nan(Nrows, Ncols);
    grid_x = nan(Nrows, Ncols);
    grid_v(nodeID > 0) = v_s(nodeID(nodeID > 0));
    grid_x(nodeID > 0) = x_s(nodeID(nodeID > 0));

    snapshots(s).t      = tOut(s);
    snapshots(s).grid_v = grid_v;
    snapshots(s).grid_x = grid_x;

    % --- v(t) figure ---
    fh = figure('Visible','off','Color','w');
    imagesc(grid_v, 'AlphaData', ~isnan(grid_v));
    set(gca,'Color',[0.15 0.15 0.15]);
    axis image off; colormap(parula); colorbar;
    hold on;
    plot(entry_c, entry_r, 'p', 'MarkerSize', 16, 'MarkerFaceColor', [0 1 0.3], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);
    plot(exit_cols, bottom_r*ones(size(exit_cols)), 's', 'MarkerSize', 8, ...
        'MarkerFaceColor', [1 0.15 0.15], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
    legend({'entry','exit (food source)'}, 'Location','southoutside', ...
        'Orientation','horizontal', 'TextColor','k', 'Color','w');
    hold off;
    title(sprintf('Potential v(t),  t = %.3g', tOut(s)));
    fname_v = fullfile(outDir, sprintf('snapshot_v_t%03d_%.3f.png', s, tOut(s)));
    exportgraphics(fh, fname_v, 'Resolution', 150);
    close(fh);

    % --- x(t) = v(t) - v_bar figure (deviation from equilibrium) ---
    fh = figure('Visible','off','Color','w');
    imagesc(grid_x, 'AlphaData', ~isnan(grid_x));
    set(gca,'Color',[0.15 0.15 0.15]);
    axis image off; colormap(redblue_cmap()); colorbar;
    clim_ = max(abs(x_s)); if clim_==0; clim_=1; end
    clim([-clim_, clim_]);
    hold on;
    plot(entry_c, entry_r, 'p', 'MarkerSize', 16, 'MarkerFaceColor', [0 1 0.3], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);
    plot(exit_cols, bottom_r*ones(size(exit_cols)), 's', 'MarkerSize', 8, ...
        'MarkerFaceColor', [1 0.15 0.15], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
    legend({'entry','exit (food source)'}, 'Location','southoutside', ...
        'Orientation','horizontal', 'TextColor','k', 'Color','w');
    hold off;
    title(sprintf('Deviation x(t) = v(t)-v_{bar},  t = %.3g', tOut(s)));
    fname_x = fullfile(outDir, sprintf('snapshot_x_t%03d_%.3f.png', s, tOut(s)));
    exportgraphics(fh, fname_x, 'Resolution', 150);
    close(fh);
end

save(fullfile(outDir,'physarum_simulation.mat'), ...
     'snapshots','tOut','vOut','xOut','vbar','B','C','VT','alpha','beta', ...
     'nodeID','wallMask','entry_node','entry_r','entry_c','exit_nodes','exit_cols','bottom_r');

fprintf('Saved %d snapshots (v & x, .png) and simulation data (.mat) to:\n  %s\n', nSnap, outDir);
fprintf('||x(t_end)|| = %.3g  (should be small: convergence to v_bar)\n', norm(xOut(end,:)));

%% ========================================================================
%                           LOCAL FUNCTIONS
%% ========================================================================

function [wallMask, Nrows, Ncols] = generate_labyrinth_maze(Nrows, Ncols, braidProb)
% Builds a genuine labyrinth: single-cell-wide winding corridors and
% dead ends, generated with the classic recursive-backtracker ("depth-
% first search") perfect-maze algorithm, optionally "braided" with a
% few extra passages to add short loops.

    % --- snap to valid odd dimensions -----------------------------------
    % cellRows cells need 2*cellRows-1 fine-grid rows (cells at rows
    % 1,3,5,...,2*cellRows-1, with a wall/passage row between each pair).
    cellRows = max(2, round((Nrows+1)/2));
    cellCols = max(2, round((Ncols+1)/2));
    Nrows = 2*cellRows - 1;
    Ncols = 2*cellCols - 1;

    wallMask = true(Nrows, Ncols);          % start fully walled
    visited  = false(cellRows, cellCols);

    % carve out every cell centre
    for i = 1:cellRows
        for j = 1:cellCols
            wallMask(2*i-1, 2*j-1) = false;
        end
    end

    % --- recursive backtracker (iterative, explicit stack) --------------
    start = [1, ceil(cellCols/2)];
    stack = start;
    visited(start(1), start(2)) = true;

    while ~isempty(stack)
        cur = stack(end,:);
        i = cur(1); j = cur(2);

        % unvisited 4-neighbours in CELL space
        nbrs = [i-1,j; i+1,j; i,j-1; i,j+1];
        valid = nbrs(:,1)>=1 & nbrs(:,1)<=cellRows & ...
                nbrs(:,2)>=1 & nbrs(:,2)<=cellCols;
        nbrs = nbrs(valid,:);
        unvis = false(size(nbrs,1),1);
        for kk = 1:size(nbrs,1)
            unvis(kk) = ~visited(nbrs(kk,1), nbrs(kk,2));
        end
        nbrs = nbrs(unvis,:);

        if isempty(nbrs)
            stack(end,:) = [];      % dead end -> backtrack
        else
            pick = nbrs(randi(size(nbrs,1)),:);
            % knock down the wall on the fine grid:
            wr = i + pick(1) - 1; wc = j + pick(2) - 1;
            wallMask(wr, wc) = false;
            visited(pick(1), pick(2)) = true;
            stack(end+1,:) = pick; %#ok<AGROW>
        end
    end

    % --- optional braiding: remove a few extra walls at dead end to create short loops
    if braidProb > 0
        for i = 1:cellRows
            for j = 1:cellCols
                gr = 2*i-1; gc = 2*j-1;
                openCount = 0;
                candidates = zeros(0,2);
                nbrs = [i-1,j; i+1,j; i,j-1; i,j+1];
                for kk = 1:size(nbrs,1)
                    ni = nbrs(kk,1); nj = nbrs(kk,2);
                    if ni<1 || ni>cellRows || nj<1 || nj>cellCols; continue; end
                    wr = i+ni-1; wc = j+nj-1;
                    if ~wallMask(wr,wc)
                        openCount = openCount + 1;
                    else
                        candidates(end+1,:) = [wr, wc]; %#ok<AGROW>
                    end
                end
                if openCount == 1 && ~isempty(candidates) && rand < braidProb
                    pick = candidates(randi(size(candidates,1)),:);
                    wallMask(pick(1), pick(2)) = false;
                end
            end
        end
    end

    % Note: row 1 and row Nrows each contain `cellCols` free cells to serve 
    % as entry/exit points; forcing the whole row open would break the 
    % labyrinth walls right at the boundary.
end

function [wallMask, nodeID, n] = enforce_connected(wallMask, nodeID)
% Keeps only the largest 4-connected component of free cells (B full rank)
    [Nrows, Ncols] = size(wallMask);
    freeMask = ~wallMask;
    CC = bwconncomp_manual(freeMask);
    [~, biggest] = max(cellfun(@numel, CC));
    keepMask = false(Nrows,Ncols);
    keepMask(CC{biggest}) = true;

    wallMask = ~keepMask;
    nodeID = zeros(Nrows,Ncols);
    freeIdx = find(keepMask);
    nodeID(freeIdx) = 1:numel(freeIdx);
    n = numel(freeIdx);
end

function CC = bwconncomp_manual(freeMask)
% Minimal 4-connectivity connected-components labelling (BFS)
    [Nrows, Ncols] = size(freeMask);
    visited = false(Nrows,Ncols);
    CC = {};
    for r = 1:Nrows
        for c = 1:Ncols
            if freeMask(r,c) && ~visited(r,c)
                stack = [r,c];
                visited(r,c) = true;
                comp = [];
                while ~isempty(stack)
                    cur = stack(end,:); stack(end,:) = [];
                    comp(end+1) = sub2ind([Nrows,Ncols], cur(1), cur(2)); 
                    nbrs = [cur(1)-1,cur(2); cur(1)+1,cur(2); cur(1),cur(2)-1; cur(1),cur(2)+1];
                    for kk = 1:4
                        rr = nbrs(kk,1); cc = nbrs(kk,2);
                        if rr>=1 && rr<=Nrows && cc>=1 && cc<=Ncols && freeMask(rr,cc) && ~visited(rr,cc)
                            visited(rr,cc) = true;
                            stack(end+1,:) = [rr,cc]; 
                        end
                    end
                end
                CC{end+1} = comp;
            end
        end
    end
end

function c = pick_free_in_row(wallMask, r, cPreferred)
% Returns a free column in row r, starting from cPreferred and searching
% outward if that cell happens to be a wall.
    Ncols = size(wallMask,2);
    if ~wallMask(r, cPreferred)
        c = cPreferred; return;
    end
    for offset = 1:Ncols
        for cand = [cPreferred-offset, cPreferred+offset]
            if cand >= 1 && cand <= Ncols && ~wallMask(r,cand)
                c = cand; return;
            end
        end
    end
    error('No free cell found on row %d.', r);
end

function [value, isterminal, direction] = steadyStateEvent(t, v, odefun)
% Stops the integration once the RHS has become negligible, i.e. v(t) has
% practically reached the steady state v_bar.
    f = odefun(t, v);
    value = norm(f) - 1e-6*max(1,norm(v));
    isterminal = 1;
    direction  = -1;
end

function cmap = redblue_cmap()
% Simple diverging blue-white-red colormap .
    n2 = 128;
    top = [linspace(1,0.85,n2)', linspace(0,0.1,n2)', linspace(0,0.1,n2)'];
    bot = [linspace(0.1,1,n2)', linspace(0.1,0,n2)', linspace(0.85,1,n2)'];
    cmap = [flipud(bot); top];
end

function y = charFunc_M(u, VT, alpha, beta)
% Piecewise-linear approximation of the ideal threshold characteristic
% function.
    y = beta.*u - 0.5*(beta-alpha).*( abs(u+VT) - abs(u-VT) );
end

function J = odeJacobian(v, B, VT, alpha, beta)
% Analytical Jacobian of the RHS f(v) = -(B*M(B'*v) - dbar) w.r.t. v,
% supplied to ode15s (with the mass matrix A).
    u = B.' * v;
    dM = alpha + (beta-alpha) * (abs(u) > VT);   % per-link derivative
    J = -B * spdiags(dM, 0, numel(dM), numel(dM)) * B.';
end
