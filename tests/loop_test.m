% loop_test.m
% Simulation useful to reproduce the time varying behaviour of the dynamic
% system running the mathematical model associated with said system used to
% explain why the lightnings always choose a specific single path.
%
% Author: David Palma, PhD, DE <david.palma@uniud.it>

clear
close 
clc

%% Add 'lib' and 'src' folders and subfolders to path
addpath(genpath('../lib'));
addpath(genpath('../src'));

%% Test
% setup the properties
r = 100; % rows
c = 30; % columns
T = 10; % total time
debug = 1; % debug mode
th = 1; % 0 = monomial, 1 = piecewise linear, 2= linear

sim_cnt = 0; % counter

% run the model

i = 0.8;  % This is already \delta*2 (use i = 0 for a homogeneous environment)

    tic
    mode = 0; % mode: ground (0), end-to-end (1)
    sim_cnt = sim_cnt + 1;

    fprintf('Test no. %d (mode: end-to-end)\n', sim_cnt);

    obj = dynamic_system_simulation(r, c, mode, i, T, debug, th);
    toc

