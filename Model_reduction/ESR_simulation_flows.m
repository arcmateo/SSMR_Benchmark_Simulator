% Ethanol Steam Reformer (with finite differences)
% Model in terms of Concentrations

clear; close all; clc;

% PARAMETERS ---------- 
% Same as in SS_function.m

np = 200; % Number of points (spatial discretization)

% Load initial conditions
ss_filename = 'SS_files\SS_u_1_np_200.mat';
load(ss_filename);

options = odeset('MaxStep', 1, 'NonNegative', 1:1600); % Options for the Solver

% Configure solver and launch simulation

t_interv = [0 2]; % [min] Simulation time

tic
[t1,x1] = ode15s(@(t,x)ESR_flows(t,x,u_ss), t_interv, x0_1, options);
toc

disp("First stage done")

figure(1)
plot(1:np, x1(end, 601:800), linewidth=2)
ylabel('H2 molar flow [mol/min]')
xlabel('Position')
grid on

figure(2)
plot(t1, x1(:,800), linewidth=2)
ylabel('H2 molar flow [mol/min]')
xlabel('Time [min]')
grid on


