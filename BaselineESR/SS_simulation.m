% Ethanol Steam Reformer 
% Steady-state simulation for initial conditions
clear; close all; clc;

% PARAMETERS ---------- 
% Same as in SS_function.m
factor = 1;

np = 600; % Number of points (spatial discretization)

T0 = 873.15; % [K] Inlet temperature

ns = 7; % Number of species

nr = 4; % Number of reactions

L = 0.23; % [m] Reactor length

L2 = 0.076; % [m] 2nd stage length

L1 = L - L2; % [m] 1st stage length

deltaz1 = L1 / np; % delta_z 1st stage
deltaz2 = L2 / np; % delta_z 2nd stage

u_ss = [0.0021; 0.0099]; % [mol/min] steady-state inputs
% order: [C2H5OH, H2O]

% SIMULATION of SS for the 1st stage 
x_ss_1 = [u_ss; zeros(ns-2,1); T0]; % Set boundary conditions (k = 0)

options = odeset('MaxStep', 1, 'NonNegative', 1:8); % Options for the Solver

tic
[z1,x0_1] = ode15s(@(z,x)SS_function(z,x), (0:deltaz1:L1), x_ss_1, options);
toc

x0_1 = x0_1(2:end,:); % (without z = 0)
x0_1 = x0_1(:); % Only 1 column vector for initial conditions!

% Save solution into a file to load later on
filename = ['SS_files\SS_u_', num2str(factor), '_np_', num2str(np), '.mat'];
save(filename, 'x0_1', 'u_ss', 'np');
disp(['Stored file: ', filename])

figure(1)
plot(1:np, x0_1(3*np+1:4*np), linewidth=2)
title('Hydrogen molar flow along the reactor')
ylabel('Molar flow [mol/min]')
xlabel('Position')
grid on