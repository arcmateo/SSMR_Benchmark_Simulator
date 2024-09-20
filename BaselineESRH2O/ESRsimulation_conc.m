% Ethanol Steam Reformer (with finite differences)
% Model in terms of Concentrations

clear; close all; clc;

Parameters % Load the parameters

% Load initial conditions
ss_filename = ['SS_files3\SS_u_1_np_',num2str(np),'.mat'];
load(ss_filename);
%x0_1c = Flows_to_Conc(x0_1);
% Transforms initial conditions from Flows to Concentrations

options = odeset('RelTol', 1e-4,'AbsTol', 1e-5,'MaxStep', 1,...
    'NonNegative', 1:8*np); % Options for the Solver

% Configure solver and launch simulation

t_interv = [0 1]; % [min] Simulation time

tic
[t1,x1] = ode15s(@(t,x)ESR_conc(t,x,u_ss), t_interv, x0_1, options);
toc

disp("First stage done")

figure(1)
plot(1:np, x1(end, 3*np+1:4*np), linewidth=2)
title('Hydrogen concentration profile')
ylabel('Concentration [mol/m3]')
xlabel('Position')
grid on

figure(2)
plot(t1, x1(:,4*np), linewidth=2)
title('Hydrogen concentration at the reactor outlet')
ylabel('Concentration [mol/m3]')
xlabel('Time [min]')
grid on
