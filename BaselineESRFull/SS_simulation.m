% Ethanol Steam Reformer 
% Steady-state simulation for initial conditions
clear; close all; clc;

Parameters % Load the parameters

factor = 1;

u_ss = [0.0021; 0.0099]; % [mol/min] steady-state inputs
% order: [C2H5OH, H2O]

% SIMULATION of SS for the 1st stage 
x_ss_1 = [u_ss; zeros(ns-2,1); T0]; % Set boundary conditions (k = 0)

options = odeset('MaxStep', 1, 'NonNegative', 1:8); % Options for the Solver

tic
[z1,x0_1] = ode15s(@(z,x)SS_function1stg(z,x), (0:deltaz1:L1), x_ss_1, options);
toc

x0_1 = x0_1(2:end,:); % (without z = 0)
x0_1 = x0_1(:); % Only 1 column vector for initial conditions!


% SIMULATION of SS for the 1st stage 
x_ss_2 = x0_1(np:np:end); % Set boundary conditions (k = np)
% End of first stage

tic
[z2,x0_2] = ode15s(@(z,x)SS_function2stg(z,x), (0:deltaz2:L2), x_ss_2, options);
toc

x0_2 = x0_2(2:end,:); % (without z = 0 of second stage))
x0_2 = x0_2(:); % only 1 column vector for initial conditions!


% Save solution into a file to load later on
filename = ['SS_files\SS_u_',num2str(factor),'_np_', num2str(np),'.mat'];
save(filename, 'x0_1', 'x0_2', 'u_ss', 'np');
disp(['Stored file: ', filename])

figure(1)
title('Molar flow rates across the reactor')
hold on
plot(1:2*np, [x0_1(0*np+1:1*np); x0_2(0*np+1:1*np)], linewidth=2)
plot(1:2*np, [x0_1(1*np+1:2*np); x0_2(1*np+1:2*np)], linewidth=2)
plot(1:2*np, [x0_1(2*np+1:3*np); x0_2(2*np+1:3*np)], linewidth=2)
plot(1:2*np, [x0_1(3*np+1:4*np); x0_2(3*np+1:4*np)], linewidth=2)
plot(1:2*np, [x0_1(4*np+1:5*np); x0_2(4*np+1:5*np)], linewidth=2)
plot(1:2*np, [x0_1(5*np+1:6*np); x0_2(5*np+1:6*np)], linewidth=2)
plot(1:2*np, [x0_1(6*np+1:7*np); x0_2(6*np+1:7*np)], linewidth=2)
legend('C2H5OH','H2O','CH4','H2','CO','CO2','CH3CHO');
ylabel('Molar flow rate [mol/min]')
xlabel('Position')
grid on

