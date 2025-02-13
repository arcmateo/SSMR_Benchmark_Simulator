% Ethanol Steam Reformer (with finite differences)
% Model in terms of Concentrations
clear; close all; clc;

p = Parameters(); % Load the parameters 
np = p.np; % Number of points (spatial discretization)
addpath('SS_files3','SS_filesH2O');

%% Simulation setup

initial_conditions = 0;
switch initial_conditions
    case 0 % ESR in steady-state with all compunds
        ss_filename = ['SS_u_1_np_',num2str(np),'.mat'];
        load(ss_filename);
    case 1 % ESR in steady-state full of steam
        ss_filename = ['SS_u_1_np_H2O_',num2str(np),'.mat'];
        load(ss_filename);
        x0_1c = Flows_to_Conc(x0_1,p);
        x0_2c = x0_1c;     
end
x0c = [x0_1c, x0_2c];


%%
options = odeset('RelTol', 1e-4,'AbsTol', 1e-5,'MaxStep', 0.1,...
    'NonNegative', 1:8*2*np); % Options for the Solver

% Configure solver and launch simulation

t = 0.4; % [min] Simulation overall time

t_interv = [0 t]; % [min] Simulation time interval

tic
[t,x] = ode15s(@(t,x)SSMR_function(t,x,u_ss,p), t_interv, x0c, options);
toc

figure(1)
plot(t, x(:,4*np), linewidth=2)
title('Hydrogen concentration at the reactor outlet')
ylabel('Concentration [mol/m3]')
xlabel('Time [min]')
grid on



