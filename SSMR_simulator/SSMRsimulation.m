% Ethanol Steam Reformer (with finite differences)
% Model in terms of Concentrations
clear; close all; clc;
global vel1 vel2

p = Parameters(); % Load the parameters 
np = p.np; % Number of points (spatial discretization)
A = p.A; % [m2] Reactor cross-sectional area 
Am = p.Am; % [m2] Membrane cross-sectional area

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

t = 0.8; % [min] Simulation overall time

t_interv = [0 t]; % [min] Simulation time interval

tic
[t,x] = ode15s(@(t,x)SSMR_function(t,x,u_ss,p), t_interv, x0c, options);
toc

Q_out1 = A*vel1; % [m3/min]
Q_out2 = (A-Am)*vel2; % [m3/min]
F_H2_out1 = Q_out1*x(:,4*np); % [mol/min]
F_H2_out2 = Q_out2*x(:,12*np); % [mol/min]
F_H2_pure = F_H2_out1 - F_H2_out2; % [mol/min]

figure(1)
plot(t, F_H2_pure, linewidth=2)
title('Hydrogen flow')
ylabel('Molar flow rate [mol/min]')
xlabel('Time [min]')
grid on



