% Ethanol Steam Reformer (with finite differences)
% Model in terms of Molar flows
clear; close all; clc;

p = Parameters(); % Load the parameters 
np = p.np; % Number of points (spatial discretization)
addpath('SS_files2','SS_filesH2O');

%% Simulation setup

initial_conditions = 0;
switch initial_conditions
    case 0 % ESR in steady-state with all compunds
        ss_filename = ['SS_u_1_np_',num2str(np),'.mat'];
        load(ss_filename);
    case 1 % ESR in steady-state full of steam
        ss_filename = ['SS_u_1_np_H2O_',num2str(np),'.mat'];
        load(ss_filename);
        x0_2 = x0_1;
end

%%

options = odeset('RelTol', 1e-4,'AbsTol', 1e-5, ...
    'MaxStep', 0.1, 'NonNegative', 1:8*np); % Options for the Solver

% Configure solver and launch simulation

t = 1; % [min] Simulation overall time

t_interv = [0 t]; % [min] Simulation time interval

tic
[t1,x1] = ode15s(@(t,x)ESR_flows1stg(t,x,u_ss,p), t_interv, x0_1, options);
toc
disp("First stage done")

% Second stage
x_boundary = x1(:,np:np:end);

tic
[t2,x2] = ode15s(@(t,x)ESR_flows2stg(t,x,x_boundary_interp(x_boundary,t1,t)',p),...
    t_interv, x0_2, options);
toc
disp("Second stage done")


figure(1)
plot(1:2*np, [x1(end,3*np+1:4*np), x2(end,3*np+1:4*np)], linewidth=2)
title('Hydrogen molar flow profile')
ylabel('Molar flow [mol/min]')
xlabel('Position')
grid on

figure(2)
plot(t1, x1(:,4*np), linewidth=2)
title('Hydrogen molar flow at the reactor outlet')
ylabel('Molar flow [mol/min]')
xlabel('Time [min]')
grid on

figure(3)
plot(t2, x2(:,4*np), linewidth=2)
title('Hydrogen molar flow at the retentate outlet')
ylabel('Molar flow [mol/min]')
xlabel('Time [min]')
grid on

