% Ethanol Steam Reformer (with finite differences)
% Model in terms of Concentrations
clear; close all; clc;
global F_out

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

%%

options = odeset('RelTol', 1e-4,'AbsTol', 1e-5,'MaxStep', 1,...
    'NonNegative', 1:8*np); % Options for the Solver

% Configure solver and launch simulation

t = 1; % [min] Simulation overall time

t_interv = [0 t]; % [min] Simulation time interval

tic

[t1,x1] = ode15s(@(t,x)ESR_conc1stg(t,x,u_ss,p), t_interv, x0_1c, options);

toc

disp("First stage done")

% Second stage

x_boundary = x1(:,np:np:end);

tic
[t2,x2] = ode15s(@(t,x)ESR_conc2stg(t,x,x_boundary_interp(x_boundary,t1,t)',p,F_out),...
    t_interv, x0_2c, options);
toc
disp("Second stage done")

figure(1)
plot(1:2*np, [x1(end,3*np+1:4*np), x2(end,3*np+1:4*np)], linewidth=2)
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

figure(3)
plot(t2, x2(:,4*np), linewidth=2)
title('Hydrogen concentration at the retentate outlet')
ylabel('Concentration [mol/m3]')
xlabel('Time [min]')
grid on

