% Ethanol Steam Reformer (with finite differences)
% Model in terms of Concentrations
clear; close all; clc;
global F_H2

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

t = 1; % [min] Simulation overall time
t_s = 0.1; % [min] Sampling time
time = 0:t_s:t;
y_output = zeros(size(time));

tic
u_ss = [0.0000, 0.0088];
for k = 1:length(time)

   [t,x] = ode15s(@(t,x)SSMR_function(t,x,u_ss,p), [0 t_s], x0c, options);
  
    y_output(k) = F_H2;
   
   x0c = x(end,:);
end
toc

figure(1)
plot(time, y_output, linewidth=2)
title('Pure hydrogen molar flow')
ylabel('Molar flow rate [mol/min]')
xlabel('Time [min]')
grid on

x0c = x(end,:);
mode = 3;

filename = ['ICH2O\Mode',num2str(mode),'_np', num2str(np),'.mat'];
save(filename, 'x0c', 'u_ss', 'np');
disp(['Stored file: ', filename])

