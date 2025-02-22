% Ethanol Steam Reformer (with finite differences)
% Model in terms of Concentrations
clear; close all; clc;
global F_H2

addpath('ICFull','ICH2O');

%% Simulation setup

% Select the number of points (spatial discretization): 50 or 200
np = 50;

% Select the normal operating conditions: Mode 1, 2 or 3
Mode = 1;

switch Mode
   case 1
   P_in = 4.0; % [bar]
   T_in = 773.15; % [K]
   case 2
   P_in = 6.0; % [bar]
   T_in = 823.15; % [K]
   case 3
   P_in = 8.0; % [bar]
   T_in = 873.15; % [K]
end

p = Parameters(P_in, T_in, np); % Load the parameters 

% Select the initial conditions:
% 0 = steady state with all compunds, 1 =  steady state full of steam
initial_conditions = 0; 

switch initial_conditions
    case 0 
        ss_filename = ['Mode',num2str(Mode),'_np',num2str(np),'.mat'];
    case 1 
        ss_filename = ['Mode',num2str(Mode),'_np',num2str(np),'_H2O','.mat'];
end
load(ss_filename); 


%%
options = odeset('RelTol', 1e-4,'AbsTol', 1e-5,'MaxStep', 0.1,...
    'NonNegative', 1:8*2*np); % Options for the Solver

% Configure solver and launch simulation

% Select the overall simulation time
t = 30; % [min] between 10 and 30 min are recommended

% Select the sampling time
t_s = 0.1; % [min] 

% Select the control law, open loop = 0, PID = 1
control_law = 1;

time = 0:t_s:t;
y_output = zeros(size(time));
u_output = zeros(size(time));
sp = 9.75e-4; % [mol/min]
y_sp = Profile(sp, time, t_s);

switch control_law
   case 0
      tic
      for k = 1:length(time)
      [t,x] = ode15s(@(t,x)SSMR_function(t,x,u_ss,p), [0 t_s], x0c, options);
      y_output(k) = F_H2;
      u_output(k) = u_ss(1);
      x0c = x(end,:);
      end
      toc
      figure(1)
      plot(time, y_output, 'b', 'LineWidth', 1.5);
      xlabel('Time (min)'); 
      ylabel('Pure hydrogen molar flow [mol/min]');
      grid on
   case 1
      ku = 0.376;
      tao = 0.25;
      %Ziegler-Nichols method
      kp = 2.1*(0.6*ku);
      ki = 0.006*((1.2*ku)/tao);
      kd = 0.001*(3*ku*tao)/40;
      integral_error = 0;
      prev_error = 0;
      tic
      for k = 1:length(time)
      u_output(k) = u_ss(1);
      [t,x] = ode15s(@(t,x)SSMR_function(t,x,u_ss,p), [0 t_s], x0c, options);
      y_output(k) = F_H2;
      error = y_sp(k) - y_output(k);
      integral_error = integral_error + error*t_s;
      derivative_error = (error - prev_error)/t_s;
      deltau = kp*error + ki*integral_error + kd*derivative_error;
      u_ss(1) = u_ss(1) + deltau;
      if u_ss(1) <= 0.0018
         u_ss(1) = 0.0018;
      end
      if u_ss(1) >= 0.0024
         u_ss(1) = 0.0024;
      end
      prev_error = error;
      x0c = x(end,:);
      end
      toc
      figure;
      subplot(2,1,1);
      plot(time, y_output, 'b', 'LineWidth', 1.5);
      hold on;
      stairs(time, y_sp, 'r--', 'LineWidth', 1.0);
      xlabel('Time (min)'); 
      ylabel('Pure H2 flow [mol/min]');
      title('System response');
      legend('Pure H2 flow', 'Setpoint');
      grid on  
      subplot(2,1,2);
      stairs(time, u_output, 'k', 'LineWidth', 1.5);
      xlabel('Time (min)'); 
      ylabel('Inlet ethanol flow [mol/min]');
      title('Control Input');
      grid on
end

