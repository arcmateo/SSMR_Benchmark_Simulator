% Staged-separation membrane reactor (with finite differences)
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
   ss = 2.27354e-4; % [mol/min]
   case 2
   P_in = 6.0; % [bar]
   T_in = 823.15; % [K]
   ss = 2.76379e-4; % [mol/min]
   case 3
   P_in = 8.0; % [bar]
   T_in = 873.15; % [K]
   ss = 5.81357e-4; % [mol/min]
end

p = Parameters(P_in, T_in, np); % Load the parameters 

% Select the initial conditions:
%   0 = steady state, reactor contains all compounds
%   1 = steady state, reactor contains only steam
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


% Select the overall simulation time
t = 30; % [min] - recommended: between 10 and 30 min

% Select the sampling time
t_s = 0.1; % [min] 

% Select the control law:
% 0 = open loop 
% 1 = PID 
control_law = 1;

% Select the set-point profile type
% 1 = set-point profile 1
% 2 = set-point profile 2
type = 1;  

time = 0:t_s:t;
y_output = zeros(size(time));
u_output = zeros(size(time));
y_sp = Profile(ss, time, t_s, type);

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
      figure;
      subplot(2,1,1);
      plot(time, y_output, 'b', 'LineWidth', 2.0);
      set(gca,'FontSize', 16)
      hold on;
      xlabel('Time (min)', fontsize = 16); 
      ylabel('Pure H_{2} flow (mol/min)', fontsize = 16);
      title('System response', fontsize = 16);
      grid on  
      subplot(2,1,2);
      stairs(time, u_output, 'k', 'LineWidth', 2.0);
      set(gca,'FontSize', 16)
      xlabel('Time (min)', fontsize = 16); 
      ylabel('Inlet ethanol flow (mol/min)', fontsize = 16);
      title('Control input', fontsize = 16);
      grid on

   case 1
      ku = 0.1265;
      tao = 0.25;
      %Ziegler-Nichols method
      kp = 18*(0.6*ku);
      ki = 0.01*((1.2*ku)/tao);
      kd = (3*ku*tao)/40;
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
      plot(time, y_output, 'b', 'LineWidth', 2.0);
      set(gca,'FontSize', 16)
      hold on;
      stairs(time, y_sp, 'r--', 'LineWidth', 2.0);
      xlabel('Time (min)',fontsize = 16); 
      ylabel('Pure H_{2} flow (mol/min)',fontsize = 16); 
      title('System response',fontsize = 16);
      legend('Pure H_{2} molar flow', 'Set-point',fontsize = 16);
      grid on  
      subplot(2,1,2);
      stairs(time, u_output, 'k', 'LineWidth', 2.0);
      set(gca,'FontSize', 16)
      xlabel('Time (min)',fontsize = 16);
      ylabel('Inlet ethanol flow (mol/min)',fontsize = 16);
      title('Control input',fontsize = 16);
      grid on
end

