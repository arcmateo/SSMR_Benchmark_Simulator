% Ethanol Steam Reformer (with finite differences)
% Model in terms of Concentrations
clear; close all; clc;
global F_H2

p = Parameters(); % Load the parameters 
np = p.np; % Number of points (spatial discretization)

addpath('ICFull');

%% Simulation setup

ss_filename = ['Mode1','_np',num2str(np),'.mat'];
load(ss_filename);

%%
options = odeset('RelTol', 1e-4,'AbsTol', 1e-5,'MaxStep', 0.1,...
    'NonNegative', 1:8*2*np); % Options for the Solver

% Configure solver and launch simulation

t = 30; % [min] Simulation overall time
t_s = 0.1; % [min] Sampling time
time = 0:t_s:t;
u_before = u_ss(1);
sp = 9.75437e-4; % [mol/min]
y_sp = ones(size(time))*sp; % [mol/min]

y_output = zeros(size(time));
u_output = zeros(size(time));

ku = 0.376;
tao = 0.25;

%Ziegler-Nichols method
kp = 2.1*(0.6*ku);
ki = 0.01*((1.2*ku)/tao);
kd = 0.001*(3*ku*tao)/40;

integral_error = 0;
prev_error = 0;

tic
for k = 1:length(time)

   [t,x] = ode15s(@(t,x)SSMR_function(t,x,u_ss,p), [0 t_s], x0c, options);
  
    y_output(k) = F_H2;

   if k >= 1 && k <= 100
      y_sp(k) = sp;
   elseif k > 100 && k <= 150
      y_sp(k) = sp*1.05;
   elseif k > 150 && k <= 200
      y_sp(k) = sp;
   elseif k > 200 && k <= 250
      y_sp(k) = sp*0.9;
   elseif k > 250 && k <= 301
      y_sp(k) = sp;
   end

   error = y_sp(k) - y_output(k);

   integral_error = integral_error + error*t_s;
   derivative_error = (error - prev_error)/t_s;
   deltau = kp*error + ki*integral_error + kd*derivative_error;

   u = u_before + deltau;

   if u <= 0.0018
    u = 0.0018;
   end

   if u >= 0.0024
      u = 0.0024;
   end

   u_ss(1) = u;
   u_before = u;
   u_output(k) = u_before;

   x0c = x(end,:);
   prev_error = error;

end
toc

figure;
subplot(2,1,1);
plot(time, y_output, 'b', 'LineWidth', 1.5);
hold on;
stairs(time, y_sp, 'r--', 'LineWidth', 1.0);
xlabel('Time (s)'); ylabel('Output y');
title('System Response');
legend('y(t)', 'Setpoint');
grid on

subplot(2,1,2);
stairs(time, u_output, 'k', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Control Input u(t)');
title('Control Signal');
grid on

