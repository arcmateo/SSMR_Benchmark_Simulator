function deltau = control(control_law, t_s, y_output, y_sp)

switch control_law

   case 0 % PID control
      ku = 0.1265;
      tao = 0.25;
      %Ziegler-Nichols method
      kp = 18*(0.6*ku);
      ki = 0.01*((1.2*ku)/tao);
      kd = (3*ku*tao)/40;
      persistent integral_error
      persistent prev_error
      if isempty(integral_error) && isempty(prev_error) 
        integral_error = 0;  % Initialize only once
        prev_error = 0;  % Initialize only once
      end
      error = y_sp - y_output;
      integral_error = integral_error + error*t_s;
      derivative_error = (error - prev_error)/t_s;
      deltau = kp*error + ki*integral_error + kd*derivative_error;
      prev_error = error;
     
end


end

