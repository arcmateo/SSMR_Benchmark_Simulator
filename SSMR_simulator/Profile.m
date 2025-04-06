function y_sp = Profile(ss, time, t_s, type, Mode)

y_sp = ones(size(time))*ss;

switch type

   case 0
      y_sp = ones(size(time))*ss;

   case 1
      for k = 1:length(time)
         t = k*t_s;
         if t  > 10 && t <= 15
            y_sp(k) = ss*1.15;
         elseif t  > 20 && t <= 25
            y_sp(k) = ss*0.85;
         end
      end

   case 2
      if Mode == 1
         for k = 1:length(time)
         t = k*t_s;
         if t > 4 && t <= 6
            y_sp(k) = 1.7052e-5*t + 1.5915e-4;
         elseif t > 6 && t <= 9
            y_sp(k) = ss*1.15;
         elseif t > 9 && t <= 11
            y_sp(k) = -3.4103e-5*t + 5.6838e-4;
         elseif t > 11 && t <= 14
            y_sp(k) = ss*0.85;
         elseif t > 14 && t <= 17
            y_sp(k) = 2.6525e-5*t - 1.7810e-4;             
         elseif t > 17 && t <= 22
            y_sp(k) = ss*1.20;
         elseif t > 22 && t <= 25
            y_sp(k) = -1.5157e-5*t + 6.0628e-4; 
         end
        end
     end

      if Mode == 2
         for k = 1:length(time)
            t = k*t_s;
            if t > 4 && t <= 6
               y_sp(k) = 2.0729e-5*t + 1.9346e-4;
            elseif t > 6 && t <= 9
               y_sp(k) = ss*1.15;
            elseif t > 9 && t <= 11
               y_sp(k) = -4.1457e-5*t + 6.9095e-4;
            elseif t > 11 && t <= 14
               y_sp(k) = ss*0.85;
            elseif t > 14 && t <= 17
               y_sp(k) = 3.2244e-5*t - 2.1649e-4;
            elseif t > 17 && t <= 22
               y_sp(k) = ss*1.20;
            elseif t > 22 && t <= 25
               y_sp(k) = -1.8425e-5*t + 7.37e-4; 
            end
         end
      end

      if Mode == 3
         for k = 1:length(time)
            t = k*t_s;
            if t > 4 && t <= 6
               y_sp(k) = 4.3602e-5*t + 4.0695e-4;
            elseif t > 6 && t <= 9
               y_sp(k) = ss*1.15;
            elseif t > 9 && t <= 11
               y_sp(k) = -8.7204e-5*t + 1.45339655e-3;
            elseif t > 11 && t <= 14
               y_sp(k) = ss*0.85;
            elseif t > 14 && t <= 17
               y_sp(k) = 6.7825e-5*t - 4.5540e-4;
            elseif t > 17 && t <= 22
               y_sp(k) = ss*1.20;
            elseif t > 22 && t <= 25
               y_sp(k) = -3.8757e-5*t + 1.5502824e-3; 
            end
         end
      end

end

end
