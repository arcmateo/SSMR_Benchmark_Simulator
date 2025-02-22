function y_sp = Profile(sp, time, t_s)

y_sp = ones(size(time))*sp;

for k = 1:length(time)

   t = k*t_s;

   if t  > 10 && t <= 15
      y_sp(k) = sp*1.05;
   elseif t  > 20 && t <= 25
      y_sp(k) = sp*0.9;
   end
end

end

