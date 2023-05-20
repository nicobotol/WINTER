function [kp_schedule, ki_schedule] = run_blade_gains(blade, i)
 % This functions sets different blade pitch gain constant according what 
 % is necessary to simulate 

  if i == 1 % with gain scheduling (don't do anything)
    kp_schedule = blade.kp_schedule;
    ki_schedule = blade.ki_schedule;
  elseif i == 2 % without scheduling (max value)
    kp_schedule = 2;
    ki_schedule = 0.9;
  elseif i == 3 % without scheduling (min value)
    kp_schedule = 0.25;
    ki_schedule = 0.2;
  elseif i == 4 % Olimpo's gain scheduling
    kp_schedule = blade.kp_schedule_report; % from Olimpo's
    ki_schedule = blade.ki_schedule_report; % from Olimpo's
  end