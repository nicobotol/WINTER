function [kp_schedule, ki_schedule] = run_stall_regulation(blade, i)
% This functions sets different blade pitch gain constant according what 
% is necessary to simulate: gain scheduling or stall regulation 

  if i == 1 % with gain scheduling
    kp_schedule = blade.kp_schedule;
    ki_schedule = blade.ki_schedule;
  elseif i == 2 % without scheduling
    kp_schedule = 0;
    ki_schedule = 0;
  end
end