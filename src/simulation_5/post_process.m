function [RMS_errors] = post_process(out_cell, wind, omega_rated, generator, simulation)
% This function performs the post processing opoerations of the obtained data

RMS_errors =  cell(1, wind.WS_len);

%% Errors on the obtained quantities wrt their expected values
for i = 1:wind.WS_len
  % rotaional speed
  [RMS_errors{i}.omega_R] = compute_RMS_error(out_cell{i}.omega_R, simulation.post_process_time(i), omega_rated);

  % rotational speed
  % t_start = out_cell{i}.omega_R.Time(end) - simulation.post_process_time(i); % [s]
  % [~, s_start] = min(abs(out_cell{i}.omega_R.Time - t_start)); % sample from where to start
  % omega_R = out_cell{i}.omega_R.Data(s_start:end);
  % omega_R_set_point = ones(length(omega_R), 1)*omega_rated;
  % RMS_errors{i}.omega_R = RMS_error(omega_R, omega_R_set_point);
  
  % generatator input power
  [RMS_errors{i}.P_G] = compute_RMS_error(out_cell{i}.P_G, simulation.post_process_time(i), generator.P_rated);
  % t_start = out_cell{i}.P_G.Time(end) - simulation.post_process_time(i); % [s]
  % [~, s_start] = min(abs(out_cell{i}.P_G.Time - t_start)); % sample from where to start
  % P_G = out_cell{i}.P_G.Data(s_start:end);
  % P_G_set_point = ones(length(P_G), 1)*rotor.P_rated;
  % RMS_errors{i}.P_G = RMS_error(P_G, P_G_set_point);
  
  % generatator input torque
  [RMS_errors{i}.T_G] = compute_RMS_error(out_cell{i}.T_G, simulation.post_process_time(i), generator.P_rated/generator.omega_rated);
end

%% RMS errors
function [RMS] = compute_RMS_error(data_set, pp_time, ref_value)
  t_start = data_set.Time(end) - pp_time; % [s]
  [~, s_start] = min(abs(data_set.Time - t_start)); % sample from where to start
  data = data_set.Data(s_start:end);
  data_setpoint = ones(length(data), 1)*ref_value;
  RMS = RMS_error(data,data_setpoint);
end



end