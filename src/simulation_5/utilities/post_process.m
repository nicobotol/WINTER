function [RMS_errors] = post_process(out_cell, wind, omega_rated, generator, simulation)
% This function performs the post processing opoerations of the obtained data

RMS_errors =  cell(1, wind.WS_len);

%% Errors on the obtained quantities wrt their expected values and their normalized values
for i = 1:wind.WS_len
  % rotaional speed
  [RMS_errors{i}.omega_R] = compute_RMS_error(out_cell{i}.omega_R, simulation.post_process_time(i), omega_rated);
  RMS_errors{i}.omega_R_norm = RMS_errors{i}.omega_R/omega_rated*100;

  % generatator input power
  [RMS_errors{i}.P_G] = compute_RMS_error(out_cell{i}.P_G, simulation.post_process_time(i), generator.P_rated);
  RMS_errors{i}.P_G_norm = RMS_errors{i}.P_G/generator.P_rated*100;
  
  % generatator input torque
  [RMS_errors{i}.T_G] = compute_RMS_error(out_cell{i}.T_G, simulation.post_process_time(i), generator.P_rated/generator.omega_rated);
  RMS_errors{i}.T_G_norm = RMS_errors{i}.T_G/(generator.P_rated/generator.omega_rated)*100;
end

%% Write the RMS table to a file
write_RMS_table(RMS_errors);

end

%% RMS errors
function [RMS] = compute_RMS_error(data_set, pp_time, ref_value)
  t_start = data_set.Time(end) - pp_time; % [s]
  [~, s_start] = min(abs(data_set.Time - t_start)); % sample from where to start
  data = data_set.Data(s_start:end);
  data_setpoint = ones(length(data), 1)*ref_value;
  RMS = RMS_error(data,data_setpoint);
end

function write_RMS_table(RMS)
  % This function writes the RMS table to a file

  % Open the file
  fid = fopen('RMS_table.txt', 'w');
  for i=1:length(RMS)
    fprintf(fid, '%5.2f & %5.2f & %5.2f  & %5.2f  & %5.2f  & %5.2f \n', RMS{i}.omega_R, RMS{i}.omega_R_norm, RMS{i}.P_G/1e6, RMS{i}.P_G_norm, RMS{i}.T_G/1e6, RMS{i}.T_G_norm);
  end
  fclose(fid);

end