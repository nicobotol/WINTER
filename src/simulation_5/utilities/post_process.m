function [RMS_errors] = post_process(out_cell, wind, omega_rated, generator, simulation, lookup_static_values, lookup_Pitch, lookup_P_GE, IMM, date_fig)
% This function performs the post processing opoerations of the obtained data

RMS_errors =  cell(1, wind.WS_len);

%% Errors on the obtained quantities wrt their expected values and their normalized values
if simulation.type == 5
  for i = 1:wind.WS_len
    % rotaional speed
    RMS_errors{i}.var{1}.Name = 'omega_R';
    rated = lookup_static_values(1:2,:); % wind speed and rotational speed
    [RMS_errors{i}.var{1}.Data] = compute_RMS_error(out_cell{i}, 'omega_R', simulation.post_process_time(i), rated);
    normalizer = interp1(lookup_static_values(1, :), lookup_static_values(2, :), wind.mean(i));
    RMS_errors{i}.var{1}.norm = RMS_errors{i}.var{1}.Data/normalizer*100;

    % generatator input power
    RMS_errors{i}.var{2}.Name = 'P_G';
    rated = [lookup_static_values(1,:); lookup_static_values(6,:)]; % wind speed and generator mechancial power
    [RMS_errors{i}.var{2}.Data] = compute_RMS_error(out_cell{i}, 'P_G', simulation.post_process_time(i), rated)/1e6;
    normalizer = interp1(lookup_static_values(1, :), lookup_static_values(6, :), wind.mean(i))/1e6;
    RMS_errors{i}.var{2}.norm = RMS_errors{i}.var{2}.Data/normalizer*100;
    
    % generatator input torque
    RMS_errors{i}.var{3}.Name = 'T_G';
    rated = [lookup_static_values(1,:); lookup_static_values(6,:)./lookup_static_values(2,:)]; % wind speed and generator mechancial power
    [RMS_errors{i}.var{3}.Data] = compute_RMS_error(out_cell{i}, 'T_G', simulation.post_process_time(i), rated)/1e6;
    normalizer = interp1(lookup_static_values(1, :), lookup_static_values(6, :), wind.mean(i)) / interp1(lookup_static_values(1, :), lookup_static_values(2, :), wind.mean(i))/1e6;
    RMS_errors{i}.var{3}.norm = RMS_errors{i}.var{3}.Data/normalizer*100;

    % pitch angle
    RMS_errors{i}.var{4}.Name = 'pitch';
    rated = [lookup_Pitch(1, :); lookup_Pitch(3, :)];
    [RMS_errors{i}.var{4}.Data] = compute_RMS_error(out_cell{i}, 'pitch', simulation.post_process_time(i), rated);
    normalizer = interp1(lookup_Pitch(1, :), lookup_Pitch(3, :), wind.mean(i));
    RMS_errors{i}.var{4}.norm = RMS_errors{i}.var{4}.Data/normalizer*100;
  end
elseif simulation.type == 10
  
  compute_energy_develop(out_cell, simulation, IMM, wind, simulation.stop_time(1)*[1], 'energy_K_opt_comp', date_fig);

elseif simulation.type == 11
  
  compute_comp_gains(out_cell, simulation, IMM, wind, simulation.stop_time(1)*[0.25 0.5 0.75 1], 'energy_com_const_gains', generator);

elseif simulation.type == 12
    % compute the mean and std value of the gain for all the cases
    compute_mean_std_develop(out_cell, simulation, IMM, wind, generator, date_fig);

    % compute the energy 
    compute_energy_develop(out_cell, simulation, IMM, wind, simulation.stop_time(1)*[1], 'energy', date_fig);
else
  for i = 1:wind.WS_len
    % rotaional speed
    [RMS_errors{i}.omega_R] = compute_RMS_error(out_cell{i}, 'omega_R', simulation.post_process_time(i), omega_rated);
    RMS_errors{i}.omega_R_norm = RMS_errors{i}.omega_R/omega_rated*100;

    % generatator input power
    [RMS_errors{i}.P_G] = compute_RMS_error(out_cell{i}, 'P_G', simulation.post_process_time(i), generator.P_rated);
    RMS_errors{i}.P_G_norm = RMS_errors{i}.P_G/generator.P_rated*100;
    
    % generatator input torque
    [RMS_errors{i}.T_G] = compute_RMS_error(out_cell{i}, 'T_G', simulation.post_process_time(i), generator.P_rated/generator.omega_rated);
    RMS_errors{i}.T_G_norm = RMS_errors{i}.T_G/(generator.P_rated/generator.omega_rated)*100;

    % generatator electrical output power
    rated = [lookup_P_GE(1,:); lookup_P_GE(3,:)];
    [RMS_errors{i}.P_GE] = compute_RMS_error(out_cell{i}, 'P_GE', simulation.post_process_time(i), rated);
    RMS_errors{i}.P_GE_norm = RMS_errors{i}.P_GE/lookup_P_GE(3,end)*100;
  end
end

%% Write the RMS table to a file
if simulation.type == 5
  mode = 'a';
  write_RMS_table(RMS_errors,mode);
else
  mode = 'w';
end


end

%% RMS errors
function [RMS] = compute_RMS_error(cell, series, pp_time, ref_value)
  t_start = cell.(series).Time(end) - pp_time; % [s]
  [~, s_start] = min(abs(cell.(series).Time - t_start)); % sample from where to start
  data = cell.(series).Data(s_start:end);
  if length(ref_value) == 1 % case of constant reference
    data_setpoint = ones(length(data), 1)*ref_value;
  else % case of time varying reference
    
    V0 = cell.wind.Data(s_start:end); % wind speed
    data_setpoint = interp1(ref_value(1,:), ref_value(2,:), V0); % interpolate the reference value
  end
  RMS = RMS_error(data,data_setpoint);
end

function write_RMS_table(RMS, mode)
  % This function writes the RMS table to a file

  % Open the file
  fid = fopen('RMS_table.txt', mode);
  for i=1:length(RMS)
    format = [];
    data = [];
    for j=1:length(RMS{i}.var)
      format = [format, '%5.2f & %5.2f &'];
      data = [data, RMS{i}.var{j}.Data, RMS{i}.var{j}.norm];
    end 
    format(end) = '';
    format = [format, '\\\\ \n'];
    data(isnan(data)) = 0.0;
    data(isinf(data)) = 0.0;
    fprintf(fid, format, data);
  end
  fclose(fid);

end