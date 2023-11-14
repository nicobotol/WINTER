function energy = compute_comp_gains(out_store, simulation, IMM, wind, vector, name, generator)
  % This function computes the integral of the power in a given time
  fileID = fopen(['../../report/macro/',name,'.tex'],'w');

  for j=1:length(vector)
    simulation.post_process_time = vector(j)*ones(length(simulation.post_process_time),1);
    energy_R = zeros(wind.WS_len, 1);
    diff_energy_R = zeros(wind.WS_len, 1);
    energy_GE = zeros(wind.WS_len, 1);
    diff_energy_GE = zeros(wind.WS_len, 1);
    
    % energy from the simulation with the reference k_opt_GE
    t_start = simulation.stop_time(4) - simulation.post_process_time(4);
    [~, s_start] = min(abs(out_store{4}.P_GE.Time - t_start)); % element where to start the integration
    energy_ref_GE = trapz(out_store{4}.P_GE.Time(s_start:end), out_store{4}.P_GE.Data(s_start:end));
    energy_ref_R = trapz(out_store{4}.P_R.Time(s_start:end), out_store{4}.P_R.Data(s_start:end));
    
    for i=1:wind.WS_len
      t_start = simulation.stop_time(i) - simulation.post_process_time(i);
      [~, s_start] = min(abs(out_store{i}.P_GE.Time - t_start)); % element where to start the integration
      energy_GE(i) = trapz(out_store{i}.P_GE.Time(s_start:end), out_store{i}.P_GE.Data(s_start:end));
      energy_R(i) = trapz(out_store{i}.P_R.Time(s_start:end), out_store{i}.P_R.Data(s_start:end));
      diff_energy_R(i) = -(energy_R(i) - energy_ref_R)./energy_ref_R*100;
      diff_energy_GE(i) = -(energy_GE(i) - energy_ref_GE)./energy_ref_GE*100;
    end

    % write a latex table on a .txt file. In the first column there should be the value insiede IMM.sigma_gain, in the second the odd values in energy, while in the third the even values in energy
    fprintf(fileID, '\\multirow{%d}{*}{%.2f}', length(wind.WS_len), vector(j));
    for i=1:wind.WS_len
      fprintf(fileID,' & %d & %.2f & %.3f & %.3f & %.3f & %.3f\\\\ \n', [i, generator.K_opt_sensitivity(i), energy_R(i)/1e9, diff_energy_R(i), energy_GE(i)/1e9, diff_energy_GE(i)]');
    end
    fprintf(fileID, '\\midrule');
    fprintf(fileID, '\n \n');

  end
  fclose(fileID);
end