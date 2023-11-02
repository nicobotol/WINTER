function energy = compute_energy_K_opt_sensitivity(out_store, simulation, IMM, wind)
  % This function computes the integral of the power in a given time
  fileID = fopen('../../report/macro/energy_sensitivity.tex','w');
  vector = [100 150 200 250];

  for j=1:length(vector)
    simulation.post_process_time = vector(j)*ones(length(simulation.post_process_time),1);
    energy = zeros(wind.WS_len, 1);
    for i=1:wind.WS_len
      t_start = simulation.stop_time(i) - simulation.post_process_time(i);
      [~, s_start] = min(abs(out_store{i}.P_GE.Time - t_start)); % element where to start the integration
      energy_GE(i) = trapz(out_store{i}.P_GE.Time(s_start:end), out_store{i}.P_GE.Data(s_start:end));
      energy_R(i) = trapz(out_store{i}.P_R.Time(s_start:end), out_store{i}.P_R.Data(s_start:end));
      diff_energy_R(i) = (energy_R(1) - energy_R(i))./energy_R(1)*100;
      diff_energy_GE(i) = (energy_GE(1) - energy_GE(i))./energy_GE(1)*100;
    end

    % write a latex table on a .txt file. In the first column there should be the value insiede IMM.sigma_gain, in the second the odd values in energy, while in the third the even values in energy
    % fprintf(fileID, '\\multirow{%d}{*}{%.2f}', length(wind.mean)/2, vector(j));
    idx =  {"IMM", "7178000", "9117000", "11007427", "11055000", "12994000", "14932000"};
    for i=1:wind.WS_len
      a  = idx(i);
      fprintf(fileID,'%s & %.2f & %.3f & %.2f & %.3f & %.2f\\\\ \n', [idx{i}, wind.mean(i), energy_R(i)/1e9, diff_energy_R(i), energy_GE(i)/1e9, diff_energy_GE(i)]');
    end
    fprintf(fileID, '\\midrule');
    fprintf(fileID, '\n \n');

  end
  fclose(fileID);
end