function energy = compute_energy_develop(out_store, simulation, IMM, wind, vector, name, date_fig)
  % This function computes the integral of the power in a given time
  fileID = fopen(strcat('../../report/macro/',date_fig, name,'.tex'),'w');

  for j=1:length(vector)
    simulation.post_process_time = vector(j)*ones(length(simulation.post_process_time),1);
    energy = zeros(wind.WS_len, 1);
    for i=1:wind.WS_len
      t_start = simulation.stop_time(i) - simulation.post_process_time(i);
      [~, s_start] = min(abs(out_store{i}.P_GE.Time - t_start)); % element where to start the integration
      energy_GE(i) = trapz(out_store{i}.P_GE.Time(s_start:end), out_store{i}.P_GE.Data(s_start:end));
      energy_R(i) = trapz(out_store{i}.P_R.Time(s_start:end), out_store{i}.P_R.Data(s_start:end));
    end
    for i=1:wind.WS_len/2
      diff_energy(2*i - 1) = -(energy_R(2*i) - energy_R(2*i - 1))./energy_R(2*i - 1)*100;
      diff_energy(2*i) = -(energy_GE(2*i) - energy_GE(2*i - 1))./energy_GE(2*i - 1)*100;
    end
    % write a latex table on a .txt file. In the first column there should be the value insiede IMM.sigma_gain, in the second the odd values in energy, while in the third the even values in energy
   
    fprintf(fileID, '\\multirow{%d}{*}{%.2f}', length(wind.mean)/2, vector(j));
    for i=1:wind.WS_len/2
    fprintf(fileID,' & %.2f & %.4f & %.4f & %.2f & %.4f & %.4f & %.2f\\\\ \n', [wind.mean(2*i), energy_R(2*i-1)/1e9, energy_R(2*i)/1e9, diff_energy(2*i-1), energy_GE(2*i-1)/1e9, energy_GE(2*i)/1e9, diff_energy(2*i)]');
    end
    fprintf(fileID, '\\midrule');
    fprintf(fileID, '\n \n');
    
  end
  fclose(fileID);
end