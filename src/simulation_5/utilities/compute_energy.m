function compute_energy(out_store, simulation, IMM, wind)
  % This function computes the integral of the power in a given time
  energy = zeros(wind.WS_len, 1);

  for i=1:wind.WS_len
    t_start = simulation.stop_time(i) - simulation.post_process_time(i); % time from when to start the integration
    [~, s_start] = min(abs(out_store{i}.P_GE.Time - t_start)); % element where to start the integration
    energy(i) = trapz(out_store{i}.P_GE.Time(s_start:end), out_store{i}.P_GE.Data(s_start:end));
  end
  
  % write a latex table on a .txt file. In the first column there should be the value insiede IMM.sigma_gain, in the second the odd values in energy, while in the third the even values in energy
  fileID = fopen('energy.txt','w');
  for i=1:wind.WS_len/2
    fprintf(fileID,'%.2f & %.6e & %.6e \\\\ \n', [wind.mean(2*i), energy(2*i-1), energy(2*i)]');
  end
  fclose(fileID);

end