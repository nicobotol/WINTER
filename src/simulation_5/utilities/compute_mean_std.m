function compute_mean_std(out_cell, simulation, IMM, wind)

  for i=1:wind.WS_len/2
    t_start = simulation.stop_time(2*i-1) - simulation.post_process_time(2*i-1); % time from when to start the integration
    [~, s_start] = min(abs(out_cell{2*i-1}.P_GE.Time - t_start)); % element where to start the integration
    data =  out_cell{2*i-1}.K_opt.Data(s_start:end);
    m_vector(i) = mean(data);
    s_vector(i) = std(data);
  end
  fileID = fopen('mean_std_IMM.txt','w');
  for i=1:wind.WS_len/2
    fprintf(fileID,'%.2f & %.3e & %.3e \\\\ \n', [IMM.sigma_gain(2*i), m_vector(i), s_vector(i)]');
  end

end