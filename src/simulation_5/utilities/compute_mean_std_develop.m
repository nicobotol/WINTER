function compute_mean_std_develop(out_cell, simulation, IMM, wind, generator, date_fig)
  fileID = fopen(strcat('../../report/macro/',date_fig,'mean_std_IMM.tex'),'w'); % file to write
  vector = [250]; % different time length to average

  for j=1:length(vector) % loop on every time to investigate
    simulation.post_process_time = vector(j)*ones(length(simulation.post_process_time),1);

    for i=1:wind.WS_len/2 % loop on every windspeed
      t_start = simulation.stop_time(2*i-1) - simulation.post_process_time(2*i-1); % time from when to start the integration
      [~, s_start] = min(abs(out_cell{2*i-1}.P_GE.Time - t_start)); % element where to start the integration
      data =  out_cell{2*i-1}.K_opt.Data(s_start:end);
      m_vector(i) = mean(data);
      s_vector(i) = std(data);
    end
    m_error = (m_vector - generator.K_opt_GE)/generator.K_opt_GE*100;
    fprintf(fileID, '\\multirow{%d}{*}{%.2f}', length(wind.mean)/2, vector(j));
    for i=1:wind.WS_len/2
      fprintf(fileID,'& %.2f & %.3e & %.3e & %.2f \\\\ \n', [wind.mean(2*i), m_vector(i), s_vector(i), m_error(i)]');
    end
    fprintf(fileID, '\\midrule');
    fprintf(fileID, '\n \n');
  end

  fclose(fileID);
end