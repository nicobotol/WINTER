% plot the probabilities of each model and the estimated gain

if simulation.model == 5

  t_start = out_store{1}.model_x_est.Time(end) - simulation.plot_time(1); % [s]
  [~, s_start] = min(abs(out_store{1}.model_x_est.Time - t_start)); % sample from where to start
  time = out_store{1}.model_x_est.Time(s_start:end);
  fig = figure('Color', 'w');hold on;box on;grid on;
  for i=1:IMM.n_models
    tmp = out_store{1}.model_x_est.Data(1,i,s_start:end);
    leg(i) = ['Mod. ', num2str(i)];
    plot(time, reshape(tmp, 1, size(tmp,3)), 'LineWidth', line_width, 'DisplayName', leg(i));
  end
  plot(time, out_store{1}.omega_tilde.Data(s_start:end), ':', 'LineWidth', line_width, 'DisplayName', 'Estimation');
  plot(time, out_store{1}.omega_R.Data(s_start:end), '--', 'LineWidth', line_width, 'DisplayName', 'Simulation');
  yline(omega_rated, '--', 'LineWidth', line_width, 'DisplayName', 'Rated');
  ylabel('$\omega_R$ [rad/s]')
  xlabel('Time [s]')
  legend()
  title('Comparison of $\omega_R$')
  ylim([0.9 1.2]*omega_rated)
  % zoom
  axes('position',[1/5 0.52 .30 .30])
  box on;hold on;grid on;
  tz_start = 30; % [s]
  [~, sz_start] = min(abs(out_store{1}.model_x_est.Time - tz_start)); % sample from where to start
  tz_stop = 30.5; % [s]
  [~, sz_stop] = min(abs(out_store{1}.model_x_est.Time - tz_stop)); % sample from where to start
  for i=1:IMM.n_models
    tmp = out_store{1}.model_x_est.Data(1,i,sz_start:sz_stop);
    time_zoom = out_store{1}.model_x_est.Time(sz_start:sz_stop);
    plot(time_zoom, reshape(tmp, 1, size(tmp,3)), 'LineWidth', line_width, 'DisplayName', leg(i));
  end
  xlim([30 30.5])
  plot(time_zoom, out_store{1}.omega_tilde.Data(sz_start:sz_stop), ':', 'LineWidth', line_width);
  plot(time_zoom, out_store{1}.omega_R.Data(sz_start:sz_stop), '--', 'LineWidth', line_width, 'DisplayName', 'Simulation');  
  title('Zoom')
  set(gca, 'FontSize', font_size)
  if simulation.print_figure == 1
    export_figure(fig, strcat(date_fig, 'omega_IMM.eps'), path_images);
  end


  fig = figure('Color', 'w');grid on;box on;hold on;
  plot(out_store{1}.K_opt.Time(s_start:end), out_store{1}.K_opt.Data(s_start:end), 'LineWidth',line_width, 'DisplayName', 'Estimation');
  yline(rotor.K_opt,'--', 'LineWidth',line_width, 'DisplayName', 'Reference');
  xlabel('Time [s]')
  ylabel('$\omega_R$ [rad/s]')
  legend('Location', 'SouthEast')
  ylim([0.95 1.03]*rotor.K_opt)
  title('$K_{opt}$')
  set(gca, 'FontSize', font_size)
  if simulation.print_figure == 1
    export_figure(fig, strcat(date_fig, 'omega_IMM.eps'), path_images);
  end


  fig = figure('Color', 'w');hold on;grid on;
  for i=1:IMM.n_models
    data = out_store{1}.mu.Data(i,1,s_start:end);
    data = reshape(data, 1, size(data,3));
    plot(time, data, 'LineWidth', line_width, 'DisplayName', ['Model ', num2str(i)]);
  end
  xlabel('Time [s]')
  ylabel('$\mu$ [-]')
  legend('Location', 'NorthEast')
  title('Probability')

end