% plot the probabilities of each model and the estimated gain

if simulation.model == 5

  t_start = out_store{1}.model_x_est.Time(end) - simulation.plot_time(1); % [s]
  [~, s_start] = min(abs(out_store{1}.model_x_est.Time - t_start)); % sample from where to start
  time = out_store{1}.model_x_est.Time(s_start:end);
  for j = 1:2:wind.WS_len
    fig = figure('Color', 'w');hold on;box on;grid on;
    max_tmp = 0;  min_tmp = 1e9;
    for i=1:IMM.n_models
      if numel(size(out_store{j}.model_x_est.Data)) >= 3
        tmp = zeros(1, size(out_store{1}.model_x_est.Time,1));
        tmp(:,:) = out_store{j}.model_x_est.Data(1, i, :);
        tmp = tmp(s_start:end);
      else
        tmp = out_store{j}.model_x_est.Data(s_start:end, 1);
      end
      plot(time, tmp, 'LineWidth', line_width, 'DisplayName', ['Mod. ', num2str(i)]);
      max_tmp = max(max_tmp, max(tmp));
      min_tmp = min(min_tmp, min(tmp));
    end
    plot(time, out_store{j}.omega_tilde.Data(s_start:end), ':', 'LineWidth', line_width, 'DisplayName', 'Est. IMM');
    plot(time, out_store{j}.omega_R.Data(s_start:end), '--', 'LineWidth', line_width, 'DisplayName', 'Sim. IMM');
    plot(time, out_store{j + 1}.omega_R.Data(s_start:end), '--', 'LineWidth', line_width, 'DisplayName', 'Sim. fix');
    % yline(omega_rated, '--', 'LineWidth', line_width, 'DisplayName', 'Rated');
    ylabel('$\omega_R$ [rad/s]')
    xlabel('Time [s]')
    legend()
    title(['Comparison of $\omega_R$ - $V_0 =$ ', num2str(wind.mean(j)), ' [m/s]'])
    ylim([min_tmp*0.9 max_tmp*1.5])
    % zoom
%     axes('position',[1/5 0.52 .30 .30])
%     box on;hold on;grid on;
%     tz_start = 160; % [s]
%     [~, sz_start] = min(abs(out_store{j}.model_x_est.Time - tz_start)); % sample from where to start
%     tz_stop = 160.5; % [s]
%     [~, sz_stop] = min(abs(out_store{j}.model_x_est.Time - tz_stop)); % sample from where to start
%     for i=1:IMM.n_models
%       if numel(size(out_store{j}.model_x_est.Data)) >= 3
%         tmp = zeros(1, size(out_store{1}.model_x_est.Time,1));
%         tmp(:,:) = out_store{j}.model_x_est.Data(1, i, :);
%         tmp = tmp(sz_start:sz_stop);
%       else
%         tmp = out_store{j}.model_x_est.Data(sz_start:sz_stop, 1);
%       end
%       
%       time_zoom = out_store{j}.model_x_est.Time(sz_start:sz_stop);
%       plot(time_zoom, tmp, 'LineWidth', line_width, 'DisplayName', ['Mod. ', num2str(i)]);
%     end
%     xlim([tz_start tz_stop])
%     plot(time_zoom, out_store{j}.omega_tilde.Data(sz_start:sz_stop), ':', 'LineWidth', line_width);
%     plot(time_zoom, out_store{j}.omega_R.Data(sz_start:sz_stop), '--', 'LineWidth', line_width, 'DisplayName', 'Simulation');  
%     plot(time_zoom, out_store{j + 1}.omega_R.Data(sz_start:sz_stop), '--', 'LineWidth', line_width, 'DisplayName', 'Sim. fix');
%     title('Zoom')
    set(gca, 'FontSize', font_size)
    if simulation.print_figure == 1
      export_figure(fig, strcat(date_fig, 'omega_IMM_',num2str(j), '.eps'), path_images);
    end
  end


  fig = figure('Color', 'w');grid on;box on;hold on;
  max_tmp = 0;  min_tmp = 1e9;
  for i = 1:wind.WS_len/2
    plot(out_store{2*i - 1}.K_opt.Time(s_start:end), out_store{2*i - 1}.K_opt.Data(s_start:end), 'LineWidth',line_width, 'DisplayName', ['Sim. ', num2str(i)]);
    max_tmp = max(max_tmp, max(out_store{2*i - 1}.K_opt.Data(s_start:end)));
    min_tmp = min(min_tmp, min(out_store{2*i - 1}.K_opt.Data(s_start:end)));
  end
  yline(generator.K_opt_GE,'--', 'LineWidth',line_width, 'DisplayName', 'Reference');
  xlabel('Time [s]')
  ylabel('$K_{opt} [Nms^2]$')
  legend('Location', 'SouthEast')
  % ylim([0.9*min_tmp 1.1*max_tmp])
  title('$K_{opt}$')
  set(gca, 'FontSize', font_size)
  if simulation.print_figure == 1
    export_figure(fig, strcat(date_fig, 'omega_IMM.eps'), path_images);
  end

%                   _           _     _ _ _ _         
%   _ __  _ __ ___ | |__   __ _| |__ (_) (_) |_ _   _ 
%  | '_ \| '__/ _ \| '_ \ / _` | '_ \| | | | __| | | |
%  | |_) | | | (_) | |_) | (_| | |_) | | | | |_| |_| |
%  | .__/|_|  \___/|_.__/ \__,_|_.__/|_|_|_|\__|\__, |
%  |_|                                          |___/ 
  for i=1:wind.WS_len
    fig = figure('Color', 'w');hold on;grid on;
    data = reshape(out_store{i}.mu.Data(:,1,s_start:end), IMM.n_models, []);
    bar(time, data', 'stacked')
    % for i=1:IMM.n_models
    %   data = out_store{1}.mu.Data(i,1,s_start:end);
    %   data = reshape(data, 1, size(data,3));
    %   plot(time, data, 'LineWidth', line_width, 'DisplayName', ['Model ', num2str(i)]);
    % end
    xlabel('Time [s]')
    ylabel('$\mu$ [-]')
    legend('Location', 'NorthEast')
    title('Probability')
  end


figure(); hold on;
plot(out_store{1}.T_R*1e-6);
plot(out_store{1}.T_R_no_noise*1e-6)
legend('Eff', 'No noise')
ylabel('Torque [MNm]')

% rotor and generator output power dynamic
plot_time_series2('fig_power_dynamic', out_store, 'P_R', 'P_G',  generator.P_rated, 'Time [s]', 'P [MW]', 'Rotor and generator powers', 1e6, 'Aero.', 'Gen.', date_fig, 'southeast')


figure(); hold on;
for i=1:wind.WS_len
  plot(out_store{i}.K_opt, 'DisplayName', num2str(i))
end
legend()


end