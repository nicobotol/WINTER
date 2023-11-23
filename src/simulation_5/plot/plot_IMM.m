% plot the probabilities of each model and the estimated gain

if simulation.model == 5 || simulation.model==7

  t_start = out_store{1}.model_x_est.Time(end) - simulation.plot_time(1); % [s]
  [~, s_start] = min(abs(out_store{1}.model_x_est.Time - t_start)); % sample from where to start
  time = out_store{1}.model_x_est.Time(s_start:end);
                                     
  %    ___  _ __ ___   ___  __ _  __ _ 
  %   / _ \| '_ ` _ \ / _ \/ _` |/ _` |
  %  | (_) | | | | | |  __/ (_| | (_| |
  %   \___/|_| |_| |_|\___|\__, |\__,_|
  %                        |___/       
  for j = 1:2:wind.WS_len
    fig = figure('Color', 'w');hold on;box on;grid on;
    tz_start = 160; % [s] time where to start the zoom
    tz_stop = 160.5; % [s] time where to sop the zoom
    max_tmp = 0;  min_tmp = 1e9;
    for i=1:IMM.n_models
      if numel(size(out_store{j}.model_x_est.Data)) >= 3
        tmp = zeros(1, size(out_store{1}.model_x_est.Time,1));
        tmp(:,:) = out_store{j}.model_x_est.Data(1, i, :);
        tmp = tmp(s_start:end);
      else
        tmp = out_store{j}.model_x_est.Data(s_start:end, 1);
      end
      plot(time, tmp, 'LineWidth', line_width, 'DisplayName', ['Mod. ', num2str(i)], 'Color', color(i));
      max_tmp = max(max_tmp, max(tmp));
      min_tmp = min(min_tmp, min(tmp));
    end
    plot(time, out_store{j}.omega_tilde.Data(s_start:end), ':', 'LineWidth', line_width, 'DisplayName', 'Est. IMM', 'Color', color(IMM.n_models + 1));
    plot(time, out_store{j}.omega_R.Data(s_start:end), '--', 'LineWidth', line_width, 'DisplayName', 'Sim. IMM', 'Color', color(IMM.n_models + 2));
    plot(time, out_store{j + 1}.omega_R.Data(s_start:end), '--', 'LineWidth', line_width, 'DisplayName', 'Sim. fix gain', 'Color', color(IMM.n_models + 3));
    % yline(omega_rated, '--', 'LineWidth', line_width, 'DisplayName', 'Rated');
    ylabel('$\omega_R$ [rad/s]')
    xlabel('Time [s]')
    legend('Location', 'NorthEast', 'NumColumns', 2)
    title(['Comparison of $\omega_R$, $V_0 =$ ', num2str(wind.mean(j)), ' [m/s]'])
    ylim([0.66 0.80]) % for const WS 8 [m/s]
    % ylim([0.58 1.1]) 
   
    % zoom
    axes('position',[0.22 0.52 .30 .30])
    box on;hold on;grid on;
    %     plot(time_zoom, out_store{j + 1}.omega_R.Data(sz_start:sz_stop), '--', 'LineWidth', line_width, 'DisplayName', 'Sim. fix');
    [~, sz_start] = min(abs(out_store{j}.model_x_est.Time - tz_start)); % sample from where to start
    [~, sz_stop] = min(abs(out_store{j}.model_x_est.Time - tz_stop)); % sample from where to start
    for i=1:IMM.n_models
      if numel(size(out_store{j}.model_x_est.Data)) >= 3
        tmp = zeros(1, size(out_store{1}.model_x_est.Time,1));
        tmp(:,:) = out_store{j}.model_x_est.Data(1, i, :);
        tmp = tmp(sz_start:sz_stop);
      else
        tmp = out_store{j}.model_x_est.Data(sz_start:sz_stop, 1);
      end
      
      time_zoom = out_store{j}.model_x_est.Time(sz_start:sz_stop);
      plot(time_zoom, tmp, 'LineWidth', line_width, 'DisplayName', ['Mod. ', num2str(i)], 'Color', color(i));
    end
    plot(time_zoom, out_store{j}.omega_tilde.Data(sz_start:sz_stop), ':', 'LineWidth', line_width, 'DisplayName', 'Est. IMM', 'Color', color(IMM.n_models + 1));
    plot(time_zoom, out_store{j}.omega_R.Data(sz_start:sz_stop), '--', 'LineWidth', line_width, 'Color', color(IMM.n_models + 2));  
    xlim([tz_start tz_stop])
    % plot(time_zoom, out_store{j}.omega_tilde.Data(sz_start:sz_stop), ':', 'LineWidth', line_width);
    title('Zoom')
    set(gca, 'FontSize', font_size)
    if simulation.print_figure == 1
      export_figure(fig, strcat(date_fig, 'omega_IMM_',num2str(j), '.eps'), path_images);
    end
  end


%   _  __               _   
%  | |/ /    ___  _ __ | |_ 
%  | ' /    / _ \| '_ \| __|
%  | . \   | (_) | |_) | |_ 
%  |_|\_\___\___/| .__/ \__|
%      |_____|   |_|        
  for i = 1:wind.WS_len/2
    fig = figure('Color', 'w');grid on;box on;hold on;
    max_tmp = 0;  min_tmp = 1e9;
      plot(out_store{2*i - 1}.K_opt.Time(s_start:end), out_store{2*i - 1}.K_opt.Data(s_start:end), 'LineWidth',line_width, 'DisplayName', ['Sim. ', num2str(i)]);
      max_tmp = max(max_tmp, max(out_store{2*i - 1}.K_opt.Data(s_start:end)));
      min_tmp = min(min_tmp, min(out_store{2*i - 1}.K_opt.Data(s_start:end)));
    yline(generator.K_opt_GE,'-', 'LineWidth', 2*line_width, 'DisplayName', 'Reference', 'Color', [color(5), 1]);
    yline(IMM.K_vector, '--', 'LineWidth', 1.5*line_width, 'HandleVisibility', 'Off', 'Color', [color(2), 1]) ;
    plot(NaN, NaN, '--', 'LineWidth', 1.5*line_width, 'Color', [color(2), 1],  'DisplayName', 'IMM');
    xlabel('Time [s]')
    ylabel('$K_{opt} [Nms^2]$')
    legend('Location', 'SouthEast')
    % ylim([0.9*min_tmp 1.1*max_tmp])
    title(['$K_{opt}$, $V_0$ = ', num2str(wind.mean(2*i)), ' [m/s]']);
    set(gca, 'FontSize', font_size)
    if simulation.print_figure == 1
      export_figure(fig, strcat(date_fig, 'K_opt_IMM.eps'), path_images);
    end
  end

%                   _           _     _ _ _ _         
%   _ __  _ __ ___ | |__   __ _| |__ (_) (_) |_ _   _ 
%  | '_ \| '__/ _ \| '_ \ / _` | '_ \| | | | __| | | |
%  | |_) | | | (_) | |_) | (_| | |_) | | | | |_| |_| |
%  | .__/|_|  \___/|_.__/ \__,_|_.__/|_|_|_|\__|\__, |
%  |_|                                          |___/ 
for i=1:wind.WS_len/2
  fig = figure('Color', 'w');hold on;grid on; box on;
  data = reshape(out_store{2*i-1}.mu.Data(:,1,s_start:end), IMM.n_models, []);
  bar(time, data', 'stacked')
  % for i=1:IMM.n_models
  %   data = out_store{1}.mu.Data(i,1,s_start:end);
  %   data = reshape(data, 1, size(data,3));
  %   plot(time, data, 'LineWidth', line_width, 'DisplayName', ['Model ', num2str(i)]);
  % end
  xlabel('Time [s]')
  ylabel('$\mu$ [-]')
  legend('Model 1', 'Model 2', 'Model 3','Location', 'NorthEast')
  title(['Bar plot of the probabilities of each model, $V_0$ = ', num2str(wind.mean(2*i)), ' [m/s]'])
  if simulation.print_figure == 1
    export_figure(fig, strcat(date_fig, 'probability_bar_IMM.eps'), path_images);
  end


  fig = figure('Color', 'w');hold on;grid on; box on;
  data = reshape(out_store{2*i-1}.mu.Data(:,1,s_start:end), IMM.n_models, []);
  plot(time, data', 'LineWidth', 0.5*line_width);
  xlabel('Time [s]')
  ylabel('$\mu$ [-]')
  legend('Model 1', 'Model 2', 'Model 3','Location', 'Best')
  title(['Probability, $V_0$ = ', num2str(wind.mean(2*i)), ' [m/s]'])
  if simulation.print_figure == 1
    export_figure(fig, strcat(date_fig, 'probability_IMM.eps'), path_images);
  end
end


% figure(); hold on; box on; grid on;
% plot(out_store{1}.T_R*1e-6);
% plot(out_store{1}.T_R_no_noise*1e-6)
% legend('Eff', 'No noise')
% ylabel('Torque [MNm]')

% % rotor and generator output power dynamic
% plot_time_series2('fig_power_dynamic', out_store, 'P_R', 'P_G',  generator.P_rated, 'Time [s]', 'P [MW]', 'Rotor and generator powers', 1e6, 'Aero.', 'Gen.', date_fig, 'southeast')


%        _           
%   _ __| |__   ___  
%  | '__| '_ \ / _ \ 
%  | |  | | | | (_) |
%  |_|  |_| |_|\___/ 
                   
fig = figure('Color', 'w'); hold on; box on; grid on;
for i=1:2:wind.WS_len
  plot(out_store{i}.rho_eff.Time(s_start:end), out_store{i}.rho_eff.Data(s_start:end), 'DisplayName', 'Real', 'LineWidth', 0.5*line_width, 'Color', color(1))
  plot(out_store{i}.rho.Time(s_start:end), out_store{i}.rho.Data(s_start:end), 'DisplayName', 'Estimated', 'LineWidth', 0.5*line_width, 'Color', color(3))
  plot(out_store{i}.rho_filter.Time(s_start:end), out_store{i}.rho_filter.Data(s_start:end), 'DisplayName', 'Filtered', 'LineWidth', 0.75*line_width, 'Color', color(2))
  yline(IMM.rho_vector, '--', 'LineWidth', 1.2*line_width, 'Color', [color(7), 1], 'HandleVisibility', 'Off');
  plot(NaN, NaN, '--', 'LineWidth', 1.2*line_width, 'Color', color(7), 'DisplayName', 'IMM');
  ylim([0.95*min(IMM.rho_vector), 1.05*max(IMM.rho_vector)])
  legend('Location', 'SouthEast')
  xlabel('Time [s]')
  ylabel('$\rho [\frac{kg}{m^3}]$')
  title(['Estimated and real $\rho$, $V_0$ = ', num2str(wind.mean(i)), ' [m/s]'])
end
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'rho_IMM.eps'), path_images);
end

end