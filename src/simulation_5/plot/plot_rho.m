fig = figure('Color','w'); hold on; grid on; box on;
plot(time, out_store{1}.rho.Data(s_start:end), 'LineWidth', 0.5*line_width, 'LineStyle', '-', 'Color', [0 0 0 0.5], 'DisplayName', 'True')
for i=1:2:wind.WS_len
  t_start = out_store{1}.model_x_est.Time(end) - simulation.plot_time(1); % [s]
  [~, s_start] = min(abs(out_store{1}.model_x_est.Time - t_start)); % sample from where to start
  time = out_store{1}.model_x_est.Time(s_start:end);

  plot(time, out_store{i}.rho_filter.Data(s_start:end), 'LineWidth', line_width, 'DisplayName', ['$V_0=', num2str(wind.mean(i), 3), ' [\frac{m}{s}$]']);
end
plot(NaN, NaN, 'LineWidth', 1.5*line_width, 'LineStyle', '--', 'Color', [color(9), 1], 'DisplayName', 'True')
yline(IMM.rho_vector, 'LineWidth', 1.5*line_width, 'LineStyle', '--', 'Color', [color(9), 1], 'HandleVisibility', 'Off')
xlabel('Time [s]');
ylabel('$\rho [\frac{kg}{m^3}]$');
title('Estimated and real $\rho$')
legend('Location', 'best', 'NumColumns', 4);
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'rho_comparison.eps'), path_images);
end