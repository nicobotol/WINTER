function plot_parametrization_wind_P_GEN(plot_name,out_cell,series,x_my_ref, y_my_ref,x_label,y_label,plot_title,scaling,simulation,leg,date_fig)

parameters

fig = figure('Color','w');
N = size(x_my_ref, 1); % data to be plotted
hold on
color_i = [20 3];
for j = 1:N
  plot(x_my_ref{j}(:), y_my_ref{j}(:), '--', 'LineWidth', line_width, 'Color', color(color_i(j)))
end
leg{N + 1} = ['Sim. Generator'];
leg{N + 2} = ['Sim. Rotor'];
plot(NaN, NaN, '-', 'LineWidth', line_width, 'Color', color(1), 'MarkerSize', marker_size)
plot(NaN, NaN, '-', 'LineWidth', line_width, 'Color', color(2), 'MarkerSize', marker_size)
for i=1:wind.WS_len
  % Time from where start to print
  t_start = out_cell{i}.(series).Time(end) - simulation.plot_time(i); % [s]
  [~, s_start] = min(abs(out_cell{i}.(series).Time - t_start)); % sample from where to start
  wind_resampled = out_cell{i}.wind.Data;

    if rem(i, 2)==1
      sign = '-';
      col = color(1);
    else
      sign = '-'; 
      col = color(2);
    end
      plot(wind_resampled(s_start:end), out_cell{i}.(series).Data(s_start:end)/scaling, sign,'LineWidth', line_width, 'Color', col, 'HandleVisibility','off', 'MarkerSize', marker_size);
end
legend(leg, 'Location', 'best', 'FontSize', font_size, 'interpreter','latex');
xlabel(x_label,'interpreter','latex')
ylabel(y_label,'interpreter','latex')
title(plot_title, 'Interpreter','latex')
grid on
set(gca, 'FontSize', font_size)
box on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
end