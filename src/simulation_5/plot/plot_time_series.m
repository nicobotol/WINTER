function plot_time_series(plot_name,out_cell,series,x_label,y_label,...
  plot_title,scaling,date_fig)

parameters

leg = cell(1, wind.WS_len);
fig = figure('Color','w');
hold on
for i = 1:wind.WS_len
  % Time from where start to print
  t_start = out_cell{i}.(series).Time(end) - simulation.plot_time(i); % [s]
  [~, s_start] = min(abs(out_cell{i}.(series).Time - t_start)); % sample from where to start

  plot(out_cell{i}.(series).Time(s_start:end), out_cell{i}.(series).Data(s_start:end)/scaling, ...
    'LineWidth', line_width);
  leg{i} = ['Sim. ', num2str(i)];
end
legend(leg, 'Location', 'northwest', 'FontSize', font_size,...
  'interpreter','latex');
xlabel(x_label,'interpreter','latex')
ylabel(y_label,'interpreter','latex')
title(plot_title, 'Interpreter','latex')
grid on
box on
set(gca, 'FontSize', font_size)
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
end