function plot_time_series(plot_name,out_cell,series,x_label,y_label,...
  plot_title,scaling,date)

parameters

leg = cell(1, wind.WS_len);
fig = figure('Position', get(0, 'Screensize'), 'Color','w');
hold on
for i = 1:wind.WS_len
  plot(out_cell{i}.(series).Time, out_cell{i}.(series).Data/scaling, ...
    'LineWidth', line_width);
  leg{i} = ['Sim. ', num2str(i)];
end
legend(leg, 'Location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
xlabel(x_label,'interpreter','latex')
ylabel(y_label,'interpreter','latex')
title(plot_title, 'Interpreter','latex')
grid on
set(gca, 'FontSize', font_size)
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, plot_name,'.png');
  export_fig('fig', fig_name);
end