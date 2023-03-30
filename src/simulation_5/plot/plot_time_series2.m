function plot_time_series2(plot_name,out_cell,series1,series2,y_line, ...
  x_label,y_label,plot_title,scaling,leg1,leg2,date_fig)

parameters

leg = cell(2*wind.WS_len, 1);
fig = figure('Color','w');
hold on
for i = 1:wind.WS_len
  plot(out_cell{i}.(series1).Time, out_cell{i}.(series1).Data/scaling, ...
    'LineWidth',line_width*0.6,'Color',colors_vect(i, :),'LineStyle','--');
  plot(out_cell{i}.(series2).Time, out_cell{i}.(series2).Data/scaling, ...
    'LineWidth', line_width*1.2, 'Color', colors_vect(i, :));
  leg{2*i - 1} = [leg1,' sim. ', num2str(i)];
  leg{2*i} = [leg2,' sim. ', num2str(i)];
end
if strcmp(y_line, 'none') == 0
  yline(y_line/scaling, 'LineStyle', '-.', 'LineWidth', line_width,...
    'Color', colors_vect(i + 1, :));
  leg{2*wind.WS_len + 1} = ['Rated'];
end
legend(leg, 'Location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
xlabel(x_label,'interpreter','latex')
ylabel(y_label,'interpreter','latex')
title(plot_title, 'Interpreter','latex')
grid on
set(gca, 'FontSize', font_size)
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.png');
  export_fig('fig', fig_name);
end