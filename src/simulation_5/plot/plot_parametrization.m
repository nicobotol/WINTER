function plot_parametrization(plot_name,out_cell,series_x,series_y, ...
  x_label,y_label,plot_title,scaling_x,scaling_y,date_fig)

parameters

leg = cell(1, wind.WS_len);
fig = figure('Color','w');
hold on
for i=1:wind.WS_len
  plot(out_cell{i}.(series_x).Data/scaling_x, out_cell{i}.(series_y).Data/scaling_y, ...
    'LineWidth', line_width, 'Color', colors_vect(i,:));
  leg{i} = ['Sim. ', num2str(i)];
end
% plot(x_my_ref, y_my_ref, '--', 'LineWidth', line_width, ...
%   'Color', colors_vect(i+1,:))
% leg{wind.WS_len + 1} =  ['Computed ref.'];
legend(leg, 'Location', 'best', 'FontSize', font_size,...
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