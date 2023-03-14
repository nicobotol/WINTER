function plot_static(plot_name, x_data, y_data, x_ref_DTU, y_ref_DTU, ...
  x_my_ref, y_my_ref, x_label, y_label, plot_title, scaling, date_fig)

parameters

leg = cell(1, wind.WS_len + 2);

fig = figure('Position', get(0, 'Screensize'), 'Color','w');
hold on
for i = 1:wind.WS_len
  plot(x_data(i), y_data(i)/scaling, 'o', 'LineWidth', line_width)
  leg{i} = ['Sim.', num2str(i)];
end
plot(x_ref_DTU, y_ref_DTU, 'LineWidth', line_width, ...
  'Color', colors_vect(i + 1,:),'LineStyle','--'); % DTU reference
plot(x_my_ref, y_my_ref, 'LineWidth',line_width, ...
  'Color', colors_vect(i + 2,:),'LineStyle','-.'); % ref. with my TSR
leg{wind.WS_len + 1} =  ['DTU 10MW ref.'];
leg{wind.WS_len + 2} =  ['Computed ref.'];
xlabel(x_label,'interpreter','latex')
ylabel(y_label,'interpreter','latex')
title(plot_title)
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.png');
  export_fig('fig', fig_name);
end