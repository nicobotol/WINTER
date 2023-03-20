function plot_parametrization_wind(plot_name,out_cell,series,x_my_ref, ...
  y_my_ref,x_label,y_label,plot_title,scaling,date_fig)

parameters

leg = cell(1, wind.WS_len + 1);
fig = figure('Position', get(0, 'Screensize'), 'Color','w');
hold on
for i=1:wind.WS_len
  % Resample the wind speed
  series_length = out_cell{i}.(series).TimeInfo.Length;
  wind_resampled = zeros(series_length, 1);
  wind_resampled = interp1(out_cell{i}.wind.Time,out_cell{i}.wind.Data, ...
    out_cell{i}.(series).Time);

  plot(wind_resampled, out_cell{i}.(series).Data/scaling, ...
    'LineWidth', line_width, 'Color', colors_vect(i,:));
  leg{i} = ['Sim. ', num2str(i)];
end
plot(x_my_ref, y_my_ref, '--', 'LineWidth', line_width, ...
  'Color', colors_vect(i+1,:))
leg{wind.WS_len + 1} =  ['Computed ref.'];
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