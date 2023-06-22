function plot_parametrization_windN(plot_name,out_cell,series,x_my_ref, ...
  y_my_ref,x_label,y_label,plot_title,scaling,simulation,leg,date_fig)

parameters

fig = figure('Color','w');
N = size(x_my_ref, 1); % data to be plotted
hold on
for j = 1:N
  plot(x_my_ref{j}(:), y_my_ref{j}(:), '--', 'LineWidth', line_width, ...
    'Color', colors_vect(wind.WS_len+j,:))
end
for i=1:wind.WS_len
  % Time from where start to print
  t_start = out_cell{i}.(series).Time(end) - simulation.plot_time(i); % [s]
  [~, s_start] = min(abs(out_cell{i}.(series).Time - t_start)); % sample from where to start
  
  % Resample the wind speed
  series_length = out_cell{i}.(series).TimeInfo.Length;
  wind_resampled = zeros(series_length, 1);
  wind_resampled = interp1(out_cell{i}.wind.Time,out_cell{i}.wind.Data, ...
    out_cell{i}.(series).Time);

  plot(wind_resampled(s_start:end), out_cell{i}.(series).Data(s_start:end)/scaling, ...
    'LineWidth', line_width, 'Color', colors_vect(i,:));
  leg{i + j} = ['Sim. ', num2str(i)];
end
legend(leg, 'Location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
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