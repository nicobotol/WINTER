function plot_parametrization_wind_zoom_P_GE(plot_name,out_cell,series,x_my_ref, ...
  y_my_ref,x_label,y_label,plot_title,scaling,simulation,date_fig)

  parameters
  
  min_zoom_V0 = 11.13;
  max_zoom_V0 = 11.15;
  axis_pos_x = 2.8/5;
  axis_pos_y = 0.4/1.2;
  axis_scale_x = 0.3;
  axis_scale_y = 0.3;

  leg = cell(1, wind.WS_len + 1);
  fig = figure('Color','w');
  hold on
  plot(x_my_ref, y_my_ref, '--', 'LineWidth', line_width, ...
    'Color', colors_vect(wind.WS_len+1,:))
  leg{1} =  ['Computed ref.'];
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
    leg{i + 1} = ['Sim. ', num2str(i)];
      
  end
  i=1;
  series_length = out_cell{i}.(series).TimeInfo.Length;
  wind_resampled = zeros(series_length, 1);
  wind_resampled = interp1(out_cell{i}.wind.Time,out_cell{i}.wind.Data, ...
    out_cell{i}.(series).Time);
  wind_resampled2 = interp1(out_cell{2}.wind.Time,out_cell{2}.wind.Data, ...
    out_cell{2}.(series).Time);
  [~, zoom_start] = min(abs(wind_resampled - min_zoom_V0));
  [~, zoom_end] = min(abs(wind_resampled2 - max_zoom_V0));
  [zoom_y_min] = out_cell{i}.(series).Data(zoom_start);
  [zoom_y_max] = out_cell{2}.(series).Data(zoom_end);
  rectangle('Position', [min_zoom_V0 zoom_y_min/1e6 ...
  max_zoom_V0-min_zoom_V0  (zoom_y_max-zoom_y_min)/1e6])
  line([min_zoom_V0 axis_pos_x*25], [zoom_y_min/1e6 axis_pos_y*12], 'Color', 'k')
  line([max_zoom_V0 (axis_pos_x+axis_scale_x)*25], [zoom_y_max/1e6 (axis_pos_y+axis_scale_y)*12], 'Color', 'k')
  legend(leg, 'Location', 'northwest', 'FontSize', font_size,...
  'interpreter','latex');
  xlabel(x_label,'interpreter','latex')
  ylabel(y_label,'interpreter','latex')
  title(plot_title, 'Interpreter','latex')
  ylim([0 12])
  grid on
  set(gca, 'FontSize', font_size)
  box on

  axes('position',[axis_pos_x axis_pos_y axis_scale_x axis_scale_y])

  for i=1:wind.WS_len
    plot(wind_resampled(zoom_start:zoom_end), out_cell{i}.(series).Data(zoom_start:zoom_end)/scaling, ...
    'LineWidth', line_width, 'Color', colors_vect(i,:)); %axis tight
    hold on
    % plot(out_cell{1}.T_G.Time(N_start:N_end), ...
    %   out_cell{1}.T_G.Data(N_start:N_end), 'LineWidth', 1.5, ...
    %   'MarkerSize',6, 'Color', colors_vect(2,:)); %axis tight
  end
  xlim([min_zoom_V0 max_zoom_V0])
  ylim([zoom_y_min zoom_y_max]/1e6);
  set(gca, 'FontSize', font_size)
  grid on
  if simulation.print_figure == 1
    fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
    export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);
  %   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
  end
end