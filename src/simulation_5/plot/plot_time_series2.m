function plot_time_series2(plot_name,out_cell,series1,series2,y_line, ...
  x_label,y_label,plot_title,scaling,leg1,leg2,date_fig, leg_loc)

parameters

leg = cell(2*wind.WS_len, 1);
fig = figure('Color','w');
hold on
y_lim_min = 1e10;
y_lim_max = -1e10;

for i = 1:wind.WS_len
    % Time from where start to print
    t_start = out_cell{i}.(series1).Time(end) - simulation.plot_time(i); % [s]
    [~, s_start] = min(abs(out_cell{i}.(series1).Time - t_start)); % sample from where to start

  plot(out_cell{i}.(series1).Time(s_start:end), out_cell{i}.(series1).Data(s_start:end)/scaling, ...
    'LineWidth',line_width*0.6,'Color',color(i),'LineStyle','--');
  plot(out_cell{i}.(series2).Time(s_start:end), out_cell{i}.(series2).Data(s_start:end)/scaling, ...
    'LineWidth', line_width*1.2, 'Color', color(i));
  leg{2*i - 1} = [leg1,' sim. ', num2str(i)];
  leg{2*i} = [leg2,' sim. ', num2str(i)];

  y_min = min(out_cell{i}.(series1).Data(s_start:end)/scaling);
  y_max = max(out_cell{i}.(series1).Data(s_start:end)/scaling);
  if y_min < y_lim_min
    y_lim_min = y_min;
  end
  if y_max > y_lim_max
    y_lim_max = y_max;
  end

end
if strcmp(y_line, 'none') == 0
  yline(y_line/scaling, 'LineStyle', '-.', 'LineWidth', line_width,...
    'Color', color(i+1));
  leg{2*wind.WS_len + 1} = ['Rated'];
end
legend(leg, 'Location', leg_loc, 'FontSize', font_size,...
  'interpreter','latex','NumColumns',wind.WS_len);
xlabel(x_label,'interpreter','latex')
ylabel(y_label,'interpreter','latex')
title(plot_title, 'Interpreter','latex')
grid on
box on
set(gca, 'FontSize', font_size)
ylim([0.95*y_lim_min 1.05*y_lim_max])

if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
end