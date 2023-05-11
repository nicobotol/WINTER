function plot_time_dataN(plot_name,out_cell,data,y_line, ...
  x_label,y_label,plot_title,scaling,leg,date_fig)

line_style = ["-", "--", ":", "-."];

parameters
N = size(data, 2); % data to be plotted
% leg = cell(N*wind.WS_len, 1);
fig = figure('Color','w');
hold on
for i = 1:wind.WS_len
    
  for j = 1:N
    % Time from where start to print
    t_start = out_cell{i}.(data(j)).Time(end) - simulation.plot_time(i); % [s]
    [~, s_start] = min(abs(out_cell{i}.(data(j)).Time - t_start)); % sample from where to start
    plot(out_cell{i}.(data(j)).Time(s_start:end), out_cell{i}.(data(j)).Data(s_start:end)/scaling, ...
      'LineWidth',line_width*0.6,'Color',colors_vect(i, :), 'LineStyle',line_style(j));
%     leg{N*i - N*j} = [leg(j),' sim. ', num2str(i)];
  end
end
if strcmp(y_line, 'none') == 0
  yline(y_line/scaling, 'LineStyle', '-.', 'LineWidth', line_width,...
    'Color', colors_vect(i + 1, :));
  leg{N*wind.WS_len + 1} = ['Rated'];
end
legend(leg, 'Location', 'south', 'FontSize', font_size,...
  'interpreter','latex','NumColumns',wind.WS_len);
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