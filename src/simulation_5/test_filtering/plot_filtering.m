data = string(2);
data = ["K_opt_no_filter", "K_opt"];
leg = string(2);
leg = ["Filtered", "Not filtered"];
plot_K('fig_K', out_store, data, generator.K_opt_GE, 'Time [s]', 'K [$Nms^2$]', 'Effect of the filter on the K gain', 1, line_width, leg, date_fig, 'southeast')

function plot_K(plot_name, out_cell, data, y_line,  x_label, y_label, plot_title, scaling, font_scale, leg, date_fig, leg_loc)
line_style = ["-", "-", ":", "-."];

parameters
t_start = 100; % time from when plot data
t_stop = 150; % time when stop plot data

N = size(data, 2); % data to be plotted
% leg = cell(N*wind.WS_len, 1);
fig = figure('Color','w');
hold on
for i = 1:wind.WS_len
    
  for j = 1:N
    [~, s_start] = min(abs(out_cell{i}.(data(j)).Time - t_start)); % sample from where to start
    [~, s_stop] = min(abs(out_cell{i}.(data(j)).Time - t_stop)); % sample from where to start
    plot(out_cell{i}.(data(j)).Time(s_start:s_stop), out_cell{i}.(data(j)).Data(s_start:s_stop)/scaling, 'LineWidth', font_scale, 'Color', color(j), 'LineStyle',line_style(j));
%     leg{N*i - N*j} = [leg(j),' sim. ', num2str(i)];
  end
end
if strcmp(y_line, 'none') == 0
  yline(y_line/scaling, 'LineStyle', '-.', 'LineWidth', line_width*1.7, 'Color', color(5));
  leg{N*wind.WS_len + 1} = ['Rated'];
end
legend(leg, 'Location', leg_loc, 'FontSize', font_size, 'interpreter', 'latex', 'NumColumns', wind.WS_len);
xlabel(x_label, 'interpreter', 'latex')
ylabel(y_label, 'interpreter', 'latex')
title(plot_title, 'Interpreter', 'latex')
grid on
box on
set(gca, 'FontSize', font_size)
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
end
end