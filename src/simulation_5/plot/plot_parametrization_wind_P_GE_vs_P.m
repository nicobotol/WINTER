function plot_parametrization_wind_P_GE_vs_P(plot_name,out_cell,series,x_my_ref,y_my_ref,x_label,y_label,plot_title,scaling,simulation,date_fig)

parameters

leg = cell(1, wind.WS_len + 1);
fig = figure('Color','w');

hold on
% plot([NaN, NaN], [NaN, NaN], 'ko')
% plot([NaN, NaN], [NaN, NaN], 'kx')
% plot([NaN, NaN], [NaN, NaN], 'kdiamond')
% plot([NaN, NaN], [NaN, NaN], 'ks')
% plot(x_my_ref, y_my_ref, '--', 'LineWidth', line_width, 'Color', color(wind.WS_len+1))
% leg{1} =  ['Computed ref.'];
k=0;
for i=1:2:wind.WS_len
  k=k+1;
  subplot(2,3,k)
  hold on
  % Time from where start to print
  t_start = out_cell{i}.P_GE.Time(end) - simulation.plot_time(i); % [s]
  [~, s_start] = min(abs(out_cell{i}.P_GE.Time - t_start)); % sample from where to start
  t_start_ii = out_cell{i+1}.P_GE.Time(end) - simulation.plot_time(i+1); % [s]
  [~, s_start_ii] = min(abs(out_cell{i+1}.P_GE.Time - t_start)); % sample from where to start
  
  % Resample the wind speed
  series_length = out_cell{i}.P_GE.TimeInfo.Length;
  wind_resampled = zeros(series_length, 1);
  wind_resampled = interp1(out_cell{i}.wind.Time,out_cell{i}.wind.Data, ...
    out_cell{i}.P_GE.Time);
  series_length_ii = out_cell{i+1}.P_GE.TimeInfo.Length;
  wind_resampled_ii = zeros(series_length_ii, 1);
  wind_resampled_ii = interp1(out_cell{i+1}.wind.Time,out_cell{i+1}.wind.Data, out_cell{i+1}.P_GE.Time);
  
  plot(wind_resampled_ii(s_start_ii:end), out_cell{i+1}.P_G.Data(s_start_ii:end)/scaling,'s', 'LineWidth', line_width, 'Color', color(4)); % rotor mechanical power
  plot(wind_resampled(s_start:end), out_cell{i}.P_G.Data(s_start:end)/scaling, 'diamond','LineWidth',  line_width, 'Color', color(3)); % generator mechanical power
  plot(wind_resampled(s_start:end), out_cell{i}.P_GE.Data(s_start:end)/scaling, 'o','LineWidth', line_width, 'Color', color(1)); % generator electrical power
  plot(wind_resampled_ii(s_start_ii:end), out_cell{i+1}.P_GE.Data(s_start_ii:end)/scaling,'x', 'LineWidth',  line_width, 'Color', color(2)); % rotor electrical power
  grid on
  set(gca, 'FontSize', font_size)
  box on
  if k == 1 || k == 4
    ylabel(y_label,'interpreter','latex')
  end
  ylim([0.995*min(out_cell{i}.P_GE.Data(s_start:end)/scaling), 1.005*max(out_cell{i}.P_G.Data(s_start:end)/scaling)])
end
xlabel(x_label,'interpreter','latex')
legends_name = {'Rotor $P_{in}$','Gen. $P_{in}$','Gen. $P_{out}$','Rotor $P_{out}$'};
legend(legends_name, 'position', [0.75,0.238193265647898,0.127430735291504,0.058136925199618]);
% legend('Rotor $P_{in}$','Gen. $P_{in}$','Gen. $P_{out}$','Rotor $P_{out}$', 'Location', 'best', 'FontSize', font_size,'interpreter','latex');
sgtitle(plot_title, 'Interpreter','latex', 'FontSize', font_size);
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
end