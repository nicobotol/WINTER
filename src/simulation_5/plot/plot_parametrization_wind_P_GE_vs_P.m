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
line_type = {':', '-', '--', '-.'};
k=0;
for i=1:2:wind.WS_len
  % k=k+1;
  if wind.WS_len ~= 2
    k=k+1;
    subplot(2,3,k)
    line_type = {'x', 'o', 'x', 'o'};
  end
  hold on
  % Time from where start to print
  t_start = out_cell{i}.P_GE.Time(end) - simulation.plot_time(i); % [s]
  [~, s_start] = min(abs(out_cell{i}.P_GE.Time - t_start)); % sample from where to start
  t_start_ii = out_cell{i+1}.P_GE.Time(end) - simulation.plot_time(i+1); % [s]
  [~, s_start_ii] = min(abs(out_cell{i+1}.P_GE.Time - t_start)); % sample from where to start
  
  wind_resampled = out_cell{i}.wind.Data;
  wind_resampled_ii = out_cell{i+1}.wind.Data;
  
  plot(wind_resampled_ii(s_start_ii:end), out_cell{i+1}.P_R.Data(s_start_ii:end)/scaling, line_type{1}, 'LineWidth', line_width, 'Color', color(1), 'MarkerSize', marker_size); % rotor mechanical power
  plot(wind_resampled_ii(s_start_ii:end), out_cell{i+1}.P_GE.Data(s_start_ii:end)/scaling, line_type{2}, 'LineWidth',  line_width, 'Color', color(1), 'MarkerSize', marker_size); % rotor electrical power
  plot(wind_resampled(s_start:end), out_cell{i}.P_R.Data(s_start:end)/scaling, line_type{3},'LineWidth',  line_width, 'Color', color(2), 'MarkerSize', marker_size); % generator mechanical power
  plot(wind_resampled(s_start:end), out_cell{i}.P_GE.Data(s_start:end)/scaling, line_type{4},'LineWidth', line_width, 'Color', color(2), 'MarkerSize', marker_size); % generator electrical power

  P_R_R(k) = mean(out_cell{i+1}.P_R.Data(s_start_ii:end)); % rotor mechanical power
  P_GE_R(k) = mean(out_cell{i+1}.P_GE.Data(s_start_ii:end)); % rotor electrical power
  P_R_G(k) = mean(out_cell{i}.P_R.Data(s_start:end)); % generator mechanical power
  P_GE_G(k) = mean(out_cell{i}.P_GE.Data(s_start:end)); % generator electrical power

  E_R(k) =  (P_R_G(k) - P_R_R(k))/P_R_G(k)*100;
  E_GE(k) = (P_GE_G(k) - P_GE_R(k))/P_GE_G(k)*100;

  grid on
  set(gca, 'FontSize', 0.9*font_size)
  box on
  if wind.WS_len ~= 2
    if (k == 1 || k == 4) 
      ylabel(y_label,'interpreter','latex')
    end
  else
    ylabel(y_label,'interpreter','latex')
  end
    ylim([0.99*min(out_cell{i}.P_GE.Data(s_start:end)/scaling), 1.005*max(out_cell{i}.P_R.Data(s_start:end)/scaling)])

end
xlabel(x_label,'interpreter','latex')
legends_name = {'Rotor $P_{R}$','Rotor $P_{GE}$','Gen. $P_{R}$','Gen. $P_{GE}$'};
if wind.WS_len ~= 2
  legend(legends_name, 'position', [0.75,0.238193265647898,0.127430735291504,0.058136925199618]);
else
  legend(legends_name, 'location', 'southeast');
end
% legend('Rotor $P_{in}$','Gen. $P_{in}$','Gen. $P_{out}$','Rotor $P_{out}$', 'Location', 'best', 'FontSize', 0.9*font_size,'interpreter','latex');
sgtitle(plot_title, 'Interpreter','latex', 'FontSize', 0.9*font_size);
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
end

name = string();
text = string();
if B_eq < 1
  name = strcat('../../report/macro/', date_fig, 'powers_error_no_B.tex');
  text = 'Without $B_{eq}$';
else 
  name = strcat('../../report/macro/', date_fig, 'powers_error.tex');
  text = 'With $B_{eq}$';
end

fileID = fopen(name,'w');
fprintf(fileID, '\\multirow{5}{*}{%s} \n', text);
for i=1:wind.WS_len/2
fprintf(fileID,' & %.2f & %.3f & %.3f & %.2f & %.3f & %.3f & %.2f\\\\ \n', [wind.mean(2*i), P_R_G(i)/1e6, P_R_R(i)/1e6, E_R(i), P_GE_G(i)/1e6, P_GE_R(i)/1e6, E_GE(i)]');
end
fclose(fileID);
