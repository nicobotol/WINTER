function plot_step_response(out_cell, date_fig)

  parameters
  N_start = ceil(step_t_start/simulation.time_step_L); % item where start 
  N_end = ceil(step_t_stop/simulation.time_step_L);    % item where stop 
  dt_ref = diff(out_cell{1}.T_G_reference.Time);
  N_start2 = ceil(step_t_start/dt_ref(1)); % item where start 
  [~, N_end2]=min(abs(out_cell{1}.T_G_reference.Time-step_t_stop));

  fig_gen_step = figure('Position', get(0, 'Screensize'),'Color','w');
  plot(out_cell{1}.T_G_reference.Time, out_cell{1}.T_G_reference.Data,'--','LineWidth', line_width);
  hold on
  plot(out_cell{1}.T_G.Time, out_cell{1}.T_G.Data, 'LineWidth', line_width*0.75);
  legend('Reference', 'Response', 'location', 'northwest', 'FontSize', font_size, 'interpreter','latex');
  grid on;
  title('Generator q-axis current step response')
  xlabel('Time [s]', 'interpreter','latex')
  ylabel('$I_q$ [A]', 'interpreter','latex')
  set(gca, 'FontSize', font_size)
  rectangle('Position', [step_t_start step_y_min step_t_stop-step_t_start  step_y_max-step_y_min])
  ylim([0, 1.2])

  axes('position',[2.8/5 0.4/1.2 .30 .30])
  plot(out_cell{1}.T_G.Time(N_start:N_end), out_cell{1}.T_G.Data(N_start:N_end), 'LineWidth', 1.5, 'MarkerSize',6, 'Color', colors_vect(2,:)); %axis tight
  hold on
  plot(out_cell{1}.T_G_reference.Time(N_start2:N_end2), out_cell{1}.T_G_reference.Data(N_start2:N_end2), '--', 'LineWidth', 1.5,'MarkerSize',6,'Color',colors_vect(1,:)); %axis tight
  xlim([step_t_start step_t_stop])
  ylim([step_y_min step_y_max]);
  grid on
  grid minor
  title('Overshoot of the response')
  set(gca, 'FontSize', font_size)

if simulation.print_figure == 1
  plot_name = 'generator_step_input';
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig_gen_step, strcat(date_fig, plot_name, '.eps'), path_images);
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
end
end