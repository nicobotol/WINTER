function plot_step_response(out_cell, date_fig)

  parameters
  N_start = ceil(step_t_start/simulation.time_step_L); % item where start 
  N_end = ceil(step_t_stop/simulation.time_step_L);    % item where stop 

  fig_gen_step = figure('Position', get(0, 'Screensize'),'Color','w');
  plot(out_cell{1}.T_G_reference.Time, out_cell{1}.T_G_reference.Data, ...
    '--','LineWidth', line_width);
  hold on
  plot(out_cell{1}.T_G.Time, out_cell{1}.T_G.Data, 'LineWidth', ...
    line_width*0.75);
  legend('Reference', 'Response', 'location', 'best', 'FontSize', ...
    font_size, 'interpreter','latex');
  grid on;
  title('Generator torque step response')
  xlabel('Time [s]', 'interpreter','latex')
  ylabel('T [Nm]', 'interpreter','latex')
  set(gca, 'FontSize', font_size)
  rectangle('Position', [step_t_start step_y_min ...
    step_t_stop-step_t_start  step_y_max-step_y_min])

  axes('position',[2.8/5 0.4/1.2 .30 .30])
  plot(out_cell{1}.T_G_reference.Time(N_start:N_end), ...
    out_cell{1}.T_G_reference.Data(N_start:N_end), '--', ...
    'LineWidth', 1.5,'MarkerSize',6,'Color',colors_vect(1,:)); %axis tight
  hold on
  plot(out_cell{1}.T_G.Time(N_start:N_end), ...
    out_cell{1}.T_G.Data(N_start:N_end), 'LineWidth', 1.5, ...
    'MarkerSize',6, 'Color', colors_vect(2,:)); %axis tight
  xlim([step_t_start step_t_stop])
  ylim([step_y_min step_y_max]);
  grid on
  title('Overshoot of the response')
  set(gca, 'FontSize', font_size)

  if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, 'fig_gen_step','.svg');
  export_fig('fig_gen_step', fig_name);
  end
end