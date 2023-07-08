%% Plot of the results
% This file is used for print the results of the simulation
% close all

addpath('\..')
addpath('lookup\')


[~, pos_4] = min(abs(lookup_static_values(1, :) - V0_cut_in)/ ...
  velocity_spacing); % position of the element coresponding to 4 [m/s]
[~, pos_25] = min(abs(lookup_static_values(1, :) - V0_cut_out)/ ...
  velocity_spacing); % position of the element coresponding to 25 [m/s]

switch simulation.type
  case 1 % constant wind speed
    % plot all the static results
    plot_all_static

    % plot the time dependent graphs
    plot_all_dynamic

    plot_all_parametrizations_wind

  case 2 % ramp
    plot_all_dynamic

  case 3 % generated wind series
    plot_all_dynamic

  case 4 % generator step response
    plot_step_response(out_store, date_fig)

  case 5 % generated WS and parametrization plot   
    plot_all_dynamic
    plot_all_parametrizations_wind
    plot_all_parametrizations

  case {6, 7, 8} % ramp and parametrization plot
    plot_all_dynamic
    plot_all_parametrizations_wind
%     plot_all_parametrizations
    [~, ~, ~] = power_check(out_store, I_eq, B_eq, wind.WS_len,...
      line_width, date_fig, "power_check", simulation, font_size, path_images);

  case 9 % different pitch mechanis dynamic
    plot_all_dynamic

  case 10 % control based on K_opt and K_opt_GE
    plot_all_parametrizations_wind_P_GE
    [E_R, E_G, E_GE] = power_check(out_store, I_eq, B_eq, wind.WS_len,...
      line_width, date_fig, "power_check", simulation, font_size, path_images);
end