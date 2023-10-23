%% Plot of the results
% This file is used for print the results of the simulation
% close all

addpath('\..')
addpath('lookup\')


[~, pos_4] = min(abs(lookup_static_values(1, :) - V0_cut_in)/  velocity_spacing); % position of the element coresponding to 4 [m/s]
[~, pos_25] = min(abs(lookup_static_values(1, :) - V0_cut_out)/  velocity_spacing); % position of the element coresponding to 25 [m/s]
[~, pos_ramp_start] = min(abs(lookup_static_values(1, :) - wind.ramp_WS_start(1))/  velocity_spacing); % position of the element coresponding to 25 [m/s]
[~, pos_ramp_stop] = min(abs(lookup_static_values(1, :) - wind.ramp_WS_stop(1))/  velocity_spacing); % position of the element coresponding to 25 [m/s]

switch simulation.type
  case 1 % constant wind speed
    % plot all the static results
    % plot_all_static
    
    % plot the time dependent graphs
    plot_all_dynamic
    plot_IMM
    % plot_all_parametrizations_wind

    % [~, ~, ~] = power_check(out_store, I_eq, B_eq, wind.WS_len, line_width, date_fig, "power_check", simulation, font_size, path_images);

  case 2 % ramp
    plot_all_dynamic

  case 3 % generated wind series
    plot_all_dynamic

  case 4 % generator step response
    plot_step_response(out_store, date_fig)

  case 5 % generated WS and parametrization plot   
    plot_time_series("fig_wind_TS",out_store, 'wind', 'Time [s]', 'Wind speed [m/s]', 'Wind speed time serie', 1, date_fig, 'northwest');
    % plot_all_dynamic
    plot_all_parametrizations_wind
    % plot_all_parametrizations

  case {6, 7, 8} % ramp and parametrization plot
    % plot_all_dynamic
    % wind speed
    % plot_time_series("fig_wind_TS",out_store, 'wind', 'Time [s]', 'Wind speed [m/s]', 'Wind speed time serie', 1, date_fig, 'northwest');
    plot_all_parametrizations_wind
    plot_IMM
    % plot_all_parametrizations
    % [~, ~, ~] = power_check(out_store, I_eq, B_eq, wind.WS_len, line_width, date_fig, "power_check", simulation, font_size, path_images);
    % rotor power parametrization

  case 9 % different pitch mechanis dynamic
    plot_all_dynamic

  case 10 % control based on K_opt and K_opt_GE
    plot_all_parametrizations_wind_P_GE
    % plot_all_dynamic;
    plot_parametrization_wind_P_GE_vs_P('comparison_control_laws',out_store,'P_GE', lookup_P_GE(1,:), lookup_P_GE(2,:)/1e6,'Wind speed [m/s]','P [MW]','Comparison of the control laws',1e6,simulation,date_fig)
    % [E_R, E_G, E_GE] = power_check(out_store, I_eq, B_eq, wind.WS_len,line_width, date_fig, "power_check", simulation, font_size, path_images);

  case 11 % control based on different K_opt and the 
    plot_parametrization_wind_zoom_P_GE('fig_electrical_power_param_zoom',out_store,'P_GE', lookup_P_GE(1,:), lookup_P_GE(3,:)/1e6,'Wind speed [m/s]','$P_{GE}$ [MW]','Generator electrical power', 1e6,simulation,date_fig)

end