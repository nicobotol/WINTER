%% Plot of the results
% This file is used for print the results of the simulation
close all

date = string(datetime('now','TimeZone','local','Format', ...
        'y_MM_d_HH_mm_ss'));  % save the date to identify the figures
pos_25=V0_cut_out/velocity_spacing; %pos. of 25 [m/s] in static_values(1,:)
pos_4 = V0_cut_in/velocity_spacing; %pos. of 4 [m/s] in static_values(1,:)

switch simulation.type
  case 1 % constant wind speed
 
    % plot all the static simulations, averaging the response in the last
    % part
    plot_all_static
    
    % plot the time dependent graphs
    plot_all_dynamic

  case 2 % ramp
    plot_all_dynamic

  case 3 % generated wind series
    plot_all_dynamic

  case 4 % generator step response
    plot_step_response(out_cell)

  case 5 % generated WS and parametrization plot
    plot_all_dynamic
    
    plot_all_parametrizations

  case 6 % ramp and parametrization plot
    plot_all_dynamic

    plot_all_parametrizations
    
end
