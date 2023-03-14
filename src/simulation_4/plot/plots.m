%% Plot of the results
% This file is used for print the results of the simulation
close all

addpath('\..')
addpath('lookup\')

date_fig = string(datetime('now','TimeZone','local','Format', ...
        'y_MM_d_HH_mm_ss'));  % save the date to identify the figures
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

  case 2 % ramp
    plot_all_dynamic

  case 3 % generated wind series
    plot_dynamic

  case 4 % generator step response
    plot_step_response(out_store, date_fig)

  case 5 % generated WS and parametrization plot
    plot_all_dynamic
    
    plot_all_parametrizations

  case 6 % ramp and parametrization plot
    plot_dynamic

    plot_all_parametrizations
    
end