% This function plots different parametrizations with the wind speed n the
% x axes

addpath('\..')
addpath('lookup\')

% generator input power parametrization
plot_parametrization_wind_P_GE('fig_power_param',out_store,'P_G', ...
 lookup_static_values(1,pos_4:pos_25), ...
 lookup_static_values(6,pos_4:pos_25)/1e6,'Wind speed [m/s]','P [MW]', ...
 'Generator input power',1e6,simulation,date_fig)

% generator output power parametrization
plot_parametrization_wind_P_GE('fig_electrical_power_param',out_store,'P_GE', ...
lookup_P_GE(1,:), ...
lookup_P_GE(2,:)/1e6,'Wind speed [m/s]','P [MW]', ...
'Generator electrical power',1e6,simulation,date_fig)

% generator torque parametrization
plot_parametrization_wind_P_GE('fig_torque_param',out_store,'T_G', ...
  lookup_static_values(1,pos_4:pos_25), ...
  lookup_static_values(5,pos_4:pos_25)/1e6,'Wind speed [m/s]','T [MNm]', ...
  'Generator torque',1e6,simulation,date_fig)

% pitch parametrization
plot_parametrization_wind_P_GE('fig_pitch_param',out_store,'pitch', ...
  lookup_Pitch(1,pos_4:pos_25),lookup_Pitch(3,pos_4:pos_25)*180/pi, ...
  'Wind speed [m/s]','$\theta$ [deg.]','Pitch angle',pi/180,simulation,date_fig)

% pitch parametrization with map made based on P_GE
x_my_ref = cell(2);
y_my_ref = cell(2);
x_my_ref{1} = lookup_Pitch(1,pos_4:pos_25);
y_my_ref{1} = lookup_Pitch(3,pos_4:pos_25)*180/pi; 
x_my_ref{2} = lookup_pitch_P_GE(1,:);
y_my_ref{2} = lookup_pitch_P_GE(2,:)*180/pi;
leg = string(2);
leg = ["Ref. rotor", "Ref. generator"];
plot_parametrization_wind_P_GEN('fig_pitch_param',out_store,'pitch', x_my_ref,y_my_ref, 'Wind speed [m/s]','$\theta$ [deg.]','Pitch angle',pi/180,simulation,leg,date_fig)

% rotational speed
plot_parametrization_wind_P_GE('fig_omega_param',out_store,'omega_R', lookup_static_values(1,pos_4:pos_25), lookup_static_values(2,pos_4:pos_25)*30/pi, 'Wind speed [m/s]','$\omega$ [rpm]','Rotor rotational speed', pi/30,simulation,date_fig)

  if simulation.type == 10 % comparison K_opt and K_opt_GE
%     plot_parametrization_wind_P_GE_zoom_P_GE('fig_electrical_power_param_zoom',out_store,'P_GE', lookup_P_GE(1,:), lookup_P_GE(2,:)/1e6,'Wind speed [m/s]','$P_{GE}$ [MW]','Generator electrical power', 1e6,simulation,date_fig)

%     plot_parametrization_wind_P_GE_zoom_P_GE('fig_power_param_zoom',out_store,'P_G', lookup_static_values(1,pos_4:pos_25), ...
%     lookup_static_values(6,pos_4:pos_25)/1e6,'Wind speed [m/s]','P [MW]', ...
%     'Generator input power', 1e6,simulation,date_fig)
% 
%     plot_parametrization_wind_P_GE_zoom_P_GE('fig_power_param_zoom',out_store,'P_R', lookup_static_values(1,pos_4:pos_25), ...
%     lookup_static_values(6,pos_4:pos_25)/1e6,'Wind speed [m/s]','P [MW]', ...
%     'Generator input power', 1e6,simulation,date_fig)
  end

  cP1 = out_store{1}.cP.Data;
  cP2=out_store{2}.cP.Data;
  t1 =  out_store{1}.cP.Time;
  t2 = out_store{2}.cP.Time;
  figure();hold on;
  plot(t1, cP1,'DisplayName','Sim1');
  plot(t2, cP2,'DisplayName','Sim2');
  legend('Location','best')
