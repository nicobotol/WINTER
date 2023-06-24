% wind speed
plot_time_series("fig_wind_TS",out_store, 'wind', 'Time [s]', ...
  'Wind speed [m/s]', 'Wind speed time serie', 1, date_fig, 'northwest');

% pitch dynamic
plot_time_series('fig_pitch_dynamic', out_store, 'pitch','Time [s]',...
   '$\theta \ [deg.]$', 'Pitch angle time serie', pi/180, date_fig, 'northwest'); 

% rotor omega dynamic
plot_time_series('fig_omega_dynamic',out_store,'omega_R','Time [s]',...
  '$\omega_{r}$ [rad/s]', 'Rotor rotational speed', 1, date_fig, 'southwest'); 

% power dynamic
plot_time_series2('fig_power_dynamic',out_store,'P_R','P_G', ...
  rotor.P_rated,'Time [s]','P [MW]','Rotor and generator powers', ...
  1e6,'Aero.','Gen.',date_fig, 'south')

% torque dynamic
plot_time_series2('fig_torque_dynamic',out_store,'T_G_reference','T_G', ...
  'none','Time [s]','T [MNm]','Generator torque', ...
  1e6,'Ref.','',date_fig, 'southeast')

% generator power dynamic
plot_time_series2('fig_generator_power_dynamic',out_store,'P_GE','P_G', ...
  rotor.P_rated,'Time [s]','P [MW]','Generator input and output power', ...
  1e6,'Electro','Mech.',date_fig,'southeast');

for i=1:wind.WS_len % sum the electrical and joule losses
  out_store{i}.P_check = [];
  out_store{i}.P_check.Data = out_store{i}.P_GE.Data + ...
    out_store{i}.P_GJoule.Data + out_store{i}.P_GInductance.Data;
  out_store{i}.P_check.Time = out_store{i}.P_GE.Time;
end

% check on the generator powers
data = string(4);
data = ["P_GE","P_GJoule", "P_G"];
leg = string(4);
leg = ["Electro","Joule loss", "Input"];
plot_time_seriesN('fig_generator_power_check',out_store, data, ...
  rotor.P_rated,'Time [s]','P [MW]','Generator powers', ...
  1e6,leg,date_fig, 'east')

% generator efficiency
for i=1:wind.WS_len % sum the electrical and joule losses
  out_store{i}.generator_eta = [];
  out_store{i}.generator_eta.Data = out_store{i}.P_GE.Data ./...
  out_store{i}.P_G.Data;
  out_store{i}.generator_eta.Time = out_store{i}.P_GE.Time;
end
% plot_time_series("fig_generator_eta",out_store, 'generator_eta', 'Time [s]', ...
%   'Generator efficiency [-]', 'Generator efficiency', 1, date_fig, 'southeast');