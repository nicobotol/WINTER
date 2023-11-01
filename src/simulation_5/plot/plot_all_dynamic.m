% wind speed
plot_time_series("fig_wind_TS",out_store, 'wind', 'Time [s]',  'Wind speed [m/s]', 'Wind speed time serie', 1, date_fig, 'northwest');

% pitch dynamic
plot_time_series('fig_pitch_dynamic', out_store, 'pitch','Time [s]','$\theta \ [deg.]$', 'Pitch angle time serie', pi/180, date_fig, 'southeast'); 

% rotor omega dynamic
data = string(1);
data = ["omega_R"];
leg = string(4);
leg = [ "Sim. 1", "Sim. 2", "Sim. 3", "Sim. 4"];
plot_time_seriesN('fig_omega_dynamic', out_store, data, omega_rated, 'Time [s]','$\omega_R$ [rad/s]', 'Rotor rotational speed', 1, line_width, leg, date_fig, 'southeast')

% power dynamic
plot_time_series2('fig_power_dynamic', out_store, 'P_R', 'P_G',  generator.P_rated, 'Time [s]', 'P [MW]', 'Rotor and generator powers', 1e6, 'Aero.', 'Gen.', date_fig, 'southeast')

% torque dynamic
plot_time_series2('fig_torque_dynamic', out_store,'T_G_reference', 'T_G', 'none', 'Time [s]', 'T [MNm]', 'Generator torque', 1e6, 'Ref.', '', date_fig, 'southeast')

% generator power dynamic
plot_time_series2('fig_generator_power_dynamic', out_store,'P_GE', 'P_G', generator.P_rated, 'Time [s]', 'P [MW]', 'Generator input and output power', 1e6, 'Electro', 'Mech.', date_fig, 'southeast');

% rotor and generator output power dynamic
plot_time_series3('fig_rotor_generator_power_dynamic', out_store,'P_R', 'P_GE', generator.P_rated, 'Time [s]', 'P [MW]', 'Rotor input $P_R$ and Generator output $P_{GE}$ power', 1e6, '$P_{R}$', '$P_{GE}$', date_fig, 'bestoutside', 2);

if simulation.type ~= 10
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
  plot_time_seriesN('fig_generator_power_check', out_store, data, generator.P_rated, 'Time [s]', 'P [MW]', 'Generator powers', 1e6, line_width*0.6, leg, date_fig, 'east')
  
  % generator efficiency
  for i=1:wind.WS_len % sum the electrical and joule losses
    out_store{i}.generator_eta = [];
    out_store{i}.generator_eta.Data = out_store{i}.P_GE.Data./out_store{i}.P_G.Data;
    out_store{i}.generator_eta.Time = out_store{i}.P_GE.Time;
  end
end
% plot_time_series("fig_generator_eta",out_store, 'generator_eta', 'Time [s]', ...
%   'Generator efficiency [-]', 'Generator efficiency', 1, date_fig, 'southeast');