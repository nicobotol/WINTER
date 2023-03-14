% This function plots the static results, averaging the response in the
% final part of the simulation

omega_R = zeros(wind.WS_len, 1);
P_G = zeros(wind.WS_len, 1);
pitch = zeros(wind.WS_len, 1);
for i=1:wind.WS_len
  print_time = simulation.stop_time(i) - simulation.plot_time(i);
  [~, pos_start] = min(out_store{i}.omega_R.Time - print_time);

  omega_R(i) = mean(out_store{i}.omega_R.Data(pos_start:end))*30/pi;
  
  P_G(i) = mean(out_store{i}.P_G.Data(pos_start:end));
  
  pitch(i) = rad2deg(mean(out_store{i}.pitch.Data(end - pos_start:end)));
end

% rotor rotational speed
plot_static('fig_omega', wind.mean, omega_R, reference(:,1), ...
  reference(:,3), lookup_static_values(1, pos_4:pos_25), ...
  lookup_static_values(2, pos_4:pos_25)*30/pi, 'Time [s]', ...
  '$\omega_{r}$ [rpm]', 'Rotor rotational speed', 1, date_fig)

% generator power
plot_static('fig_power', wind.mean, P_G, reference(:,1), ...
  reference(:,6)/1e3, lookup_static_values(1, pos_4:pos_25), ...
  lookup_static_values(7, pos_4:pos_25)/1e6, 'Time [s]', ...
  'P [MW]', 'Generator power', 1e6, date_fig) 

% pitch angle
plot_static('fig_pitch', wind.mean, pitch*180/pi, reference(:,1), ...
  reference(:,2), lookup_Pitch(1, pos_4:pos_25), ...
  lookup_Pitch(3, pos_4:pos_25)*180/pi, 'Time [s]', ...
  '$\theta \ [deg.]$', 'Collective blade pitch angle', 180/pi, date_fig)
