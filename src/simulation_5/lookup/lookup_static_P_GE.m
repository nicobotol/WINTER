% This function computes the static values for the electro power loading the rotor one

parameters;
Rs = generator.Rs; % [Ohm] stator resistance
p = generator.p; % [-] number of pole pairs
Lambda = generator.Lambda; % [Vs] flux linkage

V0_v = lookup_static_values(1, :);
omega_R = lookup_static_values(2, :);
P_R = lookup_static_values(6, :);

omega_v = lambda_GE*11.4/rotor.R; % [rad/s] rotor speed

P_G = P_R - B_eq*omega_v.^2; % [W] mechancial power to the generator
iq = 2/3*P_G./(Lambda.*omega_v*p);  % [A] generator current
uq =  omega_v*p*Lambda - Rs*iq;     % [V] generator voltage
P_electro = 1.5*uq.*iq;           % [W] electrical power
P_joule = 1.5*Rs*iq.^2;           % [W] Joule losses

B_eq = 0;
P_G = P_R - B_eq*omega_v.^2; % [W] mechancial power to the generator
iq = 2/3*P_G./(Lambda.*omega_v*p);  % [A] generator current
uq =  omega_v*p*Lambda - Rs*iq;     % [V] generator voltage
P_electro_no_B = 1.5*uq.*iq;           % [W] electrical power
P_joule_no_B = 1.5*Rs*iq.^2;           % [W] Joule losses

% static electro curves
fig_static_electro_power = figure('Color','w');
hold on
plot(V0_v, P_R/1e6, 'LineWidth', line_width, 'DisplayName', '$P_R$', 'Color', colors_vect(1,:));
plot(V0_v, P_electro/1e6, 'LineWidth', line_width, 'DisplayName', '$P_{GE}$', 'Color', colors_vect(3,:));
plot(V0_v, P_joule/1e6, '--','LineWidth', line_width, 'DisplayName', '$P_{joule}$', 'Color', colors_vect(3,:));
plot(V0_v, P_electro_no_B/1e6, 'LineWidth', line_width, 'DisplayName', '$P_{GE}$ B=0$[\frac{kgm^2}{s}]$', 'Color', colors_vect(2,:));
plot(V0_v, P_joule_no_B/1e6, '--', 'LineWidth', line_width, 'DisplayName', '$P_{joule}$ B=0$[\frac{kgm^2}{s}]$', 'Color', colors_vect(2,:));
legend('Location', 'northwest');
xlabel('$V_{0}$ [m/s]')
ylabel('P [MW]')
grid on
box on
if simulation.print_figure == 1
  export_figure(fig_static_electro_power, '\fig_static_electro_power.eps', path_images);
end