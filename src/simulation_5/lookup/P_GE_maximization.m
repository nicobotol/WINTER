% clear 
% close all
% clc

parameters;
Rs = generator.Rs; % [Ohm] stator resistance
A = rotor.A; % [m^2] rotor swept area
B = B_eq; % [kgm^2] transmission equiuvalent damping
p = generator.p; % [-] number of pole pairs
Lambda = generator.Lambda; % [Wb] generator flux linkage
eta = generator.eta; % [-] generator efficiency
R = rotor.R; % [m]

%% Below rated wind speed
omega_b = 0.3:0.001:1.2; % [rad/s] rotational speeds to test
pitch = [0:0.001:1]*pi/180; % [rad] pitch angle to be tested
V0_b = 4:0.01:V0_rated+0.001; % [m/s] velocity vector below rated
cP_max = zeros(length(V0_b), 1);
theta_opt = zeros(length(V0_b), 1);
omega_GE = zeros(length(V0_b), 1);
P_GE =  zeros(length(pitch), length(omega_b)); % [W] matrix of powers
P_GE_b = zeros(length(V0_b), 1);  % [W] power map below rated
for i = 1:length(V0_b)
    V0 = V0_b(i);
    lambda_b = omega_b'*R/V0;
    cp = interp2(lambda_vector, pitch_vector, lookup_cP, lambda_b, pitch);
    P_GE = -Rs*(A*V0^3.*cp*rho + 2*B*omega_b.^2).^2./(9*omega_b.^2*p^2*Lambda^2) + A.*cp*rho*V0^3*eta/2; % [W] generator power output

    max_tmp = max(P_GE, [], 2);       
    [~, theta_pos] = max(max_tmp);              % omega corresponding to the maximum                  
    [P_GE_b(i), omega_pos] = max(P_GE(theta_pos, :)); % max cP
    omega_GE(i) = omega_b(omega_pos);            % TSR for cP_max
    theta_opt(i) = pitch(theta_pos);            % pitch for cP_max
end

% Recompute some values
omega_rotor =  lambda_opt.*V0_b/R; % [rad/s] omega based on the maximization of the rotor power
lambda_GE =  omega_GE'*R./V0_b;  % optimal TSR based on generator maximization 
lambda_GE_mean = mean(lambda_GE);  % mean TSR based on the generator
cp_GE = interp2(lambda_vector, pitch_vector, lookup_cP, lambda_GE, theta_opt'); % cp for the maximization of the generator
cp_GE_mean = mean(cp_GE);          % mean cp for the generator maximization

%% Static generator power curve
V0_v = 4:1:25;                    % windspeed [m/s]
lambda = zeros(length(V0_v), 1);  % TSR
P_f = zeros(length(V0_v), 1);     % [W] power
for v=1:length(V0_v)
  V0 = V0_v(v);
  if V0 < V0_rated
    lambda(v) = lambda_opt;
  else
    lambda(v) = omega_rated*rotor.R/V0;
  end
end
omega = lambda.*V0_v'/rotor.R;    % [rad/s] rotional speed
stall(:) = interp1(lookup_Pitch(1, :), lookup_Pitch(2, :), V0_v); % [rad]
feather(:) = interp1(lookup_Pitch(1, :), lookup_Pitch(3, :), V0_v); % [rad]
cP_f(:) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, feather'); % [-] feather power coefficient
P_f(:) = 0.5*rotor.A*rho.*cP_f.*V0_v.^3; % [W] power to the rotor side

iq = 2/3*P_f./(Lambda.*omega*p);  % [A] generator current
uq =  omega*p*Lambda - Rs*iq;     % [V] generator voltage
P_electro = 1.5*uq.*iq;           % [W] electrical power
P_joule = 1.5*Rs*iq.^2;           % [W] Joule losses
P_electro_eta = 1.5*uq.*iq*eta;   % [W] electrical power
P_joule_eta = 1.5*Rs*(iq*eta).^2; % [W] Joule losses


%% Maximization of the genertaor power curve above the rated windspeed
V0_a = V0_rated:0.5:25;                   % [m/s] wind speed above rated one 
theta_p = [-1:0.1:25]*pi/180;       % [rad] range of angle where to search for the optimum 
lambda_a = zeros(length(V0_a), 1);  % [-] TSR
cp_a = zeros(length(theta_p), length(V0_a));
Pg_out_a = zeros(length(theta_p), length(V0_a)); % [W] power above the rated
P_GE_a = zeros(length(V0_a), 1);

for i=1:length(V0_a)
  lambda_a(i) = omega_rated*rotor.R/V0_a(i); % TSR
  cp_a(:, i) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda_a(i), theta_p');
  Pg_out_a(:, i) = -(4*Rs*B^2*omega_rated.^2)/(9*p^2*Lambda^2) - (4*Rs*A.*cp_a(:, i)'*rho.*V0_a(i).^3*B)/(9*p^2*Lambda^2) + A*cp_a(:, i)'*rho*V0_a(i).^3*eta/2 - Rs*A^2*cp_a(:, i)'.^2*rho^2*V0_a(i).^6./(9*p^2*Lambda^2*omega_rated.^2); % [W] output power at the generator
end

[~, max_p_a] = min(abs(Pg_out_a - rotor.P_rated), [], 1); % extract the maximum of the power
for i=1:length(V0_a)
  P_GE_a(i) = Pg_out_a(max_p_a(i), i); % [W] maximum power 
end
feather_a = theta_p(max_p_a); % [rad] angle for feather pitching 

%% Plots
% Below rated
figure('Color', 'w')
hold on; box on; grid on;
plot(V0_b, omega_GE, 'LineWidth',line_width, 'DisplayName','Generator ')
plot(V0_b, omega_rotor, 'LineWidth',line_width, 'DisplayName','Rotor')
title('Rotor rotational speed')
xlabel('$V_0$ [m/s]')
ylabel('$\omega$ [rad/s]')
legend('location', 'northwest')

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, theta_opt*180/pi, 'LineWidth',line_width)
title('Pitch angle')
xlabel('$V_0$ [m/s]')
ylabel('$\theta$ [deg]')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_pitch_GE', '.eps'), path_images); 
end

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, lambda_GE, 'LineWidth',line_width, 'DisplayName','Generator','Color',colors_vect(1,:))
yline(lambda_GE_mean, '--', 'LineWidth',line_width, 'DisplayName','Mean','Color',colors_vect(2,:))
yline(lambda_opt, '--','LineWidth',line_width, 'DisplayName','Rotor','Color',colors_vect(3,:))
title('Tip speed ratio')
xlabel('$V_0$ [m/s]')
ylabel('$\lambda$ [-]')
legend('location', 'east')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_lambda_GE', '.eps'), path_images); 
end

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, cp_GE, 'LineWidth',line_width, 'DisplayName','Generator')
yline(cp_GE_mean, '--', 'LineWidth',line_width, 'DisplayName','Mean','Color',colors_vect(2,:))
yline(cp_max, '--', 'LineWidth',line_width, 'DisplayName','Rotor','Color',colors_vect(3,:))
title('Power coefficient')
xlabel('$V_0$ [m/s]')
ylabel('$c_P$ [-]')
legend('location', 'northeast')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_cP_GE', '.eps'), path_images); 
end

% cp_max_gen = -((-9*Lambda^2*eta*p^2 + 8*B*Rs).*lambda.^3)./(4*A*rho*R^3.*omega*Rs);
% figure()
% hold on; box on; grid on;
% for i=1:length(lambda)
%     plot(omega, cp_max_gen(i, :), 'LineWidth', line_width, 'DisplayName',['$\lambda$ = ',num2str(lambda(i))])
% end
% yline(cp_max, '--', 'DisplayName', '$c_{P, max}^{aero}$', 'LineWidth', line_width)
% ylim([0,1])
% xlabel('$\omega$ [rad/s]')
% ylabel('$c_P$ [-]')
% legend('Location','northeast')

% static electro curves
fig_static_electro_power = figure('Color','w');
hold on
plot(V0_v, P_f/1e6, 'LineWidth', line_width, 'DisplayName', 'Mech.', 'Color', colors_vect(1,:));
plot(V0_v, P_electro/1e6, 'LineWidth', line_width, 'DisplayName', 'Electro $\eta$=1', 'Color', colors_vect(3,:));
plot(V0_v, P_joule/1e6, '--','LineWidth', line_width, 'DisplayName', 'Joule loss $\eta$=1', 'Color', colors_vect(3,:));
plot(V0_v, P_electro_eta/1e6, 'LineWidth', line_width, 'DisplayName', ['Electro $\eta$=', num2str(generator.eta)], 'Color', colors_vect(2,:));
plot(V0_v, P_joule_eta/1e6, '--', 'LineWidth', line_width, 'DisplayName', ['Joule loss $\eta$=', num2str(generator.eta)], 'Color', colors_vect(2,:));
legend('Location', 'northwest');
xlabel('$V_{0}$ [m/s]')
ylabel('P [MW]')
grid on
box on
% export_figure(fig_static_electro_power, '\fig_static_electro_power.eps', path_images);

% % optimized electro curve
% figure()
% hold on
% for i=1:length(V0)
%   plot(omega, Pg_out(:, i)/1e6, 'DisplayName', ['$V_{0}$=', num2str(V0(i))])
% end
% plot(omega(max_p), max_v/1e6, '--', 'DisplayName','Static')
% grid on
% box on
% ylabel('$P_{G, out}$ [MW]')
% xlabel('$\omega$ [rad/s]')
% legend('location', 'eastoutside')

fig_generator_new_map = figure('Color', 'w');
hold on
grid on
plot(V0_v, P_f/1e6,'LineWidth', line_width, 'DisplayName', 'Mechanical', 'Color', colors_vect(1,:))
plot(V0_b, P_GE_b/1e6, 'LineWidth', line_width, 'DisplayName', 'Electro with control', 'Color', colors_vect(2,:))
plot(V0_a, P_GE_a/1e6, 'LineWidth', line_width, 'HandleVisibility','off', 'Color', colors_vect(2,:))
plot(V0_v, P_electro_eta/1e6, 'LineWidth', line_width, 'DisplayName', ['Electro $\eta$=', num2str(generator.eta)], 'Color', colors_vect(3,:));
grid on
box on
xlabel('$V_0$ [m/s]')
ylabel('P [MW]')
legend('location', 'southeast')
title('Powers')
if simulation.print_figure == 1
  export_figure(fig_generator_new_map, '\fig_generator_new_map.eps', path_images);
end
% Power coefficients as function of V0
% figure()
% hold on
% for i=1:length(V0)
%   plot(V0(i), cp(max_p(i), i), 'o', 'Color','b')
% end
% plot(V0_v, cP_f, '--')
% xlabel('V0')
% ylabel('$c_P$ [-]')

% Power as function of the pitch parameteric in V0
% figure()
% hold on
% for i=1:length(V0_a)
%   plot(theta_p*180/pi, Pg_out_a(:, i)/1e6, 'DisplayName', ['$V_{0}$=', num2str(V0_a(i))])
% end
% plot(theta_p(max_p_a)*180/pi, P_GE_a/1e6)
% yline(rotor.P_rated/1e6, '--r')
% xlabel('V0')
% ylabel('P [MW]')

% Pitch angle as function of V0
fig_new_pitch_map = figure('Color','w');
hold on
plot(V0_v, feather*180/pi, 'LineWidth', line_width, 'DisplayName', 'Limit on the rotor')
plot(V0_a, feather_a*180/pi, 'LineWidth', line_width, 'DisplayName', 'Limit on the generator')
xlabel('$V_{0}$ [m/s]')
ylabel('Pitch angle [deg]')
title('Comparison of the angle for the pitching to feather')
legend('location', 'northwest')
grid on
box on
if simulation.print_figure == 1
  export_figure(fig_new_pitch_map, '\fig_new_pitch_map.eps', path_images);
end
%% Save the data in a lookup table
rated_values_P_GE(1) = lambda_GE_mean;
rated_values_P_GE(2) = cp_GE_mean;
save('lookup\rated_values_P_GE.mat', "rated_values_P_GE");

% pitch to feather
lookup_pitch_P_GE(1, :) = V0_a;
lookup_pitch_P_GE(2, :) = feather_a;
save('lookup\lookup_pitch_P_GE.mat', "lookup_pitch_P_GE");

% generator electrical power
lookup_P_GE(1, :) = [V0_b, V0_a];
lookup_P_GE(2, :) = [P_GE_b', P_GE_a'];
save('lookup\lookup_P_GE.mat', "lookup_P_GE");
