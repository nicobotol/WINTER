clear 
close all
clc

parameters;
Rs = generator.Rs; % [Ohm] stator resistance
A = rotor.A; % [m^2] rotor swept area
B = B_eq; % [kgm^2] transmission equiuvalent damping
p = generator.p; % [-] number of pole pairs
Lambda = generator.Lambda; % [Wb] generator flux linkage
eta = generator.eta; % [-] generator efficiency
R = rotor.R; % [m]

%% Below rated wind speed
physical(1) = generator.Rs; % [Ohm] stator resistance
physical(2) = rotor.A; % [m^2] rotor swept area
physical(3) = B_eq; % [kgm^2] transmission equiuvalent damping
physical(4) = generator.p; % [-] number of pole pairs
physical(5) = generator.Lambda; % [Wb] generator flux linkage
physical(6) = generator.eta; % [-] generator efficiency
physical(7) = rotor.R; % [m]
physical(8) = rho; % [kg/m^3]

% parameters without the damping
physical_no_B = physical;
physical_no_B(3) = 0;

% V0_b = 4:0.01:6-0.01;
V0_b = [4:0.05:V0_rated+0.001]; % [m/s] wind speed 

% Maximize the power output by selecting the proper combination of TSR and pitch angle
% The variable x containts x = (TSR, pitch angle)
lb = [7.5, -5*pi/180]; % variable lower bound
ub = [10, 5*pi/180]; % variable upper bound
P = zeros(length(V0_b),1);
min_v = zeros(length(V0_b),2); % minimization values (lambda, theta)
x0 = [7.5, 0]; % initial guess value
min_v_no_B = zeros(length(V0_b),2); % minimization values (lambda, theta)
options = optimoptions(@fmincon,'Display', 'off');
for i=1:length(V0_b)
  V0 = V0_b(i);    
  [min_v(i, :), P(i, :)] = fmincon(@(x)compute_P_GE(x, physical, lambda_vector, pitch_vector, lookup_cP, V0), x0, [], [], [], [], lb, ub, [], options);
  [min_v_no_B(i, :), P_no_B(i, :)] = fmincon(@(x)compute_P_GE(x, physical_no_B, lambda_vector, pitch_vector, lookup_cP, V0), x0, [], [], [], [], lb, ub, [],options);
end

% Recompute some values
theta_opt = min_v(:, 2); % pitch coming from minimization
theta_opt_no_B = min_v_no_B(:, 2); % pitch coming from minimization
P_GE_b = -P; % power coming from minimization (change sign because negative)
P_GE_b_no_B = -P_no_B; % power coming from minimization (change sign because negative)
omega_rotor = lambda_opt.*V0_b/R; % [rad/s] omega based on the maximization of the rotor power
omega_GE = min_v(:,1)'.*V0_b/R;
omega_GE_no_B = min_v_no_B(:,1)'.*V0_b/R;
lambda_GE =  min_v(:, 1);  % optimal TSR based on generator maximization 
lambda_GE_no_B =  min_v_no_B(:, 1);  % optimal TSR based on generator maximization 
lambda_GE_mean = mean(lambda_GE);  % mean TSR based on the generator
lambda_GE_mean_no_B = mean(lambda_GE_no_B);  % mean TSR based on the generator
cp_GE = interp2(lambda_vector, pitch_vector, lookup_cP, min_v(:, 1), min_v(:, 2)); % cp for the maximization of the generator
cp_GE_no_B = interp2(lambda_vector, pitch_vector, lookup_cP, min_v_no_B(:, 1), min_v_no_B(:, 2)); % cp for the maximization of the generator
cp_GE_mean = mean(cp_GE);          % mean cp for the generator maximization
cp_GE_mean_no_B = mean(cp_GE_no_B);          % mean cp for the generator maximization
cp_E = P_GE_b./(0.5*A*rho*V0_b'.^3); % electrical generator power normalized by the wind one

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
fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, omega_GE, 'LineWidth',line_width, 'DisplayName','Gen.')
plot(V0_b, omega_GE_no_B, 'LineWidth',line_width, 'DisplayName','Gen. B=0 $[\frac{kgm^2}{s}]$')
plot(V0_b, omega_rotor,'--', 'LineWidth',line_width, 'DisplayName','Rotor', 'Color', colors_vect(5,:))
title('Rotor rotational speed')
xlabel('$V_0$ [m/s]')
ylabel('$\omega$ [rad/s]')
legend('location', 'northwest')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_omega_GE', '.eps'), path_images); 
end

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, theta_opt*180/pi, 'LineWidth',line_width, 'DisplayName','Gen.')
plot(V0_b, theta_opt_no_B*180/pi, 'LineWidth',line_width, 'DisplayName','Gen. B=0$[\frac{kgm^2}{s}]$')
yline(0, '--', 'Color', colors_vect(5,:),'LineWidth',line_width, 'FontSize', font_size, 'DisplayName','Rotor')
title('Pitch angle')
xlabel('$V_0$ [m/s]')
ylabel('$\theta$ [deg]')
legend('location','east')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_pitch_GE', '.eps'), path_images); 
end

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, lambda_GE, 'LineWidth',line_width, 'DisplayName','Gen.','Color',colors_vect(1,:))
plot(V0_b, lambda_GE_no_B, 'LineWidth',line_width, 'DisplayName','Gen. B=0$[\frac{kgm^2}{s}]$','Color',colors_vect(2,:))
yl_1 = yline(lambda_GE_mean, '--', 'Mean', 'LineWidth',line_width, 'HandleVisibility', 'off', 'Color',colors_vect(1,:), 'FontSize', font_size);
yl_1.LabelHorizontalAlignment = 'left';
yl_2 = yline(lambda_GE_mean_no_B, '--', 'Mean', 'LineWidth',line_width, 'HandleVisibility', 'off', 'Color',colors_vect(2,:), 'FontSize', font_size);
yl_2.LabelHorizontalAlignment = 'left';
yl_3 = yline(lambda_opt, '--', 'Rotor','LineWidth',line_width, 'HandleVisibility', 'off','Color',colors_vect(5,:), 'FontSize', font_size);
yl_3.LabelHorizontalAlignment = 'left';
title('Tip speed ratio')
xlabel('$V_0$ [m/s]')
ylabel('$\lambda$ [-]')
legend('location', 'east')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_lambda_GE', '.eps'), path_images); 
end

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, cp_GE, 'LineWidth',line_width, 'DisplayName','Gen.')
plot(V0_b, cp_GE_no_B, 'LineWidth',line_width, 'DisplayName','Gen. B=0$[\frac{kgm^2}{s}]$')
y1 = yline(cp_GE_mean, '--', 'Mean', 'LineWidth',line_width, 'HandleVisibility','off','Color',colors_vect(1,:), 'FontSize', font_size);
y2 = yline(cp_GE_mean_no_B, '--', 'Mean', 'LineWidth',line_width,'HandleVisibility','off','Color',colors_vect(2,:), 'FontSize', font_size);
y1.LabelHorizontalAlignment = 'left';
y1.LabelVerticalAlignment = 'bottom';
y2.LabelHorizontalAlignment = 'left';
y2.LabelVerticalAlignment = 'bottom';
yline(cp_max, '--', 'LineWidth',line_width, 'DisplayName','Rotor','Color',colors_vect(5,:))
title('Power coefficient')
xlabel('$V_0$ [m/s]')
ylabel('$c_P$ [-]')
legend('location', 'southeast')
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
% ylabel('$P_{GE}$ [MW]')
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

% Pitch angle as function of V0
fig_cpE = figure('Color','w');
hold on
plot(V0_b, cp_E, 'LineWidth', line_width, 'DisplayName', '')
xlabel('$V_{0}$ [m/s]')
ylabel('$c_{P,E} [-]$')
title('Electrical power coefficient')
legend('location', 'northwest')
grid on
box on
if simulation.print_figure == 1
  export_figure(fig_cpE, '\fig_fig_cpE.eps', path_images);
end

%% Save the data in a lookup table
% clear("lookup_P_GE")
% clear("rated_values_P_GE")
% clear("lookup_pitch_P_GE")
% 
% rated_values_P_GE(1) = lambda_GE_mean;
% rated_values_P_GE(2) = cp_GE_mean;
% rated_values_P_GE(3) = lambda_GE_mean*V0_rated/rotor.R;
% save('lookup\rated_values_P_GE.mat', "rated_values_P_GE");
% 
% % pitch to feather
% lookup_pitch_P_GE = zeros(2, length(V0_a) + 1);
% lookup_pitch_P_GE(1, :) = [0 V0_a];
% lookup_pitch_P_GE(2, :) = [0 feather_a];
% save('lookup\lookup_pitch_P_GE.mat', "lookup_pitch_P_GE");
% 
% % generator electrical power
% lookup_P_GE(1, :) = [V0_b, V0_a];
% lookup_P_GE(2, :) = [P_GE_b', P_GE_a'];
% save('lookup\lookup_P_GE.mat', "lookup_P_GE");
