clear;
close all;
clc;

% load the parameters
addpath('..\');                        
addpath('..\aerodynamic_functions');
addpath('..\aerodynamic_functions\airfoil_data')
parameters;                                       % load the parameters
load('lookup_cP_theta_lambda.mat');               % load the cP lookup 
load('lookup_cT_theta_lambda.mat');               % load the cT lookup 
load("rated_values.mat");
load('lookup_pitch.mat');
V0_rated = rated_values(1);                       % rated wind speed [m/s]
omega_rated = rated_values(2);                    % rated speed [rad/s]
lambda_opt = rated_values(4);                     % optimal TSR
cp_max = rated_values(5);                         % maximum cp

%% Compare cP in the loookup table with the results proposed by PoliMi
% redefine a range for the pitch and the lambda
pitch_r = [-3, 12]*pi/180; % [rad]
lambda_r = [2, 10];
pitch_i = 6;
lambda_i = 10;
pitch_v = linspace(pitch_r(1), pitch_r(2), pitch_i);
lambda_v = linspace(lambda_r(1), lambda_r(2), lambda_i);

% table read from the PoliMi pubblication
%           2   4    6    8    10
table_cP = [0.02 0.2  0.45 0.43 0.34;  % -3
            0.03 0.25 0.44 0.47 0.42;  % 0
            0.04 0.27 0.41 0.45 0.45;  % 3
            0.05 0.26 0.34 0.37 0.36;  % 6
            0.07 0.24 0.27 0.24 0.15;  % 9
            0.08 0.20 0.18 0.07 -0.15];% 12
%           2   4    6    8    10
table_cT = [0.12 0.37 0.72 0.95 1.15;  % -3
            0.12 0.37 0.65 0.85 0.95;  % 0
            0.12 0.37 0.55 0.68 0.78;  % 3
            0.12 0.32 0.48 0.50 0.52;  % 6
            0.12 0.28 0.32 0.28 0.23;  % 9
            0.12 0.23 0.21 0.16 -0.09];% 12
lambda_t = [2, 4, 6, 8, 10];
theta_t = [-3, 0, 3, 6, 9, 12];

LegS = cell(1, pitch_i);
cP_interp = zeros(lambda_i, 1);
cT_interp = zeros(lambda_i, 1);
fig_cP_cT_comp_polimi = figure('Color','w');
%% result from the look up table
for p = 1:pitch_i
  cP_interp = interp2(lambda_vector, pitch_vector, lookup_cP, ...
    lambda_v(:), pitch_v(p));
  cT_interp = interp2(lambda_vector, pitch_vector, lookup_cT, ...
    lambda_v(:), pitch_v(p));

  % plot
  subplot(121)
  plot(lambda_v, cP_interp, '--', 'Color', colors_vect(p, :), ...
    'LineWidth', line_width);
  hold on
  subplot(122)
  plot(lambda_v, cT_interp, '--', 'Color', colors_vect(p, :), ...
    'LineWidth', line_width);
  hold on
  LegS{p} = ['$\theta$ = ',num2str(rad2deg(pitch_v(p))),'[deg]'];
end
% results from the pubblication
for p = 1:pitch_i
    subplot(121)
    plot(lambda_t, table_cP(p, :), 'o', 'Color', colors_vect(p, :), ...
      'MarkerSize', marker_size, 'LineWidth', line_width);
    hold on
    subplot(122)
    plot(lambda_t, table_cT(p, :), 'o', 'Color', colors_vect(p, :), ...
      'MarkerSize', marker_size, 'LineWidth', line_width);
    hold on
end
legend(LegS, 'location', 'best');
subplot(121)
xlabel('$\lambda$ [-]')
ylabel('cP [-]')
set(gca, 'FontSize', font_size)
title('Power coefficient')
grid on
subplot(122)
xlabel('$\lambda$ [-]')
ylabel('cT [-]')
title('Thrust coefficient')
grid on
set(gca, 'FontSize', font_size)
export_figure(fig_cP_cT_comp_polimi, '\fig_cP_cT_comp_polimi.eps', path_images);

%% Comparison between the results and the DTU report pag.34
V0_v = 4:1:25;                  % windspeed [m/s]
cP_v = zeros(length(V0_v), 1);  % cP from the lookup
cP_ref = reference(:, 4);       % cP from the dtu report
cT_v = zeros(length(V0_v), 1);  % cT from the lookup 
cT_ref = reference(:, 5);       % cT from the dtu report

for v=1:length(V0_v)
  V0 = V0_v(v);
  if V0 < V0_rated
    lambda = lambda_opt;
    pitch = 0;
  else
    lambda = omega_rated*rotor.R/V0;
    pitch = interp1(lookup_Pitch(1, :), lookup_Pitch(3, :), V0);
  end

  cP_v(v) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, pitch);
  cT_v(v) = interp2(lambda_vector, pitch_vector, lookup_cT, lambda, pitch);
end

fig_cP_cT_comp = figure('Color','w');
subplot(121)
plot(V0_v, cP_v, 'LineWidth', line_width)
hold on 
plot(V0_v, cP_ref, 'o', 'LineWidth', line_width)
grid on
xlabel('Wind speed [m/s]')
ylabel('cP [-]')
title('Power coefficient')
set(gca, 'FontSize', font_size)
subplot(122)
plot(V0_v, cT_v, 'LineWidth', line_width)
hold on 
plot(V0_v, cT_ref, 'o', 'LineWidth', line_width)
xlabel('Wind speed [m/s]')
ylabel('cT [-]')
legend('Computed', 'DTU report', 'location', 'best')
title('Thrust coefficient')
grid on
set(gca, 'FontSize', font_size)
export_figure(fig_cP_cT_comp, '\fig_cP_cT_comp.eps', path_images);

%% Plot the power and thrust as functions of the wind speed
stall = zeros(length(V0_v), 1);
feather = zeros(length(V0_v), 1);
lambda = zeros(length(V0_v), 1);
cP_s = zeros(length(V0_v), 1);
cP_f = zeros(length(V0_v), 1);
cT_s = zeros(length(V0_v), 1);
cT_f = zeros(length(V0_v), 1);
P_s = zeros(length(V0_v), 1);
P_f = zeros(length(V0_v), 1);
T_s = zeros(length(V0_v), 1);
T_f = zeros(length(V0_v), 1);

% generate the vectors of angles for stall and feathering
for v=1:length(V0_v)
  V0 = V0_v(v);
  if V0 < V0_rated
    lambda(v) = lambda_opt;
  else
    lambda(v) = omega_rated*rotor.R/V0;
  end
end

omega = lambda.*V0_v'/rotor.R; % [rad/s] rotional speed

stall(:) = interp1(lookup_Pitch(1, :), lookup_Pitch(2, :), V0_v); % [rad]
feather(:) = interp1(lookup_Pitch(1, :), lookup_Pitch(3, :), V0_v);%[rad]

cP_s(:) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, stall);
cP_f(:) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, feather);
cT_s(:) = interp2(lambda_vector, pitch_vector, lookup_cT, lambda, stall);
cT_f(:) = interp2(lambda_vector, pitch_vector, lookup_cT, lambda, feather);

P_s(:) = 0.5*rotor.A*rho.*cP_s.*V0_v'.^3; % [W]
P_f(:) = 0.5*rotor.A*rho.*cP_f.*V0_v'.^3; % [W]
T_s(:) = 0.5*rotor.A*rho.*cT_s.*V0_v'.^2; % [N]
T_f(:) = 0.5*rotor.A*rho.*cT_f.*V0_v'.^2; % [N]

fig_power_vs_V0 = figure('Color','w');
plot(V0_v, P_s/1e6, 'LineWidth', line_width*2, 'Color', colors_vect(3,:));
hold on
plot(V0_v, P_f/1e6, '--', 'LineWidth', line_width, 'Color', colors_vect(1,:));
legend('Stall', 'Feather', 'Location', 'best');
grid on
xlabel('$V_0$ [m/s]')
ylabel('Power [MW]')
title('Power as function of the wind speed')
set(gca, 'FontSize', font_size)
export_figure(fig_power_vs_V0, '\fig_power_vs_V0.eps', path_images);

fig_thrust_vs_V0 = figure('Color','w');
plot(V0_v, T_s/1e6, 'LineWidth', line_width, 'Color', colors_vect(3,:));
hold on
plot(V0_v, T_f/1e6, '--', 'LineWidth', line_width, 'Color', colors_vect(1,:));
legend('Stall', 'Feather', 'Location', 'best');
grid on
xlabel('$V_0$ [m/s]')
ylabel('Thrust [MN]')
title('Thrust as function of the wind speed')
set(gca, 'FontSize', font_size)
export_figure(fig_thrust_vs_V0, '\fig_thrust_vs_V0.eps', path_images);

%% Steady state powers
iq = 2/3*P_s./(generator.Lambda.*omega*generator.p); % [A] generator current
uq =  omega*generator.p*generator.Lambda - generator.Rs*iq; % [V] generator voltage
P_electro = 1.5*uq.*iq; % [W] electrical power
P_joule = 1.5*generator.Rs*iq.^2; % [W] Joule losses
P_electro_eta = 1.5*uq.*iq*generator.eta; % [W] electrical power
P_joule_eta = 1.5*generator.Rs*(iq*generator.eta).^2; % [W] Joule losses

fig_static_electro_power = figure('Color','w');
hold on
plot(V0_v, P_s/1e6, 'LineWidth', line_width, 'DisplayName', 'Mech.', 'Color', colors_vect(1,:));
plot(V0_v, P_electro/1e6, 'LineWidth', line_width, 'DisplayName', 'Electro $\eta$=1', 'Color', colors_vect(3,:));
plot(V0_v, P_joule/1e6, '--','LineWidth', line_width, 'DisplayName', 'Joule loss $\eta$=1', 'Color', colors_vect(3,:));
plot(V0_v, P_electro_eta/1e6, 'LineWidth', line_width, 'DisplayName', ['Electro $\eta$=', num2str(generator.eta)], 'Color', colors_vect(2,:));
plot(V0_v, P_joule_eta/1e6, '--', 'LineWidth', line_width, 'DisplayName', ['Joule loss $\eta$=', num2str(generator.eta)], 'Color', colors_vect(2,:));
legend('Location', 'northwest');
xlabel('$V_0$ [m/s]')
ylabel('P [MW]')
grid on
box on
export_figure(fig_static_electro_power, '\fig_static_electro_power.eps', path_images);
