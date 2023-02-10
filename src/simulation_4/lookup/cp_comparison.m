clear;
close all;
clc;

% load the parameters
parameters
addpath('\..')
addpath("lookup");

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
colors = ['b', 'r', 'k', 'c', 'm', 'g'];
figure(1)
% result from the look up table
for p = 1:pitch_i
  cP_interp = interp2(lambda_vector, pitch_vector, lookup_cp, ...
    lambda_v(:), pitch_v(p));
  cT_interp = interp2(lambda_vector, pitch_vector, lookup_cT, ...
    lambda_v(:), pitch_v(p));

  % plot
  subplot(121)
  plot(lambda_v, cP_interp, '--', 'Color', colors(p), ...
    'LineWidth', line_width);
  hold on
  subplot(122)
  plot(lambda_v, cT_interp, '--', 'Color', colors(p), ...
    'LineWidth', line_width);
  hold on
  LegS{p} = ['\theta {',num2str(rad2deg(pitch_v(p))),'°}'];
end
% results from the pubblication
for p = 1:pitch_i
    subplot(121)
    plot(lambda_t, table_cP(p, :), 'o', 'Color', colors(p), ...
      'MarkerSize', marker_size, 'LineWidth', line_width);
    hold on
    subplot(122)
    plot(lambda_t, table_cT(p, :), 'o', 'Color', colors(p), ...
      'MarkerSize', marker_size, 'LineWidth', line_width);
    hold on
end
legend(LegS, 'location', 'best');
subplot(121)
xlabel('\lambda [-]')
ylabel('c_P [-]')
grid on
subplot(122)
xlabel('\lambda [-]')
ylabel('c_T [-]')
grid on
%% Compare the results from the analytical formula

cP = zeros(pitch_i, lambda_i);
for l = 1:lambda_i
  for p = 1:pitch_i
    lambda = lambda_v(l);
    pitch = pitch_v(p)*180/pi; % [°]
    lambda1 = 1/(1/(lambda + 0.089) - 0.035/(pitch^3 + 1));
    cP(p, l) = 0.5*(116/lambda1 - 0.4*(pitch - 5))*exp(-16.5/lambda1);
  end
end

figure(3)
for i = 1:pitch_i
  plot(lambda_v, cP(i, :));
  hold on
end

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

figure(2)
subplot(121)
plot(V0_v, cP_v, 'LineWidth', line_width)
hold on 
plot(V0_v, cP_ref, 'o', 'LineWidth', line_width)
grid on
xlabel('Wind speed [m/s]')
ylabel('cP [-]')
subplot(122)
plot(V0_v, cT_v, 'LineWidth', line_width)
hold on 
plot(V0_v, cT_ref, 'o', 'LineWidth', line_width)
xlabel('Wind speed [m/s]')
ylabel('cT [-]')
legend('Computed', 'DTU report', 'location', 'best')
grid on
