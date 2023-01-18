%% EXAMPLE OF CONTROL 
% This simulation is the most simple one: the power is based on a look-up
% table relating cp and tip speed ratio, no pitching is taken into account,
% the genreator model is discrete

clear 
close all
clc

%% Load parameters

% rotor parameters
rotor.R = 9.5;      % rotor diameter [m]
rotor.I = 15e3;     % rotor iniertia [kgm^2]
rotor.omega_r = 1;  % initial rotational speed [rad/s]

% generator parameters
generator.sincronous_velocity = 5;    % [rad/s]
generator.MG_coefficient = 2.55e5;    % generator curve coefficient

% general parameters
rho = 1.225;                % air density [kg/m^3]
stop_time = 300;            % max time to investigaste [s]

cP_vs_lambda = load('Cp_Lambda.txt');   % power coefficient's look-up table 
wind_speed = load('usim.dat');          % wind time history
sample_time = wind_speed(2,1) - wind_speed(1, 1); % WS sample time [s]

% R = 9.5;                    % rotor diameter [m]
% sincronous_velocity = 5;    % [rad/s]
% MG_coefficient = 2.55e5;    % generator curve coefficient
% I = 15e3;                   % rotor iniertia [kgm^2]
% omega_r = 1;                % initial condition rotational speed [rad/s]

%% Simulink simulation
mdl = 'winter_simulink';                        % model's name
open_system(mdl);                               % open the model
set_param(mdl, 'StopTime', num2str(stop_time)); % set simulation time
in = Simulink.SimulationInput(mdl);             % set simulation parameters
out = sim(in, 'ShowProgress','on');             % run the simulation

%% Plot of the results
figure()
plot(out.tout, out.Torques(:, 2))
hold on
plot(out.tout, out.Torques(:, 1))
legend('Rotor', 'Generator', 'Location','best')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on