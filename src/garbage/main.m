%% Simulation 3
% cP computed via an expression based on 10-MW_Direct-Drive...
% damping of the rotor dynamic taken into account


clear 
close all
clc

%% Load parameters

% general parameters
rho = 1.225;                % air density [kg/m^3]
stop_time = 300;            % max time to investigaste [s]

% rotor parameters
rotor.R = 9.5;      % rotor diameter [m]
rotor.I = 15e3;     % rotor iniertia [kgm^2]
rotor.omega_r = 1;  % initial rotational speed [rad/s]
rotor.cP_c = [0.78, 151, 0.58, 0.002, 2.14, 13.2, 20.9, -0.002, -0.008];
rotor_bus_info = Simulink.Bus.createObject(rotor); % transform tha struct into a bus
rotor_bus = evalin('base', rotor_bus_info.busName);

% generator parameters
generator.sincronous_velocity = 5;    % [rad/s]
generator.MG_coefficient = 2.55e5;    % generator curve coefficient
generator.I = 1e3;                    % generator iniertia [kgm^2]
generator.vll = 4e3;                  % rated line-to-line voltage [V]
generator.is = 1443.4;                % rated stator current [A]
generator.fe = 26.66;                 % rated stator frequency [Hz]
genarator.poles = 320;                % number of poles
genarator.Ld = 1.8e-3;                % d-axis stator inductance [H]
generator.Lq = 1.8e-3;                % q-axis stator inductance [H]
generator.Rs = 64;                    % stator resistance [ohm]
generator.Lambda = 19.49;             % magnet flux-linkage [Wb]

% gearbox parameters
gearbox.ratio = 1;  % gearbox transmission ratio 

cP_vs_lambda = load('Cp_Lambda.txt');   % power coefficient's look-up table 
wind_speed = load('usim.dat');          % wind time history
sample_time = wind_speed(2,1) - wind_speed(1, 1); % WS sample time [s]

I_eq = rotor.I/gearbox.ratio^2 + generator.I; % equivalent inertia [kgm^2]

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