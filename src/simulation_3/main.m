%% Simulation 3
% cP computed via an expression based on 10-MW_Direct-Drive...
% damping of the rotor dynamic taken into account

clear 
close all
clc

s = tf('s'); % define s a complex variable

%% Load parameters
parameters;

% transform the structs of parameters into buses for simulink
rotor_bus_info = Simulink.Bus.createObject(rotor); 
rotor_bus = evalin('base', rotor_bus_info.busName);
generator_bus_info = Simulink.Bus.createObject(generator); 
generator_bus = evalin('base', generator_bus_info.busName);
gearbox_bus_info = Simulink.Bus.createObject(gearbox); 
gearbox_bus = evalin('base', gearbox_bus_info.busName);

% wind time history
wind_speed = load('usim.dat');                    % [m/s] 
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