clear 
close all
clc

parameters

s = tf('s'); % define s as complex variable

% generator parameters
B = generator.B;
I = generator.I;
p = generator.p;
Lambda = generator.Lambda;
R = generator.Rs;
L = generator.Lq;
tau_c = generator.tau_c;

omega_r_ref = 1;  % omega reference value [rad]

Mr = 1e5;         % torque reference [Nm]

% current control values
kiq_p = 50; % proportional gain
kiq_i = 100; % integral gain
kiq_d = 0;    % deruivative gain
Yq = (B + s*I)/(1.5*(p*Lambda)^2*(1 + 2/3*(R + s*L)*(B + s*I) ...
  /(p*Lambda)^2));            % genertaor transfer function
Gc = 1/(1 + s*tau_c);         % power converter TF
Riq = kiq_p + 1/s*kiq_i + s*kiq_d;  % controller TF
iq_controller_type = 0;
if iq_controller_type == 1 % P
  Riq_num = kiq_p;
  Riq_den = 1;
else                   % Pi or PID
  Riq_num = Riq.num{1};
  Riq_den = Riq.den{1};
end

% speed control values
komega_p = 50; % proportional gain
komega_i = 100; % integral gain
komega_d = 0;    % deruivative gain
Romega = komega_p + 1/s*komega_i + s*komega_d;  % controller TF
omega_controller_type = 0;
if omega_controller_type == 1 % P
  Romega_num = komega_p;
  Romega_den = 1;
else                   % Pi or PID
  Romega_num = Romega.num{1};
  Romega_den = Romega.den{1};
end



%% Run the model
mdl = 'generator_simulink';                        % model's name
open_system(mdl);                               % open the model
set_param(mdl, 'StopTime', num2str(stop_time)); % set simulation time
in = Simulink.SimulationInput(mdl);             % set simulation parameters
out = sim(in, 'ShowProgress','on');             % run the simulation