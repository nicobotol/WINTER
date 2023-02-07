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

% reference values
omega_r_ref = 1;  % omega [rad]
Mr = 1e5;         % torque [Nm]
Iq_ref = 12e3; 
tA = 25;

Yq = (B + s*I)/(1.5*(p*Lambda)^2*(1 + 2/3*(R + s*L)*(B + s*I) ...
  /(p*Lambda)^2));                  % genertaor transfer function
Gc = 1/(1 + s*tau_c);               % power converter TF
G_noR = Yq*Gc;

opts = pidtuneOptions('PhaseMargin', 30);
[Riq] = pidtune(G_noR, 'pi', opts);
% [Riq] = pidtune(G_noR, 'pi');
% kiq_p = Riq.kp;
kiq_p = 20;
kiq_i = Riq.ki;
kiq_d = Riq.kd;
Riq = kiq_p + 1/s*kiq_i + s*kiq_d;  % controller TF
Riq_num = Riq.num{1};
Riq_den = Riq.den{1};



% [kiq_p, kiq_i, kiq_d] = PI_tuning(G_noR, tA);
% 
% Riq = kiq_p + 1/s*kiq_i + s*kiq_d;  % controller TF
% iq_controller_type = 0;
% if iq_controller_type == 1 % P
%   Riq_num = kiq_p;
%   Riq_den = 1;
% else                   % Pi or PID
%   Riq_num = Riq.num{1};
%   Riq_den = Riq.den{1};
% end
% 
% % speed control values
% komega_p = 50; % proportional gain
% komega_i = 100; % integral gain
% komega_d = 0;    % deruivative gain
% Romega = komega_p + 1/s*komega_i + s*komega_d;  % controller TF
% omega_controller_type = 0;
% if omega_controller_type == 1 % P
%   Romega_num = komega_p;
%   Romega_den = 1;
% else                   % Pi or PID
%   Romega_num = Romega.num{1};
%   Romega_den = Romega.den{1};
% end


%% Run the model
mdl = 'generator_simulink';                        % model's name
open_system(mdl);                               % open the model
set_param(mdl, 'StopTime', num2str(stop_time)); % set simulation time
in = Simulink.SimulationInput(mdl);             % set simulation parameters
out = sim(in, 'ShowProgress','on');             % run the simulation