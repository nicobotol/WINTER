clear 
close all
clc

%% Load parameters and open the model
parameters
mdl = 'generator_simulink';                     % model's name
open_system(mdl);                               % open the model

%% Define the transfer functions
s = tf('s'); % define s as complex variable

% generator parameters
B = generator.B;
I = generator.I;
p = generator.p;
Lambda = generator.Lambda;
R = generator.Rs;
L = generator.Ld;
tau_c =generator.tau_c;
iq_pm = generator.iq_pm;
iq_omegaBP = generator.iq_omegaBP;
omega_pm = generator.omega_pm;    % phase margin for the speed controller [°]
omega_omegaBP = generator.omega_omegaBP;% speed controller crossover frequency [rad/s]
n = gearbox.ratio;

% reference values
omega_r_ref = 1;  % omega [rad]
Mr = 1e5;         % torque [Nm]
Iq_ref = 5e3; 
tA = 0.5;

%% Iq controller
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2);  % generataor TF
Gc = 1/(1 + s*tau_c);               % power converter TF
Giq_noR = Yiq*Gc;                   % G2


% ST0 = slTuner('mdl', {Riq2, Romega2}); % set the regulator blocks
% addPoint(ST0, {omega_ref_signal, Iq, omega_r}); % define the signals of interest
% Req_omega = TuningGoal.LoopShape('omega_r', 0.2/s);
% Req_omega.Name = 'Omega loop';
% 
% Req_iq = TuningGoal.LoopShape('Iq', 2/s);
% Req_iq.Openings = 'omega_r';
% Req_iq.Name = 'Iq loop';
% 
% % tune the system and view the results
% ST = systune(ST0, [Req_omega, Req_iq]);
% showTunable(ST)
% viewGoal([Req_omega, Req_iq], ST)
%%
% tune current PID controller
opts = pidtuneOptions('PhaseMargin', iq_pm);
[Riq, info] = pidtune(Giq_noR, 'pi', iq_omegaBP, opts);
kiq_p = Riq.kp;
kiq_i = Riq.ki;
kiq_d = Riq.kd;
Riq = kiq_p + 1/s*kiq_i + s*kiq_d;  % Iq controller TF
info

W_iq = Giq_noR*Riq/(1 + Giq_noR*Riq); % close loop TF

%% Speed controller
Gomega_noR = 1.5*p*Lambda/(n*(s*I_eq + B_eq)); % G1

% tune speed PID controller
opts = pidtuneOptions('PhaseMargin', omega_pm);
[Romega, info] = pidtune(Giq_noR, 'pi', omega_omegaBP, opts);
komega_p = Romega.kp;
komega_i = Romega.ki;
komega_d = Romega.kd;
Romega = komega_p + 1/s*komega_i + s*komega_d;  % Iq controller TF
info

%% Run the model
set_param(mdl, 'StopTime', num2str(stop_time)); % set simulation time
in = Simulink.SimulationInput(mdl);             % set simulation parameters
out = sim(in, 'ShowProgress','on');             % run the simulation