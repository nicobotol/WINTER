function [Romega, ki, kp, kd, tau_d1] = omega_G_pid(design_method, enable_plot)
%% Tune a PID controller for the rotational speed
% G_cl -> close loop tf of the generator
parameters

% Redefine some parameters for clarity
omega_omegaBP = generator.omega_omegaBP;  % crossband frequency [rad/s]
omega_pm = generator.omega_pm;            % phase margin [rad]
B = B_eq;
I = I_eq;

%% Manual design of the controller
if design_method == 0
% Define the transfer function as symbol
syms s ki kd
G = 1/(s*I + B);                      % transfer function
[kp, ki, kd, tau_d1] = pid_tune(G, omega_omegaBP); % tune the gains

% Redefine the transfer function as 'transfer function' type 
s = tf('s');
Romega = (kp + ki/s + kd*s)/(1 + s*tau_d1);                   % regulator
GR = G*Romega;

%% Automatic design of the controller
elseif design_method == 1
s = tf('s');
G = 1/(s*I + B);            % transfer function
tau_d1 = 1/(10*omega_omegaBP);
opts = pidtuneOptions('PhaseMargin', omega_pm);
Romega_pid = pidtune(G, 'pid', omega_omegaBP, opts);
ki = Romega_pid.ki; % integral gain
kp = Romega_pid.kp; % proportional gain [1/s]
kd = Romega_pid.kd; % integral gain [s]
Romega = (Romega_pid.kp + Romega_pid.ki/s + s*Romega_pid.kd)/(1 + s*tau_d1);  
GR = G*Romega; % regulator

fprintf('ki = %f\n', Romega_pid.ki);
fprintf('kp = %f\n', Romega_pid.kp);

end

% Plot the phase margin
[~, PM] = margin(GR);
fprintf('Phase margin = %f [°]\n', PM);

% Close loop transfer function
G_cl = GR/(1 + GR);

%% Plot
if enable_plot == 1
  % Open loop, with regulator, regulator
  TF = [G_cl, GR, Romega];                                      % TF to plot
  legends = {'Not regulated', 'Regulated', 'Regulator'};  % legends names
  bode_plot(TF, legends, 'Open loop Bode plot', 'fig_bode_generator')

  % Close loop
  bode_plot(G_cl, {'Close loop'}, 'Close loop Bode plot', 'fig_bode_cl')

end

end