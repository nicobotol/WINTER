function [Yiq, Gc, Riq, GR] = PMSM_TF_pid(design_method, enable_plot)
%% Tune a PID controller for the generator
parameters

% Redefine the parameters for clarity
B = B_eq;
I = I_eq;
p = generator.p;
Lambda = generator.Lambda;
R = generator.Rs;
L = generator.Ld;
tau_c =generator.tau_c;
iq_pm = generator.iq_pm;
iq_omegaBP = generator.iq_omegaBP;

%% Manual design of the controller
if design_method == 0
% Define the transfer function as symbol
syms s ki kd
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generataor TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
[kp, ki, kd, tau_d1] = pid_tune(G, iq_omegaBP);           % tune the gains

% Redefine the transfer function as 'transfer function' type 
s = tf('s');
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
Riq = (kp + ki/s + kd*s)/(1 + s*tau_d1);                   % regulator
GR = G*Riq;

%% Automatic design of the controller
elseif design_method == 1
s = tf('s');
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
opts = pidtuneOptions('PhaseMargin', generator.iq_pm);
Riq_pid = pidtune(G, 'pid', iq_omegaBP, opts);
Riq = (Riq_pid.kp + Riq_pid.ki/s + Riq_pid.kd*s)/(1 + s/(10*iq_omegaBP));  
GR = G*Riq; % regulator

fprintf('ki = %f\n', Riq_pid.ki);
fprintf('kp = %f\n', Riq_pid.kp);
fprintf('kd = %f\n', Riq_pid.kd);
end

% Plot the phase margin
[~, PM] = margin(GR);
fprintf('Phase margin = %f [°]\n', PM);

% Close loop transfer function
G_cl = GR/(1 + GR);

%% Plot
if enable_plot == 1
  % Open loop, with regulator, regulator
  TF = [G, GR, Riq];                                      % TF to plot
  legends = {'Not regulated', 'Regulated', 'Regulator'};  % legends names
  bode_plot(TF, legends, 'Open loop Bode plot', 'fig_bode_generator')

  % Close loop
  bode_plot(G_cl, {'Close loop'}, 'Close loop Bode plot', 'fig_bode_cl')

end

end