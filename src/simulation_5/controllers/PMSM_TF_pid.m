function [Yiq, Gc, Riq, GR, G_cl, generator] = ...
  PMSM_TF_pid(design_method, enable_plot)
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
Yiq = -(B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generataor TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
[kp, ki, kd, tau_d1] = pid_tune(G, iq_omegaBP);           % tune the gains

% Redefine the transfer function as 'transfer function' type 
s = tf('s');
Yiq = -(B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
Riq = -(kp + ki/s + kd*s)/(1 + s*tau_d1);                   % regulator
GR = G*Riq;

generator.kp = kp;
generator.ki = ki;
generator.kd = kd;
generator.tau_d1 = tau_d1;

% Save the gains as latex macro
names = ["GenkpMacroMan", "GenkiMacroMan", "GenkdMacroMan", "GentaudOneMacroMan", "GenMarginMan"];

%% Automatic design of the controller
elseif design_method == 1
  s = tf('s');
  Yiq = -(B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
  Gc = 1/(1 + s*tau_c);                                  % power converter TF
  G = Yiq*Gc;  
  opts = pidtuneOptions('PhaseMargin', generator.iq_pm);
  Riq_pid = pidtune(G, 'pid', iq_omegaBP, opts);
  kp = Riq_pid.kp;
  ki = Riq_pid.ki;
  kd = Riq_pid.kd;
  tau_d1 = 1/(10*iq_omegaBP);
  
  Riq = (kp + ki/s + kd*s)/(1 + s*tau_d1);  
  GR = G*Riq; % regulator
  
  kp = -kp;
  ki = -ki;
  kd = -kd;

  fprintf('ki = %f\n', Riq_pid.ki);
  fprintf('kp = %f\n', Riq_pid.kp);
  fprintf('kd = %f\n', Riq_pid.kd);
  
  generator.kp = Riq_pid.kp;
  generator.ki = Riq_pid.ki;
  generator.kd = Riq_pid.kd;
  generator.tau_d1 = tau_d1;

  names = ["GenkpMacroAuto", "GenkiMacroAuto", "GenkdMacroAuto", "GentaudOneMacroAuto", "GenMarginAuto"];
end

% Plot the phase margin
[~, PM] = margin(GR);
fprintf('Phase margin = %f [°]\n', PM);

% Close loop transfer function
G_cl = GR/(1 + GR);
% num = G_cl.Numerator{1,1}(end-3:end-2);
% den = G_cl.Denominator{1,1}(end-4:end-2);
% num = G_cl.Numerator{1,1}(1:end-2);
% den = G_cl.Denominator{1,1}(1:end-2);
% num = G_cl.Numerator{1,1}(9);
% den = G_cl.Denominator{1,1}(7:9);
% G_cl_simply = tf(num, den)
% 
% pole(G_cl)
% pole(G_cl_simply)
% zero(G_cl)
% zero(G_cl_simply)

if B_eq == 0
  generator.ki = 79.883578;
  generator.kp = 1.484989;
  generator.kd = 0.000723;
  generator.tau_d1 = 1/(1e1*generator.iq_omegaBP);
end

%% Plot
if enable_plot == 1
  % Open loop, with regulator, regulator
  TF = [G, GR, Riq];                                      % TF to plot
  legends = {'Not regulated', 'Regulated', 'Regulator'};  % legends names
  bode_plot(TF, legends, 'Open loop Bode plot', 'fig_bode_generator')
  
  % Close loop
  %   bode_plot([G_cl, G_cl_simply], {'Close loop', 'Simplified'}, 'Close loop Bode plot', 'fig_bode_cl')
  bode_plot(G_cl, {'Close loop'}, 'Close loop Bode plot', 'fig_bode_cl')
end

file_name = "C:\Users\Niccolò\Documents\UNIVERSITA\TESI_MAGISTRALE\report\macro\macro.tex";
mode = "a";
data = [kp, ki, kd, tau_d1, PM];
write_latex_macro(file_name, names, data, mode)


end