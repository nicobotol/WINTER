function [Yiq, Gc, Riq, GR,G_cl,generator] = PMSM_TF_pi(design_method, enable_plot)
%% Tune a PI controller for the generator
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
syms s ki
Yiq = -(B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generataor TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
[ki, kp] = pi_tune(G, iq_omegaBP);                        % tune the gains

% Redefine the transfer function as 'transfer function' type 
s = tf('s');
Yiq = -(B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
Riq = -(kp + ki/s);                                        % regulator
GR = G*Riq;

generator.kp = kp;
generator.ki = ki;
generator.kd = 0;
generator.tau_d1 = 0;

%% Automatic design of the controller
elseif design_method == 1
s = tf('s');
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
opts = pidtuneOptions('PhaseMargin', generator.iq_pm);
Riq_pid = pidtune(G, 'pi', iq_omegaBP, opts);
GR = G*Riq_pid;
Riq = (Riq_pid.kp + Riq_pid.ki/s);                     % regulator



fprintf('ki = %f\n', Riq_pid.ki);
fprintf('kp = %f\n', Riq_pid.kp);

generator.kp = Riq_pid.kp;
generator.ki = Riq_pid.ki;
generator.kd = 0;
generator.tau_d1 = 0;

end

% Plot the phase margin
[~, PM] = margin(GR)
% Close loop transfer function
G_cl = GR/(1 + GR);

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

end