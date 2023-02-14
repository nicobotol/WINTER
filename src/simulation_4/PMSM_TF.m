function [Yiq, Gc, Riq] = PMSM_TF(design_method, enable_plot, name)
parameters

% Redefine the parameters for clarity
B = generator.B;
I = generator.I;
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
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generataor TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
[ki, kp] = pi_tune(G, iq_omegaBP);                        % tune the gains

% Redefine the transfer function as 'transfer function' type 
s = tf('s');
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
Riq = (kp + ki/s);                                        % regulator
GR = G*Riq;

%% Automatic design of the controller
elseif design_method == 1
s = tf('s');
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
opts = pidtuneOptions('PhaseMargin', 45);
Riq_pid = pidtune(G, 'pi', iq_omegaBP, opts);
GR = G*Riq_pid;
Riq = (Riq_pid.kp + Riq_pid.ki/s);                     % regulator
end

%% Plot
if enable_plot == 1
[magG, phaseG, woutG] = bode(G);
[magGR, phaseGR, woutGR] = bode(GR);
[magRiq, phaseRiq, woutRiq] = bode(Riq);
fig_bode_gnerator = figure('Position', get(0, 'Screensize'));
subplot(2,1,1)
semilogx(woutG, 20*log10(squeeze(magG)), 'LineWidth', line_width)
hold on
semilogx(woutGR, 20*log10(squeeze(magGR)), 'LineWidth', line_width)
semilogx(woutRiq, 20*log10(squeeze(magRiq)), 'LineWidth', line_width)
hold off
grid
subplot(2,1,2)
semilogx(woutG, squeeze(phaseG), 'LineWidth', line_width)
hold on
semilogx(woutGR, squeeze(phaseGR), 'LineWidth', line_width)
semilogx(woutRiq, squeeze(phaseRiq), 'LineWidth', line_width)
hold off
grid
sgtitle('PSMS Bode plot');
legend('Open loop', 'With regulator', 'Regulator')
end
end