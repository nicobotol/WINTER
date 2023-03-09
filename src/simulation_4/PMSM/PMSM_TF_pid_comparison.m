function [Yiq, Gc, Riq_manual, GR_manual] = PMSM_TF_pid_comparison(design_method, enable_plot)
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

% Define the transfer function as symbol
syms s ki kd
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generataor TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
[kp, ki, kd, tau_d1] = pid_tune_geometrical_mean(G, iq_omegaBP);           % tune the gains

% Redefine the transfer function as 'transfer function' type 
s = tf('s');
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
Riq_manual = (kp + ki/s + kd*s)/(1 + s*tau_d1);                   % regulator
GR_manual = G*Riq_manual;

%% Automatic design of the controller

s = tf('s');
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                  % power converter TF
G = Yiq*Gc;  
opts = pidtuneOptions('PhaseMargin', generator.iq_pm);
Riq_pid = pidtune(G, 'pid', iq_omegaBP, opts);
Riq_auto = (Riq_pid.kp + Riq_pid.ki/s + Riq_pid.kd*s)/(1 + s/(10*iq_omegaBP));  
GR_auto = G*Riq_auto; % regulator

fprintf('ki = %f\n', Riq_pid.ki);
fprintf('kp = %f\n', Riq_pid.kp);
fprintf('kd = %f\n', Riq_pid.kd);


%% Plot
if enable_plot == 1
[magG, phaseG, woutG] = bode(G);
[magGR_manual, phaseGR_manual, woutGR_manual] = bode(GR_manual);
[magGR_auto, phaseGR_auto, woutGR_auto] = bode(GR_auto);
[magRiq_manual, phaseRiq_manual, woutRiq_manual] = bode(Riq_manual);
[magRiq_auto, phaseRiq_auto, woutRiq_auto] = bode(Riq_auto);
fig_bode_generator = figure('Position', get(0, 'Screensize'), 'Color','w');
% subplot(2,1,1)
semilogx(woutG, 20*log10(squeeze(magG)), 'LineWidth', line_width)
hold on
semilogx(woutGR_manual, 20*log10(squeeze(magGR_manual)), 'LineWidth', line_width)
semilogx(woutRiq_manual, 20*log10(squeeze(magRiq_manual)), 'LineWidth', line_width)
semilogx(woutGR_auto, 20*log10(squeeze(magGR_auto)), 'LineWidth', line_width)
semilogx(woutRiq_auto, 20*log10(squeeze(magRiq_auto)), 'LineWidth', line_width)
hold off
xlabel('$\omega$ [rad/s]', 'FontSize', font_size, 'interpreter','latex')
ylabel('Mag. [dB]', 'FontSize', font_size, 'interpreter','latex')
grid on
set(gca, 'FontSize', font_size)
legend('Open loop', 'With regulator manual', 'Regulator manual', ...
  'With regulator auto', 'Regulator auto','interpreter','latex',...
  'FontSize', font_size, 'location', 'best')
% subplot(2,1,2)
% semilogx(woutG, squeeze(phaseG), 'LineWidth', line_width)
% hold on
% semilogx(woutGR_manual, squeeze(phaseGR_manual), 'LineWidth', line_width)
% semilogx(woutRiq_manual, squeeze(phaseRiq_manual), 'LineWidth', line_width)
% semilogx(woutGR_auto, squeeze(phaseGR_auto), 'LineWidth', line_width)
% semilogx(woutRiq_auto, squeeze(phaseRiq_auto), 'LineWidth', line_width)
% hold off
% grid on
% xlabel('$\omega$ \ [rad/s]', 'FontSize', font_size, 'interpreter','latex')
% ylabel('Phase [deg.]', 'FontSize', font_size, 'interpreter','latex')
% sgtitle('PMSM Bode plot', 'FontSize', font_size);
set(gca, 'FontSize', font_size)
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\fig_bode_generator','.png');
  export_fig('fig_bode_generator', fig_name);
end
end

end