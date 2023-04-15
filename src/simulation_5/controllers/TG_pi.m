function [R_TG, GR] = TG_pi(design_method, enable_plot)
%% Tune a PI controller for the generator torque setpoint
parameters

% Redefine the parameters for clarity
B = B_eq;
I = I_eq;
TG_pm = generator.TG_pm;
TG_omegaBP = generator.TG_omegaBP;

%% Manual design of the controller
if design_method == 0
% Define the transfer function as symbol
syms s
G = 1/(B + s*I);
[ki, kp] = pi_tune(G, TG_omegaBP);  % tune the gains

% Redefine the transfer function as 'transfer function' type 
s = tf('s');
G = 1/(B + s*I); 
R_TG = (kp + ki/s);                                        % regulator
GR = G*R_TG;

%% Automatic design of the controller
elseif design_method == 1
s = tf('s');
G = 1/(B + s*I);  
opts = pidtuneOptions('PhaseMargin', TG_pm);
R_TG_pi = pidtune(G, 'pi', iq_omegaBP, opts);
GR = G*R_TG_pi;
R_TG = (R_TG_pi.kp + R_TG_pi.ki/s);                     % regulator

fprintf('ki = %f\n', R_TG_pi.ki);
fprintf('kp = %f\n', R_TG_pi.kp);

end

% Plot the phase margin
[~, PM] = margin(GR);
fprintf('Phase margin = %f\n', PM);

%% Plot
if enable_plot == 1
[magG, phaseG, woutG] = bode(G);
[magGR, phaseGR, woutGR] = bode(GR);
[magR_TG, phaseR_TG, woutR_TG] = bode(R_TG);
fig_bode_generator = figure('Position', get(0, 'Screensize'));
subplot(2,1,1)
xlabel('\omega [rad/s]', 'FontSize', font_size, 'interpreter','tex')
ylabel('Mag. [dB]', 'FontSize', font_size, 'interpreter','latex')
semilogx(woutG, 20*log10(squeeze(magG)), 'LineWidth', line_width)
hold on
semilogx(woutGR, 20*log10(squeeze(magGR)), 'LineWidth', line_width)
semilogx(woutR_TG, 20*log10(squeeze(magR_TG)), 'LineWidth', line_width)
hold off
grid on
subplot(2,1,2)
semilogx(woutG, squeeze(phaseG), 'LineWidth', line_width)
hold on
semilogx(woutGR, squeeze(phaseGR), 'LineWidth', line_width)
semilogx(woutR_TG, squeeze(phaseR_TG), 'LineWidth', line_width)
hold off
grid
xlabel('\omega [rad/s]', 'FontSize', font_size, 'interpreter','tex')
ylabel('Phase [°]', 'FontSize', font_size, 'interpreter','latex')
sgtitle('PSMS Bode plot');
legend('Open loop', 'With regulator', 'Regulator','interpreter','latex')
end

end