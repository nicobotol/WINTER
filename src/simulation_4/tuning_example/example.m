clear
% close all
clc
font_size =14;
line_width = 2;
syms s ki

G =  (1 + s)/(1+s/1000)/(1+s/100)/(1+s/10);
omega = 1500;

pol = poles(G, s);                                  % pole of the TF
pol_abs = abs(pol);
pol_abs_sort = sort(pol_abs, 'ascend');
tau_p = 1/pol_abs_sort(1);                          % 1st zero
tau_d = 1/pol_abs_sort(3);                        % 2nd zero
R = ki/s*(1 + s*tau_p)*(1 + s*tau_d);  % regulator
GH = G*R;                           % series of regulator and system
GH_mag = abs(subs(GH, s, 1j*omega));% magnitude at the crossover frequency
ki = eval(solve(GH_mag == 1, ki));        % integral gain
kp = eval(ki*(tau_d + tau_p));
kd = eval(ki*tau_d*tau_p);

clear s;
s = tf('s');
G =  (1 + s)/(1+s/1000)/(1+s/100)/(1+s/10);
Riq = kp + ki/s + kd*s;
GR = G*Riq;


[magG, phaseG, woutG] = bode(G);
[magGR, phaseGR, woutGR] = bode(GR);
[magRiq, phaseRiq, woutRiq] = bode(Riq);
fig_bode_generator = figure('Position', get(0, 'Screensize'));
subplot(2,1,1)
xlabel('\omega [rad/s]', 'FontSize', font_size, 'interpreter','tex')
ylabel('Mag. [dB]', 'FontSize', font_size, 'interpreter','latex')
semilogx(woutG, 20*log10(squeeze(magG)), 'LineWidth', line_width)
hold on
semilogx(woutGR, 20*log10(squeeze(magGR)), 'LineWidth', line_width)
semilogx(woutRiq, 20*log10(squeeze(magRiq)), 'LineWidth', line_width)
hold off
grid on
subplot(2,1,2)
semilogx(woutG, squeeze(phaseG), 'LineWidth', line_width)
hold on
semilogx(woutGR, squeeze(phaseGR), 'LineWidth', line_width)
semilogx(woutRiq, squeeze(phaseRiq), 'LineWidth', line_width)
hold off
grid
xlabel('\omega [rad/s]', 'FontSize', font_size, 'interpreter','tex')
ylabel('Phase [deg.]', 'FontSize', font_size, 'interpreter','latex')
sgtitle('PSMS Bode plot');

