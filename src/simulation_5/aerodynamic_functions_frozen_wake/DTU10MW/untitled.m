I_eq = 1.56e8;
omega_zeta = 0.6;
zeta = 0.7;
n = 50;
dTdtheta0 = 2.72e6;
kk = 9; 
theta_v = 0:0.25:22;
GK = (1./(1 + theta_v./kk));
ki = I_eq*omega_zeta^2/(180/pi*dTdtheta0).*GK;
kp = 2*I_eq*zeta*omega_zeta/(180/pi*dTdtheta0)*GK;
ki_val = polyval(blade.ki_schedule_report, theta_v*pi/180);
kp_val = polyval(blade.kp_schedule_report, theta_v*pi/180);

figure(10); clf;
plot(theta_v, GK)
hold on
plot(theta_v, ki)
plot(theta_v, ki_val)
plot(theta_v, kp, '--')
plot(theta_v, kp_val, '--')
legend('GK', 'ki', ' ki ref', 'kp', 'kp ref')