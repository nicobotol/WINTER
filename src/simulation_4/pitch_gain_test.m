clear
close all
clc

parameters;

p_range = [0:1:30]*pi/180;
ki = polyval(blade.ki_schedule, p_range);
kp = polyval(blade.kp_schedule, p_range);

figure();
plot(p_range*180/pi, ki);
hold on
plot(p_range*180/pi, kp);
plot(blade.ki_tab(1,:), blade.ki_tab(2,:));
plot(blade.kp_tab(1,:), blade.kp_tab(2,:));
legend('Ki', 'Kp','Ki tab', 'Kp tab', 'location', 'best')
hold off

%%
V0 = 12:1:25;
pitch =  interp1(reference(:, 1), reference(:, 2), V0)*pi/180;
ki = polyval(blade.ki_schedule, pitch);
kp = polyval(blade.kp_schedule, pitch);

figure()
plot(V0, ki)
hold on
plot(V0, kp)
