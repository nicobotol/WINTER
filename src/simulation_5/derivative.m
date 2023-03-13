close all
clc

t_start = 0;
t_stop = 500;
dt = 5e-5;
t = t_start:dt:t_stop;
v_start = 4;
v_stop = 25;
v = v_start + (v_stop - v_start)/(t_stop - t_start)*t;
pitch = interp1(lookup_Pitch(1, :), lookup_Pitch(3, :), v)*180/pi;

deriv = gradient(pitch, dt);

figure()
plot(t, v)
hold on
plot(t, pitch)
xlabel('Time [s]')
ylabel('Wind speed [m/s]')

figure()
plot(t, deriv)
yline(blade.pitch_rate*180/pi)
xlabel('Time [s]')
ylabel('Pitch derivative [deg/s]')