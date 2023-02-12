%% Transfer functions of the PMSM generator

clc
close all
s = tf('s'); % define s as complex variable (i.e. s = j*omega)
syms ki R2 R1 Giq_noR 

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
omega_pm = generator.omega_pm;    
omega_omegaBP = generator.omega_omegaBP;
n = gearbox.ratio;

Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2);  % generataor TF
Gc = 1/(1 + s*tau_c);                                      % power converter TF
Giq_noR = Yiq*Gc;

% tune current PID controller
opts = pidtuneOptions('PhaseMargin', iq_pm);
[Riq, info] = pidtune(Giq_noR, 'pi', iq_omegaBP, opts);
kiq_p = Riq.kp;
kiq_i = Riq.ki;
kiq_d = Riq.kd;
Riq = kiq_p + 1/s*kiq_i + s*kiq_d;  % Iq controller TF

% tune in the second way
clear s
syms s
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2);  % generataor TF
Gc = 1/(1 + s*tau_c);                                      % power converter TF
Giq_noR = Yiq*Gc;
pol = poles(Giq_noR); % poles of the TF
pol_abs = abs(pol)
tau_i = pol_abs(1);
R1 = ki*(1 + s*tau_i)/s;
GH = R1*Giq_noR;
R2 = eval(norm(subs(GH, s, 1j*iq_omegaBP)));
ki_solve = eval(solve(R2 == 1, ki));
kp_solve = eval(tau_i*ki_solve);
s = tf('s');
R3 = kp_solve + ki_solve/s; 
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2);  % generataor TF
Gc = 1/(1 + s*tau_c);                                      % power converter TF
Giq_noR = Yiq*Gc;

figure(1)
bodeplot(Giq_noR, Riq*Giq_noR, R3*Giq_noR)
legend('Not controlled', 'Controlled', 'Manual','location', 'best')
grid on