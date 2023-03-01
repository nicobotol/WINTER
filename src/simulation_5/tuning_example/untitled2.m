clear 
close all
clc

syms s ki
% omega = 500*2*pi;
% I = 0.5;
% tau_c = 50e-6;
% L = 5e-3;
% R = 0.5;
% B = 0.191;
% p = 4;
% Lambda = 0.205;
% iq_omegaBP = 2*500*pi;
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
omega_pm = generator.omega_pm;    
omega_omegaBP = generator.omega_omegaBP;
n = gearbox.ratio;
omega = 10000;

%% Manual design of the controller
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generataor TF
Gc = 1/(1 + s*tau_c);                                     % power converter TF
G = Yiq*Gc;  
% pol = poles(G, s);
% tau = eval(1/abs(max(pol)))
% GH = G*(1 + s*tau)*ki/s;
% GH_mag = abs(subs(GH, s, 1j*omega));
% ki = solve(GH_mag == 1, ki );
% kp = tau*ki;
% kp_s = eval(kp);
% ki_s = eval(ki);
[ki_s, kp_s] = pi_tune(G, omega);

s = tf('s');
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generator TF
Gc = 1/(1 + s*tau_c);                                     % power converter TF
G = Yiq*Gc;  
Reg = (kp_s + ki_s/s);
GR = G*Reg;
% phase [deg]
phase = atan2d(imag(evalfr(GR, omega)),real(evalfr(GR, omega)))

%% Automatic design of the controller
opts = pidtuneOptions('PhaseMargin', 45);
Rauto1 = pidtune(G, 'pi', omega, opts);
GRauto1 = G*Rauto1;
Rauto2 = pidtune(G, 'pi', omega);
GRauto2 = G*Rauto2;

bodeplot(G, GR, GRauto1, GRauto2, Reg)
grid on
legend('No R', 'Manual', 'Auto1', 'Auto2', 'Regulator')