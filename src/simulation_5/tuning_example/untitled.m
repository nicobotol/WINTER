clear 
close all
clc
s = tf('s');

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

Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2);  % generator TF
Gc = 1/(1 + s*tau_c);                                      % power converter TF
Giq_noR = Yiq*Gc;

evalfr(Giq_noR, iq_omegaBP)
syms s ki ki1
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2);  % generator TF
Gc = 1/(1 + s*tau_c);                                      % power converter TF
Giq_noR = Yiq*Gc;
tau_reg = 1/eval(max(abs(poles(Giq_noR, s))));
Reg = ki*(1 + s*tau_reg);
GH = Reg*Giq_noR;
GH_mag = abs(subs(GH, s, 1j*iq_omegaBP));
ki = eval(solve(GH_mag == 1, ki))
kp = ki*tau_reg

% tau_reg1 = 1/eval(min(abs(poles(Giq_noR, s))));
% Reg1 = ki1*(1 + s*tau_reg1);
% GH1 = Reg1*Giq_noR;
% GH_mag1 = abs(subs(GH1, s, 1j*iq_omegaBP));
% ki1 = eval(solve(GH_mag1 == 1, ki1));
% kp1 = ki1*tau_reg1;


s = tf('s'); 
Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2);  % generator TF
Gc = 1/(1 + s*tau_c);                                      % power converter TF
Giq_noR = Yiq*Gc;
Reg = kp + ki/s;
GH = Reg*Giq_noR;
% 
% Reg1 = kp1 + ki1/s;
% GH1 = Reg1*Giq_noR;

Rauto = pidtune(Giq_noR, 'pi');


bodeplot(Giq_noR, GH, Rauto*Giq_noR, Reg)
legend('Not controlled', 'Manual', 'Auto', 'Controller')
grid on
