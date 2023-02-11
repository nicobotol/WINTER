clear 
close all
% s = tf('s');
parameters
syms s kp ki tau_iq % B I p Lambda R L tau_c 

% generator parameters
B = generator.B;
I = generator.I;
p = generator.p;
Lambda = generator.Lambda;
R = generator.Rs;
L = generator.Ld;
tau_c =generator.tau_c;
iq_pm = generator.iq_pm;
iq_omegaBP = generator.iq_omegaBP;
omega_pm = generator.omega_pm;    % phase margin for the speed controller [°]
omega_omegaBP = generator.omega_omegaBP;% speed controller crossover frequency [rad/s]
n = gearbox.ratio;

% I = 0.5;
% tau_c = 50e-6;
% L = 5e-3;
% R = 0.5;
% B = 0.191;
% p = 4;
% Lambda = 0.205;
% iq_omegaBP = 2*500*pi;

Yiq = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2); % generataor TF
Gc = 1/(1 + s*tau_c);                                     % power converter TF
Giq_noR = Yiq*Gc;                                         % G2

Riq = ki*(1 + tau_iq*s)/s;  % controller TF
GH = Riq*Giq_noR;           % close loop TF
pol = poles(GH, s);       % find the poles of the TF
max_pol = pol(1);

GH = simplify(subs(GH, tau_iq, -1/max_pol));

GH = norm(subs(GH, s, 1j*iq_omegaBP));
kiq = simplify(solve(GH == 1, ki));
eval(kiq)
kpq = kiq*(-1/max_pol);
eval(kpq)
