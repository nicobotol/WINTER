clear 
clc

parameters
P = 1e6*[0:0.01:10]; % mechanical power [W]
R = generator.Rs/10;
L = generator.Lq;
Lambda = generator.Lambda;
p = generator.p;
omega = omega_rated;

%%
iq = 2/3.*P*omega/(p*Lambda);
uq = 2/3.*P*omega/(p*Lambda)*R + omega*Lambda*p;

P_elettrica = uq.*iq;

eta = P./P_elettrica;
eta2 = 1./(4/9.*P*omega^2/(p*Lambda)^2*R+2/3*omega^2);
figure(1); clf;
plot(P, eta)
hold on
plot(P, eta2);
xlabel('Power [MW]')
ylabel('$\eta$')

%%
% uq = 3500;
% P = (uq - omega*Lambda*p)*1.5*p*Lambda/(omega*R);
% iq = 2/3*P;
