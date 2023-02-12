close all
clc
s = tf('s');

G1 = 1/(1 + 100*s);
G2 = 1/(1 + 2000*s);

figure(1)
bodeplot(G1, G2, G2*(1+2000*s)/s)
legend('100', '2000', 'ctrl', 'location', 'best')
grid on

%%
G3 = 1/((1 + 100*s^2));
G4 = G3*(1 + s*10);
pole(G3)
pole(G4)
figure(2)
bodeplot(G3, G4)
legend('NC', 'C')