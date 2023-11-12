clear; close; clc;
syms PR PGE PG Rs iq A cp rho V0 R lambda Lambda Beq n omega_R p pi

PR = 1/2*cp*rho*V0^3*pi*R^2; % rotor power
PG = PR - Beq*omega_R;
iq = 2/3*PG*n/(Lambda*p*omega_R);
PGE = PG - 3/2*Rs*iq^2;