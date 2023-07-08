clc
clear
close all

parameters;
physical(1) = generator.Rs; % [Ohm] stator resistance
physical(2) = rotor.A; % [m^2] rotor swept area
physical(3) = B_eq; % [kgm^2] transmission equiuvalent damping
physical(4) = generator.p; % [-] number of pole pairs
physical(5) = generator.Lambda; % [Wb] generator flux linkage
physical(6) = generator.eta; % [-] generator efficiency
physical(7) = rotor.R; % [m]
physical(8) = rho; % [kg/m^3]

V0_b = 4:0.1:V0_rated+0.001; % [m/s] wind speed 

% Maximize the power output by selecting the proper combination of TSR and pitch angle
% The variable x containts x = (TSR, pitch angle)
lb = [5, -10*pi/180]; % variable lower bound
ub = [10, 40*pi/180]; % variable upper bound
P = zeros(length(V0_b),1);
min_v = zeros(length(V0_b),2);
x0 = [5, 0];
for i=1:length(V0_b)
  V0 = V0_b(i);    
  [min_v(i, :), P(i, :)] = fmincon(@(x)compute_P_GE(x, physical, lambda_vector, pitch_vector, lookup_cP, V0), x0, [], [], [], [], lb, ub);
end

%% Compute the cp corresponding to the optimized parameters
cp = interp2(lambda_vector, pitch_vector, lookup_cP, min_v(:, 1), min_v(:, 2));

figure(); 
plot(V0_b, min_v(:, 1))
title('lambda')

figure(); 
plot(V0_b, min_v(:, 2)*180/pi)
title('theta')

figure();
plot(V0_b, cp)
yline(cp_max, 'r--')
title('$c_P$')
