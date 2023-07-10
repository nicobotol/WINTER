function [P_GE] = compute_P_GE(x, physical, lambda_vector, pitch_vector, lookup_cP, V0)

Rs = physical(1); % [Ohm] stator resistance
A = physical(2); % [m^2] rotor swept area
B = physical(3); % [kgm^2] transmission equiuvalent damping
p = physical(4); % [-] number of pole pairs
Lambda = physical(5); % [Wb] generator flux linkage
eta = physical(6); % [-] generator efficiency
R = physical(7); % [m]
rho = physical(8); % [kg/m^3] air density

lambda = x(1);
theta = x(2);

cp = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, theta,'cubic',0);
P_GE = -Rs*(A*cp*rho*V0^3 - 2*B*lambda^2*V0^2/R^2)^2*R^2/(9*lambda^2*V0^2*p^2*Lambda^2) + (-B*lambda*V0/R + A*cp*rho*V0^2*R/(2*lambda))*lambda*V0/R; % [W] generator power
P_GE = -P_GE; % in order to maximize the function minimize the inverse

end