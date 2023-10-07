function [rotor, generator, blade, T_R0] = initial_conditions_GE(rho, lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, V0_0, V0_rated, rated_values_P_GE, lookup_pitch_P_GE)
% This functions sets the initial conditions fo the simulation
% V0_0 -> initial windspeed [m/s]

omega_rated_GE = rated_values_P_GE(3);  % [rad/s]
lambda_GE = rated_values_P_GE(1);   % optimal TSR

% initial rotational speed [rad/s]
if V0_0 < V0_rated
  rotor.omega_R0 = lambda_GE*V0_0/rotor.R; 
else
  rotor.omega_R0 = omega_rated_GE;
end

% initial generator torque (generator side)
generator.T_G0 = generator.K_opt_GE*gearbox.ratio*rotor.omega_R0^2; % [Nm]
generator.P_G0 = generator.K_opt_GE*rotor.omega_R0^3; % [W]

% initial blade pitch angle
blade.pitch0 = interp1(lookup_pitch_P_GE(1,:), lookup_pitch_P_GE(2,:), V0_0);

% Intial aero torque
lambda0 = rotor.omega_R0*rotor.R/V0_0;
cP0 = interp2(lambda_vector, pitch_vector, lookup_cP, lambda0, blade.pitch0);

P_R0 = 0.5*rho*V0_0^3*pi*rotor.R^2*cP0; % rotor power [W]
T_R0 = P_R0/rotor.omega_R0;            % rotor torque [Nm]

end