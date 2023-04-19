function [rotor, generator, blade] = initial_conditions(rotor, blade, ...
  generator, gearbox, V0_0, rated_values, lookup_Pitch)
% This functions sets the initial conditions fo the simulation
% V0_0 -> initial windspeed [m/s]

V0_rated = rated_values(1);     % [m/s]
omega_rated = rated_values(2);  % [rad/s]
lambda_opt = rated_values(4);   % optimal TSR

% initial rotational speed [rad/s]
if V0_0 < V0_rated
  rotor.omega_R0 = lambda_opt*V0_0/rotor.R - 0.5*pi/30; 
else
  rotor.omega_R0 = omega_rated + 0.5*pi/30;
end

% initial generator torque (generator side)
generator.T_G0 = generator.K_opt*gearbox.ratio*rotor.omega_R0^2; % [Nm]
generator.P_G0 = generator.K_opt*rotor.omega_R0^3; % [W]

% initial blade pitch angle
blade.pitch0 = interp1(lookup_Pitch(1,:), lookup_Pitch(3,:), V0_0);

end