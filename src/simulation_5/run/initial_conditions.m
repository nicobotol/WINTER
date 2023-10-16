function [rotor, generator, blade, T_R0, x_est_initial, P_est_initial] = initial_conditions(rho, lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, V0_0, rated_values, lookup_Pitch, IMM)
% This functions sets the initial conditions fo the simulation
% V0_0 -> initial windspeed [m/s]

V0_rated = rated_values(1);     % [m/s]
omega_rated = rated_values(2);  % [rad/s]
lambda_opt = rated_values(4);   % optimal TSR

% initial rotational speed [rad/s]
if V0_0 < V0_rated
  rotor.omega_R0 = lambda_opt*V0_0/rotor.R; 
else
  rotor.omega_R0 = omega_rated;
end

% initial generator torque (generator side)
generator.T_G0 = generator.K_opt*gearbox.ratio*rotor.omega_R0^2; % [Nm]
generator.P_G0 = generator.K_opt*rotor.omega_R0^3; % [W]

% initial blade pitch angle
blade.pitch0 = interp1(lookup_Pitch(1,:), lookup_Pitch(3,:), V0_0);

% Intial aero torque
lambda0 = rotor.omega_R0*rotor.R/V0_0;
cP0 = interp2(lambda_vector, pitch_vector, lookup_cP, lambda0, blade.pitch0);

P_R0 = 0.5*rho*V0_0^3*pi*rotor.R^2*cP0; % rotor power [W]
T_R0 = P_R0/rotor.omega_R0;            % rotor torque [Nm]

IMM.P_est = 1e3*eye(IMM.states_len, IMM.states_len); % initial covariance matrix
IMM.x_est = [rotor.omega_R0];       % initial state estimate
for j = 1:IMM.n_models
  x_est_initial(1:IMM.states_len, j) = IMM.x_est; 
  P_est_initial(1:IMM.states_len, 1:IMM.states_len, j) = IMM.P_est; 
end

end