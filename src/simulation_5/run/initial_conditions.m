function [rotor, generator, blade, T_R0, IMM] = initial_conditions(rho, lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, V0_0, rated_values, lookup_Pitch, IMM)
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

  % Initialization of the IMM
  P_est_initial = zeros(IMM.states_len, IMM.states_len, IMM.n_models);
  x_est_initial = zeros(IMM.states_len, IMM.n_models);  
  mu_initial = zeros(IMM.states_len, IMM.n_models)';  
  if IMM.states_len == 2
    IMM.x_est = [rotor.omega_R0; V0_0]; % state
  elseif IMM.states_len == 3
    IMM.x_est = [rotor.omega_R0; V0_0; 1e-3]; % state
  end
  
  for j = 1:IMM.n_models
    model{j}.x_est = IMM.x_est; % state
    model{j}.P_est = IMM.P_est; % covariance
    model{j}.mu = 1/IMM.n_models; % mode probability
    IMM.x_est_initial(1:IMM.states_len, j) = model{j}.x_est; 
    IMM.mu_initial(j, 1) = model{j}.mu;
    IMM.P_est_initial(1:IMM.states_len, 1:IMM.states_len, j) = model{j}.P_est; 
  end

end