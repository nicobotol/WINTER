function [rotor, generator, blade, T_R0] = run_Kopt_KoptGE(rho, ...
  lambda_vector, pitch_vector, lookup_cP, rotor, blade, ...
  generator, gearbox, V0_0, rated_values, lookup_Pitch, simulation, i)
% This test runs the coparison between the genrator control below rate with
% the K_opt based on the rotor and generator maximization
  
  if simulation.model == 3
  else 
    error('Wrong model selected, select simulation.mdl =3 in parameters.m')
  end

  V0_rated = rated_values(1);     % [m/s]
  omega_rated = rated_values(2);  % [rad/s]
  lambda_opt = rated_values(4);   % optimal TSR

  if i == 1 % case of K_opt based on generator
     [rotor, generator, blade, T_R0] = initial_conditions_GE(rho, ...
        lambda_vector, pitch_vector, lookup_cP, rotor, blade, ...
        generator, gearbox, V0_0, rated_values_P_GE, lookup_pitch_P_GE);
  elseif i == 2 % case of K_opt based on rotor
    generator.K_opt_GE =  generator.K_opt;
  end

end