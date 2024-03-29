function [rotor, generator, blade, T_R0, omega_rated_GE, IMM] = run_Kopt_KoptGE(rho, lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, V0_0, rated_values, rated_values_P_GE, rated_values_P_GE_no_B, lookup_Pitch, lookup_pitch_P_GE, simulation,IMM, B_eq, i)
  % This test runs the coparison between the genrator control below rate with
  % the K_opt based on the rotor and generator maximization
  
  if simulation.model == 3
  else 
    error('Wrong model selected, select simulation.mdl =3 in parameters.m')
  end
  V0_rated =  rated_values(1);

  if rem(i, 2) == 1 % case of K_opt based on generator
    if B_eq > 1 % case with damping
    lambda_GE = rated_values_P_GE(1);
    cp_GE = rated_values_P_GE(2);
    generator.K_opt_GE = rho*pi*rotor.R^5*cp_GE/(2*lambda_GE^3);
    omega_rated_GE = rated_values_P_GE(3); % restore the rated speed
    else % case without damping
      lambda_GE = rated_values_P_GE_no_B(1);
      cp_GE = rated_values_P_GE_no_B(2);
      generator.K_opt_GE = rho*pi*rotor.R^5*cp_GE/(2*lambda_GE^3);
      omega_rated_GE = rated_values_P_GE_no_B(3); % restore the rated speed
    end
    [rotor, generator, blade, T_R0, IMM] = initial_conditions_GE(rho, lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, V0_0, V0_rated, rated_values_P_GE, lookup_pitch_P_GE, IMM);
  else % case of K_opt based on rotor
    generator.K_opt_GE = generator.K_opt;
    omega_rated_GE = rated_values(2); % restore the rated speed
    [rotor, generator, blade, T_R0, IMM] = initial_conditions(rho, lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, V0_0, rated_values, lookup_Pitch, IMM);
  end

  % if rem(i, 2) == 1
  %   K_gain = generator.K_opt_GE;
  % else
  %   K_gain = generator.K_opt;
  % end
  % [rotor, generator, blade, T_R0] = initial_conditions(rho, lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, V0_0, rated_values, lookup_Pitch);

end