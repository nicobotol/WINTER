function [generator] = run_Kopt_KoptGE(generator, simulation, i)
% This test runs the coparison between the genrator control below rate with
% the K_opt based on the rotor and generator maximization

  if simulation.model == 3
  else 
    error('Wrong model selected, select simulation.mdl =3 in parameters.m')
  end

  if i == 1 % case of K_opt based on generator
    % do nothing
  elseif i == 2 % case of K_opt based on rotor
    generator.K_opt_GE =  generator.K_opt;
  end

end