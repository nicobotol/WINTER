function [generator, IMM] = run_K_opt_sensitivity(generator, IMM, i)
  % K_opt_unchanged -> value of K_opt not rescaled

  % Rescale the gain times a constant
  if i == 1 % run the IMM
    IMM.enable = 1;
  else % run the constant   
    IMM.enable = 0;
    generator.K_opt_GE = IMM.K_vector(i - 1);
  end

end