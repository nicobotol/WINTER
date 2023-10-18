function [generator] = run_K_opt_sensitivity(generator, i)
  % K_opt_unchanged -> value of K_opt not rescaled

  % Rescale the gain times a constant
  generator.gain_K_opt = generator.K_opt_sensitivity(i);

end