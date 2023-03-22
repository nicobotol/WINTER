%% Load PMSM transfer functions and design the controller

% Controller for the PMSM generator
[Yiq, Gc, Riq, GR, G_cl] = PMSM_TF_pid(generator.design, generator.bode_plot);

% Controller for the generator rotational speed
[Romega generator_ki, generator_kp, generator_kd, generator_tau_d1] = ...
  omega_G_pid(generator.design2, generator.bode_plot, G_cl);