%% Load PMSM transfer functions and design the controller

% Controller for the PMSM generator
[Yiq, Gc, Riq, GR, G_cl, generator.kp, generator.ki, generator.kd, ... 
  generator.tau_d1] = PMSM_TF_pid(generator.design, generator.bode_plot);

% Controller for the generator rotational speed
% [R_TG, GR] = TG_pi(generator.design_TG, generator.bode_plot_TG);
