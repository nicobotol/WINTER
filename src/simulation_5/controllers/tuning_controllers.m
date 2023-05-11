%% Load PMSM transfer functions and design the controller

% Controller for the PMSM generator - q-axis
[Yiq, Gc, Riq, GR, G_cl, generator.kp, generator.ki, generator.kd, ... 
  generator.tau_d1] = PMSM_TF_pid(generator.design, generator.bode_plot);

% Controller for the PMSM generator - d-axis
[~, ~, ~, ~, ~, generator.kp_daxis, generator.ki_daxis, ...
  generator.kd_daxis, generator.tau_d1_daxis] = ...
  PMSM_daxis_TF_pid(generator.design, generator.bode_plot);

% Controller for the generator rotational speed
% [R_TG, GR] = TG_pi(generator.design_TG, generator.bode_plot_TG);
