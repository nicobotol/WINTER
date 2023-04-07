%% Load PMSM transfer functions and design the controller

% Controller for the PMSM generator
[Yiq, Gc, Riq, GR, G_cl] = PMSM_TF_pid(generator.design, generator.bode_plot);

% Controller for the generator rotational speed
% [Romega, generator_ki, generator_kp, generator_kd] = ...
%   omega_G_pi(generator.design_omega, generator.bode_plot);
