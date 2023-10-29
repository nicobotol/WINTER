function x_new = my_dynamic(x_est, T_R, param)
  
  states_len = size(x_est, 1);
  I_eq = param{1}; % Inertia
  B_eq = param{2}; % Damping
  K_opt = param{3}; % Power coefficient
  x_new = zeros(states_len,1);

  x_new(1) = x_est(1) + 1/I_eq*(T_R - K_opt*x_est(1)^2 - B_eq*x_est(1));
end