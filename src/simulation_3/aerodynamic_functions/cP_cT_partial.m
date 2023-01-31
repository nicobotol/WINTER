function [cp_partial, cT_partial, pt, pn] = cP_cT_partial(r_item_no_tip, r_vector, ...
  beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda, ...
  Theta_p, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0, omega, i_max)

% Function to compute the cp and cT of the local blade section, looping
% along all blade sections

    for i=1:r_item_no_tip % loop over the blade positions 
      r = r_vector(i);
      beta = beta_vector(i);
      thick = thick_vector(i);
      c = c_vector(i);
      sigma = sigma_function(c, B, r);
      %cp_partial = zeros(1, r_item_no_tip);
      %cT_partial = zeros(1, r_item_no_tip);
      %pt = zeros(1, r_item_no_tip);
      %pn = zeros(1, r_item_no_tip);

      % compute the a and a_prime with the iterative method
      [a, a_prime, ct, cn, ~] = induction_factor_convergence(a_guess, ...
        a_prime_guess, R, r, lambda, beta, Theta_p, B, sigma, aoa_mat, ...
        cl_mat, cd_mat, thick_prof, thick, fake_zero, i_max);
      
      cp_partial(i) = r*((1 - a)^2 + (lambda*r/R*(1 + a_prime))^2)*c*ct;
      cT_partial(i) = ((1 - a)^2 + (lambda*r/R*(1 + a_prime))^2)*c*cn;

      pt(i) = 0.5 * rho * (V0^2*(1 - a)^2 + (1 + a_prime)^2*omega^2*r^2)*c*ct; % pt coefficient
      pn(i) = 0.5 * rho * (V0^2*(1 - a)^2 + (1 + a_prime)^2*omega^2*r^2)*c*cn; % pn coefficient

    end % stop looping over blade position r
end