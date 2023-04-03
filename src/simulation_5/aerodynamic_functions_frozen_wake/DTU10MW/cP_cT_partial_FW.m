function [cp_partial, cT_partial, pt, pn] = cP_cT_partial_FW( a_store, a_prime_store,...
  r_item_no_tip, r_vector, beta_vector, thick_vector, c_vector, R, ...
  lambda, Theta_p, aoa_mat, cl_mat, cd_mat, thick_prof, rho, V0, omega)

  % Function to compute the cp and cT of the local blade section assuming frozen wake, looping
  % along all blade sections
  
  cp_partial = zeros(1, r_item_no_tip);
  cT_partial = zeros(1, r_item_no_tip);
  pt = zeros(1, r_item_no_tip);
  pn = zeros(1, r_item_no_tip);

  for i=1:r_item_no_tip % loop over the blade positions 
    r = r_vector(i);
    beta = beta_vector(i);
    thick = thick_vector(i);
    c = c_vector(i);

    a = a_store(i);
    a_prime = a_prime_store(i);

    % compute ct  and cn using the frozen a and a_prime
    [ct, cn, ~] = induction_factor_FW(a, a_prime, R, r, ...
    lambda, beta, Theta_p, aoa_mat, cl_mat, cd_mat, thick_prof, ...
    thick);
    
    cp_partial(i) = r*((1 - a)^2 + (lambda*r/R*(1 + a_prime))^2)*c*ct;
    cT_partial(i) = ((1 - a)^2 + (lambda*r/R*(1 + a_prime))^2)*c*cn;

    pt(i) = 0.5 * rho * (V0^2*(1 - a)^2 + (1 + a_prime)^2*omega^2*r^2)*c*ct; % pt coefficient
    pn(i) = 0.5 * rho * (V0^2*(1 - a)^2 + (1 + a_prime)^2*omega^2*r^2)*c*cn; % pn coefficient
  end % stop looping over blade position r
end