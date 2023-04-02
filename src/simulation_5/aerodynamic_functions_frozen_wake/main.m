clear 
close all
clc

addpath("..\")

parameters

V0_vect = 11.44;
theta_offset = 10*pi/180;
delta_theta = 0.5*pi/180;
P = cell(1, length(V0_vect)); % cell storing the power for each wind speed for each pitch 
theta_vector = cell(1, length(V0_vect));

for v=1:length(V0_vect)
  V0 = V0_vect(v);
  lambda = omega_rated*rotor.R/V0; 
  theta = interp1(lookup_Pitch(1,:), lookup_Pitch(3,:), V0);

  % Determine the values of the induction factors along the blade
  [~, ~, ~, ~, a_store, a_prime_store] = ...
    cP_cT_partial(r_item_no_tip, r_vector, beta_vector, thick_vector, ...
    c_vector, rotor.blades, a_guess, a_prime_guess, rotor.R, lambda, theta, ...
    aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0, ...
    omega_rated, i_max);


  theta_vector{v} = theta-theta_offset:delta_theta:theta+theta_offset;
  cP = zeros(length(theta_vector), 1);
  for t=1:length(theta_vector)
    theta = theta_vector(t);
    [cp_partial, ~, ~, ~] = cP_cT_partial_FW( a_store, a_prime_store,...
      r_item_no_tip, r_vector, beta_vector, thick_vector, c_vector, rotor.R, ...
      lambda, theta, aoa_mat, cl_mat, cd_mat, thick_prof, rho, V0, omega_rated);
  
    cP(t) = lambda*rotor.blades/(rotor.A*rotor.R)*...
        trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial); 
  end
  P{v} = zeros(length(theta_vector), 1);
  P{v}(:) = 0.5*rotor.A*cP*V0^3*rho;

  dPdtheta = diff(P{v})/delta_theta

end