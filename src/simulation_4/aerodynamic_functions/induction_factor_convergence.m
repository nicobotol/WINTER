function [a, a_prime, ct, cn, phi] = induction_factor_convergence(a_guess, ...
  a_prime_guess, R, r, lambda, beta, Theta_p, B, sigma, aoa_mat, cl_mat, ...
  cd_mat, thick_prof, thick, fake_zero, i_max)
% This function performs the calculations to make the induction functions
% converge
% a -> (INPUT) initial guess of axial induction factor
% a -> (OUTPUT)  axial induction factor after iterations
% a_prime -> (INPUT) intial guess of tangential induction factor
% a_prime -> (OUTPUT) tangential induction factor after iterations
% R -> rotor radius (m)
% r -> radius of the considered section (m)
% lambda -> coefficient lambda = (omega*R)/V0
% beta -> twist angle (rad)
% Theta_p -> pith angle (rad)
% cd -> cd factor
% cl -> cl factor
% ct -> ct factor
% cn -> cn factor
% c -> chord length of the actaul section (m)
% B -> # of blades (#)
% fake_zero -> numerical approximation of 0 (provided in parameters.m)
% i_max ->  max number of iterations for the convergence cycle (provided in
%           parameters.m)

a = a_guess;
a_prime = a_prime_guess;

for i = 1:i_max
  % initialize a and a_prime
  a_old = a;
  a_prime_old = a_prime;
  
  % update a and a_prime
  [a, a_prime, ct, cn, phi] = induction_factor(a, a_prime, R, r, lambda, ...
    beta, Theta_p, sigma, B, aoa_mat, cl_mat, cd_mat, thick_prof, thick);
  
  % compute the error
  epsilon = abs(a_old - a);
  epsilon_prime = abs(a_prime_old - a_prime);
  
  % brake the code if convergence
  if (epsilon < fake_zero) && (epsilon_prime < fake_zero)
    break
  end

end

end