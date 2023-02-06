function [a, a_prime, ct, cn, phi] = induction_factor(a, a_prime, R, r, ...
  lambda, beta, Theta_p, sigma, B, aoa_mat, cl_mat, cd_mat, thick_vector, ...
  thick)
%   Function to calculate a and a_prime
% a -> (INPUT) initial guess of axial induction factor
% a -> (OUTPUT)  axial induction factor after iterations
% a_prime -> (INPUT) intial guess of tangential induction factor
% a_prime -> (OUTPUT) tangential induction factor after iterations
% phi -> flow angle (rad)
% R -> rotor radius (m)
% r -> radius of the considered section (m)
% lambda -> coefficient lambda = (omega*R)/V0
% beta -> twist angle (rad)
% Theta_p -> pith angle (rad)
% cd -> cd factor
% ct -> ct factor
% cd -> cd factor
% cn -> cn factor
% c -> chord length of the actaul section (m)
% B -> # of blades (#)
% sigma -> solidity correction factor

  phi = atan(((1 - a)*R) / ((1 + a_prime) * lambda * r)); % flow angle

  Theta = beta + Theta_p; % Total pitch angle
  alpha = phi - Theta;  % angle of attack

  % function to interpolate cl and cd for the given value of alpha and t/c
  [cl, cd] = interpolate_cl_cd(aoa_mat, cl_mat, cd_mat, thick_vector, ...
    alpha, thick);

  cn = cl*cos(phi) + cd*sin(phi); 
  ct = cl*sin(phi) - cd*cos(phi); 
  
  F = tip_loss(B, R, r, phi); % Prandtl's tip loss correction 

  %cT = (1 - a)^2*cn*sigma / (sin(phi)^2);

  % better guess for a and a_prime
  if a <= 1/3 
    %cT = 4*a*F*(1 - a); 
    a = (sigma * cn * (1 - a)) / (4 * F * sin(phi)^2);
  else
    %cT = 4*a*F*(1 - 1/4*(5 - 3*a)*a);
    cT = (1 - a)^2*cn*sigma / (sin(phi)^2);
    beta_correction = 0.1;
    a_star = cT / (4*F*(1 - 1/4*(5 - 3*a)*a));
    a = beta_correction*a_star + (1 - beta_correction)*a;
  end
  
  a_prime = sigma*ct*(1 + a_prime) / (4*F*sin(phi)*cos(phi));

end