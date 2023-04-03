function [ct, cn, phi] = induction_factor_FW(a, a_prime, R, r, ...
  lambda, beta, Theta_p, aoa_mat, cl_mat, cd_mat, thick_prof, ...
  thick)
%   Function to calculate ct and cn in frozen wake
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
  [cl, cd] = interpolate_cl_cd(aoa_mat, cl_mat, cd_mat, thick_prof, ...
    alpha, thick);

  cn = cl*cos(phi) + cd*sin(phi); 
  ct = cl*sin(phi) - cd*cos(phi); 
  
end