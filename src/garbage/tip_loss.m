function [F] = tip_loss(B, R, r, phi)
% Function to compute the tip loss factor F
% F -> tip loss factor 
% B -> # of blades (#)
% R -> rotor diameter (m)
% r -> actual radius (m)
% phi -> flow angle (rad)

  F = 2/pi * acos( exp(- B * (R - r) / (2*r*sin(abs(phi)))));

end