function [cP] = cP_section(B, c, lambda, ct, a, R, phi)
% computation of the power factor for one specific section
% cP -> power factor
% a -> axial induction factor
% phi -> flow angle (rad)
% R -> rotor radius (m)
% lambda -> coefficient lambda = (omega*R)/V0
% ct -> ct factor
% c -> chord length of the actaul section (m)
% B -> # of blades (#)

cP = B*c*lambda*ct*(1 - a)^2 / (2*pi*R*sin(phi)^2);
end