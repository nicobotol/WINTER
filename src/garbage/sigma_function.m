function [sigma] = sigma_function(c, B, r)
% function to compute the solidty correction factor
% sigma -> solidity correction factor
% c -> chord length (m)
% B -> number of blades (#)
% r -> actual radius (m)

  sigma = (c*B) / (2*pi*r);

end