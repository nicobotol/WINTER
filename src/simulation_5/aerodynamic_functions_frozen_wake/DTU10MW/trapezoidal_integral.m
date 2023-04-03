function [out] = trapezoidal_integral(x_vect, y_vect)
% This function implements the integration with the trapezoidal rules

out = 0;
for i = 1:(size(x_vect,2) - 1)
  out = out + ((y_vect(i) + y_vect(i + 1))*(x_vect(i + 1) - x_vect(i))) / 2;
end
out = out + x_vect(1) * y_vect(1) / 2; % add the first part