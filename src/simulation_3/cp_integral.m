function [cp] = cp_integral(cp_partial, B, r_vector, lambda, A)
% Function used to perform the integration of the cp computed for each
% blade position

sum = 0;
% perform integration by trapezoidal method
for i=1:size(cp_partial)-1
  sum = sum + (cp_partial(i) + cp_partial(i+1)) * (r_vector(i + 1) - r_vector(i)) / 2;
end

cp = lambda * B / A * sum;

end