function prova_sigma
parameters

% omega = 0:0.01:omega_rated;
% 
% sigma_v = sigma(generator.omega1_max, generator.omega2_max, omega);
% 
% figure()
% plot(omega, sigma_v)

omega = 0:0.01:4;

sigma_v = sigma(1,2, omega);

figure()
plot(omega, sigma_v)
% plot(omega, -2*omega.^3+9*omega.^2-12*omega+5)
end

function y = sigma(x0, x1, x)
  a = [2 -3*(x0 + x1) 6*x1*x0 (x0 - 3*x1)*x0^2]/(x0 - x1)^3;  % poly coeff.
  y = zeros(length(x), 1);
  for i = 1:length(x)
    if x(i) < x0
      y(i) = 0;
    elseif x0 <= x(i) && x(i) <= x1
      y(i) = polyval(a, x(i));  
    elseif x(i) > x1
      y(i) = 1;
    end
  end

end

    