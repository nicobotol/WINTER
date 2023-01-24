clear all
close all
clc

load("Pitching.mat");

P_rated = 10e6;  % rated power [W]
R = 178.3/2; % radius [m]
V0_rated = 11.44; % rated windspeed [m/s]
rho = 1.225; % air densityt

cp_max = 0.465;
lambda_opt = 7.857;
omega_rated = V0_rated*lambda_opt/R;
cP_c = [0.78, 151, 0.58, 0.002, 2.14, 13.2, 20.9, -0.002, -0.008];

V0_vector = 4:0.5:25;
V0_vector(end) = V0_rated;
V0_vector = sort(V0_vector);

P = zeros(length(V0_vector), 1);
omega = zeros(length(V0_vector), 1);
P_app = zeros(length(V0_vector), 1);

for i=1:length(V0_vector)
  V0 = V0_vector(i);
  theta_p = interp1(Theta_p_limit(1,:), Theta_p_limit(3,:), V0)*180/pi; % pitch

  % Standard method
  if V0 < V0_rated
    P(i) = 0.5*rho*pi*R^2*V0^3*cp_max;
    omega(i) = V0*lambda_opt/R; 
    lambda = lambda_opt;
  else
    P(i) = P_rated;
    omega(i) = omega_rated;
    lambda = omega(i)*R/V0;
  end

  % Approximation method

  lambda_w = 1/(1/(lambda + cP_c(8)*theta_p) - cP_c(9)/(theta_p^3 + 1));
  cP = cP_c(1)*(cP_c(2)/lambda_w - cP_c(3)*theta_p - ...
    cP_c(4)*theta_p^cP_c(5) - cP_c(6))*exp(-cP_c(7)/lambda_w);

  P_app(i) = 0.5*rho*V0^3*pi*R^2*cP; % rotor power [W]
end

figure(1)
plot(V0_vector, P);
hold on
plot(V0_vector, P_app);
legend('normal', 'app')
