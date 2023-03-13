clc
clear
close all

addpath("..\")
% load the parameters
parameters

omega_R = zeros(velocity_item, 1);  % [rad/s]
omega_G = zeros(velocity_item, 1);  % [rad/s]
T_R = zeros(velocity_item, 1);      % [Nm]
T_G = zeros(velocity_item, 1);      % [Nm]
P_R = zeros(velocity_item, 1);      % [W]
P_G = zeros(velocity_item, 1);      % [W]

for i=1:velocity_item
  V0 = velocity_vector(i);
  if V0 < V0_rated
    % rotor speed [rad/s]
    omega_R(i) = lambda_opt*V0/rotor.R;
    % generator rotational speed [rad/s]
    omega_G(i) = omega_R(i)/gearbox.ratio; 
    % rotor torque [Nm]
    T_R(i) = rotor.K_opt*omega_R(i)^2;
    % generator torque [Nm]
    T_G(i) = generator.K_opt*omega_G(i)^2;
    % generator power [W]
    P_G(i) = generator.K_opt*omega_G(i)^3;
    % rotor power [W]
    P_R(i) = rotor.K_opt*omega_R(i)^3; 
  else
    % rotor speed [rad/s]
    omega_R(i) = omega_rated;
    % generator rotational speed [rad/s]
    omega_G(i) = omega_R(i)/gearbox.ratio; 
    % generator torque [Nm]
    T_G(i) = rotor.P_rated/omega_rated*gearbox.ratio; % generator torque [Nm]
    % generator power [W]
    P_G(i) = generator.K_opt*omega_G(i)^3; 
  end
end

% show the value
figure()
subplot(3,1,1)
plot(velocity_vector, omega_R, 'o')
hold on
plot(velocity_vector, omega_G, '--')
ylabel('[rad/s]', 'Interpreter','latex')
hold off
subplot(3,1,2)
plot(velocity_vector, T_G)
ylabel('[Nm]', 'Interpreter','latex')
subplot(3,1,3)
plot(velocity_vector, P_G)
ylabel('[W]', 'Interpreter','latex')
xlabel('[m/s]', 'Interpreter', 'latex')

% Save the reusults in a lookup table
lookup_static_values = zeros(7, velocity_item);
lookup_static_values(1, :) = velocity_vector;
lookup_static_values(2, :) = omega_R;
lookup_static_values(3, :) = omega_G;
lookup_static_values(4, :) = T_R;
lookup_static_values(5, :) = T_G;
lookup_static_values(6, :) = P_R;
lookup_static_values(7, :) = P_G;

save('lookup_static_values.mat', 'lookup_static_values');