function [z] = my_measure(x)
  % z_hat -> fake measuerement to be used in the KF

  z_hat = x(1); % rotational speed [rad/s]
  z = x(1); %+ mvnrnd([0], R)'; % rotational speed [rad/s]
end