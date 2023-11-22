  %   __  __                                                    _   
  %  |  \/  | ___  __ _ ___ _   _ _ __ ___ _ __ ___   ___ _ __ | |_ 
  %  | |\/| |/ _ \/ _` / __| | | | '__/ _ \ '_ ` _ \ / _ \ '_ \| __|
  %  | |  | |  __/ (_| \__ \ |_| | | |  __/ | | | | |  __/ | | | |_ 
  %  |_|  |_|\___|\__,_|___/\__,_|_|  \___|_| |_| |_|\___|_| |_|\__|
                                       
  function [z, z_hat] = measurement_2measure(x, param, measure_len)
    % z_hat -> fake measuerement to be used in the KF
    % measure rotational speed and wind speed

    z = zeros(measure_len, 1);
    z_hat = zeros(measure_len, 1);

    % In case of z_hat -> x = [omega; V0]
    z_hat(1) = x(1); % rotational speed [rad/s]
    z_hat(2) = x(2); % measirement speed [m/s]

    % In case of z -> x=[omega; T_G]
    z(1) = x(1); % rotational speed [rad/s]
    z(2) = x(2); % wind speed [m/s]

    z = z + my_mvnrnd(zeros(measure_len, 1), param{4}, 1)';
      
  end

  function y = my_mvnrnd(mu, cov, n)
    % Cholesky decomposition
    L = chol(cov, 'lower');

    % Generate standard normal random numbers
    z = randn(n, length(mu));

    % Transform to multivariate normal
    y = L*z' + mu;

    y = y';
  end
