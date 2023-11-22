  %   __  __                                                    _   
  %  |  \/  | ___  __ _ ___ _   _ _ __ ___ _ __ ___   ___ _ __ | |_ 
  %  | |\/| |/ _ \/ _` / __| | | | '__/ _ \ '_ ` _ \ / _ \ '_ \| __|
  %  | |  | |  __/ (_| \__ \ |_| | | |  __/ | | | | |  __/ | | | |_ 
  %  |_|  |_|\___|\__,_|___/\__,_|_|  \___|_| |_| |_|\___|_| |_|\__|
                                       
  function [z, z_hat] = measurement_3measure(x, param, measure_len)
    % z_hat -> fake measuerement to be used in the KF
    % measure rotational speed, torque, and wind speed

    z = zeros(measure_len, 1);
    z_hat = zeros(measure_len, 1);
    rho = param{2};
    R = param{8};
    B_eq = param{6};
    cp = param{9};
    rho_cp_R2 = rho*cp*R^2;
    
    % In case of z_hat -> x = [omega; V0]
    z_hat(1) = x(1); % rotational speed [rad/s]
    z_hat(2) = (1/2*pi*rho_cp_R2*x(2)^3/x(1) - B_eq*x(1))*1e-6;
    z_hat(3) = x(2); % wind speed [m/s]

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
