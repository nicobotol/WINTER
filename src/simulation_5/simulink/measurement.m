  %   __  __                                                    _   
  %  |  \/  | ___  __ _ ___ _   _ _ __ ___ _ __ ___   ___ _ __ | |_ 
  %  | |\/| |/ _ \/ _` / __| | | | '__/ _ \ '_ ` _ \ / _ \ '_ \| __|
  %  | |  | |  __/ (_| \__ \ |_| | | |  __/ | | | | |  __/ | | | |_ 
  %  |_|  |_|\___|\__,_|___/\__,_|_|  \___|_| |_| |_|\___|_| |_|\__|
                                       
  function [z, z_hat] = measurement(x, param, measure_len)
    % z_hat -> fake measuerement to be used in the KF

    z = zeros(measure_len, 1);
    z_hat = zeros(measure_len, 1);

    B_eq = param{6};
    rho_cp_R2 = param{2};
    
    if measure_len == 2 
      % In case of z_hat -> x = [omega; V0]
      z_hat(1) = x(1); % rotational speed [rad/s]
      z_hat(2) = (1/2*pi*rho_cp_R2*x(2)^3/x(1) - B_eq*x(1))*1e-6; % torque [MNm]
    
      % In case of z -> x=[omega; T_G]
      z(1) = x(1); % rotational speed [rad/s]
      z(2) = x(2); % torque [MNm]
    elseif measure_len == 3
      z_hat(1) = x(1); % rotational speed [rad/s]
      z_hat(2) = x(3)*1e-6; % current [mA]
      z_hat(3) = (1/2*pi*rho_cp_R2*x(2)^3/x(1) - B_eq*x(1))*1e-6; % torque

      z(1) = x(1); % rotational speed [rad/s]
      z(2) = x(2); % current [mA]
      z(3) = x(3); % torque [MNm]
    end


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
