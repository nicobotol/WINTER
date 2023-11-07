  %   __  __                                                    _   
  %  |  \/  | ___  __ _ ___ _   _ _ __ ___ _ __ ___   ___ _ __ | |_ 
  %  | |\/| |/ _ \/ _` / __| | | | '__/ _ \ '_ ` _ \ / _ \ '_ \| __|
  %  | |  | |  __/ (_| \__ \ |_| | | |  __/ | | | | |  __/ | | | |_ 
  %  |_|  |_|\___|\__,_|___/\__,_|_|  \___|_| |_| |_|\___|_| |_|\__|
                                       
  function [z, z_hat] = measurement(x, R)
    % z_hat -> fake measuerement to be used in the KF
    states_len = size(x, 1);
    z = zeros(states_len, 1);
    z_hat = zeros(states_len, 1);

    z_hat(1) = x(1); % rotational speed [rad/s]
    z_hat(2) = x(2); % rotational speed [rad/s]
    z(1) = x(1); % rotational speed [rad/s]
    z(2) = x(2); % current
  end