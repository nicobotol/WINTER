% This function wants to test the effectivness of the IMM for the estimation of the K_opt coefficient

close all; clear; clc;
parameters

omega = 0.6;
x_est = omega*ones(IMM.states_len, IMM.n_models);
P_est = zeros(IMM.states_len, IMM.states_len, IMM.n_models);
for i=1:IMM.n_models
  P_est(:,:,i) = 1e3*eye(IMM.states_len, IMM.states_len);
end
t_tot = 1000;
K_real = IMM.K_vector(2); % real gain
mu_vector = ones(IMM.n_models, 1)/IMM.n_models;
omega_tilde_store = [];
omega_store = [];
K_store = [];
mu_store =[];
K_opt = IMM.K_vector(2); % initial K_opt
omega_tilde = omega; % initial speed estimation

for t=1:t_tot
  T_G = K_opt*omega_tilde^2; % feedback torque
  T_R = T_G + 1e6*rand(1); % aero torque
  param{1} = I_eq;
  param{2} = B_eq;
  param{3} = K_real;
  omega = dynamic(omega, T_R, param); %+ my_mvnrnd(0, IMM.Q, 1); % real model evolution

  [omega_tilde, x_est, P_est, mu_vector, K_opt] = gain_IMM(omega, x_est, P_est, mu_vector, T_R, IMM, I_eq, B_eq);

  omega_store = [omega_store omega]; % real hystory
  omega_tilde_store = [omega_tilde_store omega_tilde];
  mu_store = [mu_store mu_vector];
  K_store = [K_store K_opt];
end

figure(); hold on; grid on;
plot(omega_store, 'LineWidth', 2)
plot(omega_tilde_store, '--', 'LineWidth', 2)
title('omega')
legend('Real', 'Est.')

figure(); hold on; grid on;
plot(K_store, 'LineWidth', 2)
yline(K_real, 'LineWidth', 2)

figure(); hold on; grid on;
for i=1:IMM.n_models
  plot(mu_store(i, :), 'DisplayName', ['Model ', num2str(i)]);
end
legend()

function [omega_tilde, x_est, P_est, mu_vector, K_opt] = gain_IMM(omega, x_est, P_est, mu_vector, T_R, IMM, I_eq, B_eq)
  % x_est -> (states_len, n_models) items containing the state estimation
  % P_est -> (states_len, states_len, n_models) containing the covariance
  
  n_models = IMM.n_models; % number of models to test
  states_len = IMM.states_len; % length of the states
  Pi = IMM.Pi; 
  R = IMM.R;
  Q = IMM.Q;
  K_vector = IMM.K_vector;
  mu = zeros(n_models);
  % Vector of parameters
  param = cell(6,1);
  param{1} = I_eq; % inertia
  param{2} = B_eq; % damping
  param{3} = 0; % K_opt
  param{4} = R; % measure noise
  param{5} = Q; % process noise 
  x_real = [omega];
  omega_tilde = 0;
  
  % Assamble the model cell
  model = cell(n_models, 1);
  for j=1:n_models
    param{3} = IMM.K_vector(j);
    model{j} = struct('K_opt', IMM.K_vector(j), 'x_est', x_est(1:states_len, j), 'P_est', P_est(1:states_len,1:states_len,j), 'z_hat', 0, 'P_prior', P_est(1:states_len,1:states_len,j), 'Lambda', 0, 'mu', mu_vector(j), 'mu_new', 0, 'x_tilde', zeros(states_len, 1), 'P_tilde', zeros(states_len));
  end
  
  % Simulate the measurement
  z = measurement(x_real, param{5}) + my_mvnrnd([0], R, 1)'; % rotational speed measurement [rad/s]
  
  % Filtering
  model = filtering(model, z, T_R, param);
  
  % Mode probability updating
  model = probability_updating(model, z, param);
  
  % State combination
  x_est_g = zeros(states_len, 1); % global state
  [x_est_g(1:states_len, 1), ~] = state_combination(model);
  omega_tilde = x_est_g(1); % estimation of the rotational speed

  [model, model_prob] = filter_interaction(model, Pi);
  [~, idx] = max(model_prob); % model with the highest probability
  K_opt = IMM.K_vector(idx);
  
  % Assemble the model cell
  for j=1:n_models
    x_est(1:states_len, j) = model{j}.x_est(1:states_len);
    P_est(1:states_len,1:states_len,j) = model{j}.P_est(1:states_len, 1:states_len);
    mu_vector(j) = model{j}.mu_new(1);
  end
  
  end
  
  %   _____ _  _______ 
  %  | ____| |/ /  ___|
  %  |  _| | ' /| |_   
  %  | |___| . \|  _|  
  %  |_____|_|\_\_|    
                     
  function [model] = EKF(model, z, T_G, param, matrices)
    % y -> measurement of the rotational speed [rad/s]
    
    x_est = model.x_est;
    P_est = model.P_est;
    states_len = size(x_est, 2);
  
    F = matrices.F;
    H = matrices.H;
    Q = matrices.Q;
    R = matrices.R;
    G = matrices.G;
  
    % Prediction
    x_est = dynamic(x_est(1:states_len), T_G, param);
    P_est = F*P_est*F' + G*Q*G';
    P_prior = P_est;
    states_len = size(x_est, 1); % number of states
  
    % Update
    [~, z_hat] = measurement(x_est, param{4}); % estimation of the measurement using the sensor's model
    Innovation = z - z_hat;
    S_Inno = H*P_est*H' + R;
    W = P_est*H'*inv(S_Inno); % kalman gain
    x_est = x_est + W*Innovation; % update state estimate
    P_est = (eye(states_len) - W*H)*P_est; % update covariance matrix
  
  
    model.x_est(1:states_len) = x_est(1:states_len);
    model.P_est(1:states_len, 1:states_len) = P_est(1:states_len, 1:states_len);
    model.z_hat = z_hat;
    model.P_prior = P_prior;
  end
  
  %   ____                              _      
  %  |  _ \ _   _ _ __   __ _ _ __ ___ (_) ___ 
  %  | | | | | | | '_ \ / _` | '_ ` _ \| |/ __|
  %  | |_| | |_| | | | | (_| | | | | | | | (__ 
  %  |____/ \__, |_| |_|\__,_|_| |_| |_|_|\___|
  %         |___/                              
  
  function x_new = dynamic(x_est, T_R, param)
  
  states_len = size(x_est, 1);
    I_eq = param{1}; % Inertia
    B_eq = param{2}; % Damping
    K_opt = param{3}; % Power coefficient
    x_new = zeros(states_len,1);
  
    x_new(1) = x_est(1) + 1/I_eq*(T_R - K_opt*x_est(1)^2 - B_eq*x_est(1));
  end
  
  %   __  __                                                    _   
  %  |  \/  | ___  __ _ ___ _   _ _ __ ___ _ __ ___   ___ _ __ | |_ 
  %  | |\/| |/ _ \/ _` / __| | | | '__/ _ \ '_ ` _ \ / _ \ '_ \| __|
  %  | |  | |  __/ (_| \__ \ |_| | | |  __/ | | | | |  __/ | | | |_ 
  %  |_|  |_|\___|\__,_|___/\__,_|_|  \___|_| |_| |_|\___|_| |_|\__|
                                       
  function [z, z_hat] = measurement(x, R)
    % z_hat -> fake measuerement to be used in the KF
  
    z_hat = x(1); % rotational speed [rad/s]
    z = x(1); %+ mvnrnd([0], R)'; % rotational speed [rad/s]
  end
  
  %    ____                            _                         _     
  %   / ___|___  _ __ ___  _ __  _   _| |_ ___   _ __ ___   __ _| |_   
  %  | |   / _ \| '_ ` _ \| '_ \| | | | __/ _ \ | '_ ` _ \ / _` | __|  
  %  | |__| (_) | | | | | | |_) | |_| | ||  __/ | | | | | | (_| | |_ _ 
  %   \____\___/|_| |_| |_| .__/ \__,_|\__\___| |_| |_| |_|\__,_|\__(_)
  %                       |_|                                          
  
  function matrices = compute_matrices(x_est, param)
  
    I_eq = param{1}; % Inertia
    B_eq = param{2}; % Density
    K_opt = param{3}; % Power coefficient
  
    F = [-K_opt/I_eq*x_est(1)-B_eq/I_eq];
    H = [1];
    G = [1];
  %   matrices = cell(4, 1);
    matrices = struct('F', F, 'H', H, 'G', G);
  end
  
  %   _____ _ _ _            _             
  %  |  ___(_) | |_ ___ _ __(_)_ __   __ _ 
  %  | |_  | | | __/ _ \ '__| | '_ \ / _` |
  %  |  _| | | | ||  __/ |  | | | | | (_| |
  %  |_|   |_|_|\__\___|_|  |_|_| |_|\__, |
  %                                  |___/ 
  
  function model = filtering(model, z, T_G, param)
    n_models = size(model, 1);
    states_len = size(model{1}.x_est, 2);
    for j=1:n_models
      param{4} = model{j}.K_opt;
      % compute the linearized matrices
      matrices = compute_matrices(model{j}.x_est(1:states_len), param);
      matrices.R = param{4};
      matrices.Q = param{5};
  
      [model{j}] = EKF(model{j}, z, T_G, param, matrices);
    end
  end
  
  %   ____            _                       _       _   _             
  %  |  _ \ _ __ ___ | |__    _   _ _ __   __| | __ _| |_(_)_ __   __ _ 
  %  | |_) | '__/ _ \| '_ \  | | | | '_ \ / _` |/ _` | __| | '_ \ / _` |
  %  |  __/| | | (_) | |_) | | |_| | |_) | (_| | (_| | |_| | | | | (_| |
  %  |_|   |_|  \___/|_.__/   \__,_| .__/ \__,_|\__,_|\__|_|_| |_|\__, |
  %                                |_|                            |___/ 
  
  function model = probability_updating(model, z, param)
    n_models = size(model, 1);
    for j=1:n_models
      % compute the linearized matrices
      matrices = compute_matrices(model{j}.x_est, param);
      matrices.R = param{4};
      matrices.Q = param{5}; 
      H = matrices.H;
      R = matrices.R;
      G = matrices.G;
      z_hat = model{j}.z_hat;
      P_prior = model{j}.P_prior;
  
      S = H*P_prior*H' + G*R*G';
      zeta = z - z_hat;
      model{j}.Lambda = 1/sqrt(2*pi*det(S))*exp(-1/2*zeta'*inv(S)*zeta);
    end
    % Compute the denominator
    den = 0;
    for j=1:n_models
      den = den + model{j}.mu*model{j}.Lambda;
    end
    % Update model probability
    for j=1:n_models
      model{j}.mu = model{j}.mu*model{j}.Lambda/den;
    end
  end
  
  %   ____  _        _                             _      
  %  / ___|| |_ __ _| |_ ___    ___ ___  _ __ ___ | |__   
  %  \___ \| __/ _` | __/ _ \  / __/ _ \| '_ ` _ \| '_ \  
  %   ___) | || (_| | ||  __/ | (_| (_) | | | | | | |_) | 
  %  |____/ \__\__,_|\__\___|  \___\___/|_| |_| |_|_.__(_)
                                                        
  function [x_est, P_est] = state_combination(model)
    n_models = size(model, 1);
    states_len = size(model{1}.x_est, 1);
    x_est = zeros(states_len,1);
    P_est = zeros(states_len,states_len);
    for j=1:n_models
      x_est(1:states_len) = x_est(1:states_len) + model{j}.mu*model{j}.x_est(1:states_len);
    end
    
    for j=1:n_models
      P_est = P_est + model{j}.mu*(model{j}.P_est + (x_est - model{j}.x_est)*(x_est - model{j}.x_est)');
    end
  
  end
  
  %   _____ _ _ _              _       _                      _   _             
  %  |  ___(_) | |_ ___ _ __  (_)_ __ | |_ ___ _ __ __ _  ___| |_(_) ___  _ __  
  %  | |_  | | | __/ _ \ '__| | | '_ \| __/ _ \ '__/ _` |/ __| __| |/ _ \| '_ \ 
  %  |  _| | | | ||  __/ |    | | | | | ||  __/ | | (_| | (__| |_| | (_) | | | |
  %  |_|   |_|_|\__\___|_|    |_|_| |_|\__\___|_|  \__,_|\___|\__|_|\___/|_| |_|
                                                                              
  function [model, model_prob] = filter_interaction(model, Pi)
    n_models = size(model, 1);
    states_len = size(model{1}.x_est, 1);
    mu = zeros(n_models, n_models);
  
    for j=1:n_models
      model{j}.mu_new = 0;
      for i=1:n_models
          model{j}.mu_new = model{j}.mu_new + Pi(i, j)*model{i}.mu; 
      end
    end
  
    for j=1:n_models
      for i=1:n_models
          mu(i, j) = Pi(i, j)*model{i}.mu/model{j}.mu_new;
      end
    end
  
    for j=1:n_models
      model{j}.x_tilde(1:states_len, 1) = zeros(states_len, 1);
      model{j}.P_tilde = zeros(states_len, states_len);
      for i=1:n_models
        model{j}.x_tilde = model{j}.x_tilde + mu(i, j)*model{i}.x_est;
      end
      for i=1:n_models
          model{j}.P_tilde = model{j}.P_tilde + mu(i, j)*(model{i}.P_est + (model{j}.x_tilde - model{i}.x_est)*(model{j}.x_tilde - model{i}.x_est)');
      end
      model{j}.x_est = model{j}.x_tilde;
      model{j}.P_est = model{j}.P_tilde;
    end
  
    % Store the value of the probability
    model_prob = zeros(n_models, 1);
    for j=1:n_models
      model_prob(j) = model{j}.mu_new;
    end
  end
  
  %    ___        _               _               _     
  %   / _ \ _   _| |_ _ __  _   _| |_    ___  ___| |_   
  %  | | | | | | | __| '_ \| | | | __|  / _ \/ __| __|  
  %  | |_| | |_| | |_| |_) | |_| | |_  |  __/\__ \ |_ _ 
  %   \___/ \__,_|\__| .__/ \__,_|\__|  \___||___/\__(_)
  %                  |_|                                
  
  function model = output_estimation(model, param)
    n_models = size(model, 1);
    R = param{3}; % Radius
    lambda = 0;
    for j=1:n_models
      x_est_1 = model{j}.x_est(1);
      x_est_2 = model{j}.x_est(2);
      lambda = R*x_est_1/x_est_2; 
    end
  end
  
  %                                       _ 
  %   _ __ _____   ___ __  _ __ _ __   __| |
  %  | '_ ` _ \ \ / / '_ \| '__| '_ \ / _` |
  %  | | | | | \ V /| | | | |  | | | | (_| |
  %  |_| |_| |_|\_/ |_| |_|_|  |_| |_|\__,_|
                                          
  function y = my_mvnrnd(mu, cov, n)
      % Cholesky decomposition
      L = chol(cov, 'lower');
  
      % Generate standard normal random numbers
      z = randn(n, length(mu));
  
      % Transform to multivariate normal
      y = L*z' + mu;
  
      y = y';
  end
  