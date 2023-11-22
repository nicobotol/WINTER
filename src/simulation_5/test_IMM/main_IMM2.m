% This function wants to test the effectivness of the IMM for the estimation of the K_opt coefficient. What is done is to replicate the simulink simulation 

close all; clear; clc;
parameters
[rotor, generator, blade, T_R0, IMM] = initial_conditions(rho,lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, 8, rated_values, lookup_Pitch, IMM);
V0 = 8;
omega = lambda_opt*V0/rotor.R;
% x_est = omega*ones(IMM.states_len, IMM.n_models);
% P_est = zeros(IMM.states_len, IMM.states_len, IMM.n_models);
x_est = IMM.x_est_initial;
P_est = IMM.P_est_initial;
for i=1:IMM.n_models
  P_est(:,:,i) = 1e3*eye(IMM.states_len, IMM.states_len);
end
t_tot = 1000;
K_real = IMM.K_vector(3); % real gain
mu_vector = ones(IMM.n_models, 1)/IMM.n_models;
omega_tilde_store = [];
omega_store = [];
K_store = [];
mu_store =[];
T_R_store = [];
T_G_store = [];
x_est_store = zeros(IMM.n_models, t_tot);
K_opt = mean(IMM.K_vector); % initial gain estimation
omega_tilde = omega; % initial speed estimation
v=8;
x_real = [0.5; 8];
for t=1:t_tot
  K_opt = generator.K_opt_GE;
  T_R = 0.5*V0^3*rho*cp_GE*rotor.A + 5e2*(rand(1) - 0.5)*2; % aero torque
  T_G = K_opt*omega^2; % feedback torque
  param = IMM.param;
  x_real = dynamic(x_real, T_R, param, dt); %+ my_mvnrnd(0, IMM.Q, 1); % real model evolution
  z = x_real(1) + my_mvnrnd(0, IMM.W, 1); % measurement 
  [x_est, P_est, mu_vector, omega_tilde, rho, K_opt] = gain_IMM(x_est, P_est, mu_vector, IMM, I_eq, dt, z, generator, B_eq);

  omega_store = [omega_store omega]; % real hystory
  omega_tilde_store = [omega_tilde_store omega_tilde];
  mu_store = [mu_store mu_vector];
  K_store = [K_store K_opt];
  x_est_store(1:IMM.states_len, 1:IMM.n_models, t) = [x_est];
  T_R_store = [T_R_store, T_R];
  T_G_store = [T_G_store, T_G];

end

%         _       _   
%   _ __ | | ___ | |_ 
%  | '_ \| |/ _ \| __|
%  | |_) | | (_) | |_ 
%  | .__/|_|\___/ \__|
%  |_|                

figure(); hold on; grid on;
plot(K_store, 'LineWidth', 2)
yline(K_real, 'LineWidth', 2)
title('$K_{opt}$ estimation')

figure(); hold on; grid on;
for i=1:IMM.n_models
  plot(mu_store(i, :), 'DisplayName', ['Model ', num2str(i)]);
end
legend()
title('Probabilities of the models')

figure(); hold on; grid on;
for i=1:IMM.n_models
  plot(x_est_store(i, :), 'DisplayName', ['Model ', num2str(i)]);
end
plot(omega_store, 'LineWidth', 2, 'DisplayName', 'Real')
plot(omega_tilde_store, '--', 'LineWidth', 2, 'DisplayName', 'Estimation mean')
legend()
title('Rotational speed')

figure(); hold on; grid on;
plot(T_R_store, 'LineWidth', 2, 'DisplayName', 'Aero torque')
plot(T_G_store, 'LineWidth', 2, 'DisplayName', 'Feedback torque')
legend()
title('Torques')

%               _         ___ __  __ __  __ 
%    __ _  __ _(_)_ __   |_ _|  \/  |  \/  |
%   / _` |/ _` | | '_ \   | || |\/| | |\/| |
%  | (_| | (_| | | | | |  | || |  | | |  | |
%   \__, |\__,_|_|_| |_| |___|_|  |_|_|  |_|
%   |___/                                   

function [x_est, P_est, mu_vector, omega_tilde, rho, K_opt] = gain_IMM(x_est, P_est, mu_vector, IMM, I_eq, dt, z, generator, B_eq)
  % x_est -> (states_len, n_models) items containing the state estimation
  % P_est -> (states_len, states_len, n_models) containing the covariance
  
  n_models = IMM.n_models; % number of models to test
  states_len = IMM.states_len; % length of the states
  measure_len = IMM.measure_len; % length of the measurement
  Pi = IMM.Pi; 
  R = IMM.W;
  Q = IMM.Q;
  K_vector = IMM.K_vector;
  mu = zeros(n_models);
  % Vector of parameters
  param = IMM.param;

  % Assamble the model cell
  model = cell(n_models, 1);
  for j=1:n_models
    param{2} = IMM.rho_vector(j); % rho
    model{j} = struct('rho', IMM.rho_vector(j), 'K_opt', IMM.K_vector(j), 'x_est', x_est(1:states_len, j), 'P_est', P_est(1:states_len,1:states_len,j), 'z_hat', zeros(measure_len, 1), 'P_prior', P_est(1:states_len,1:states_len,j), 'Lambda', 0, 'mu', mu_vector(j), 'mu_new', 0, 'x_tilde', zeros(states_len, 1), 'P_tilde', zeros(states_len));
  end
  
  % Filtering
  u = [0]; % input
  model = filtering(model, z, u, param, dt);
  
  % Mode probability updating
  model = probability_updating(model, z, param, dt);
  [rho, K_opt] = weight_models(model); % weighted mean of the filters

  % State combination
  x_est_g = zeros(states_len, 1); % global state
  [x_est_g(1:states_len, 1), ~] = state_combination(model);
  omega_tilde = x_est_g(1); % estimation of the rotational speed

  % Filter interaction
  [model, model_prob] = filter_interaction(model, Pi);
  
  % Assemble the model cell
  for j=1:n_models
    x_est(1:states_len, j) = model{j}.x_est(1:states_len);
    P_est(1:states_len,1:states_len,j) = model{j}.P_est(1:states_len, 1:states_len);
    mu_vector(j) = model{j}.mu_new(1);
  end
  
end
  
  %   _____ _ _ _            _             
  %  |  ___(_) | |_ ___ _ __(_)_ __   __ _ 
  %  | |_  | | | __/ _ \ '__| | '_ \ / _` |
  %  |  _| | | | ||  __/ |  | | | | | (_| |
  %  |_|   |_|_|\__\___|_|  |_|_| |_|\__, |
  %                                  |___/ 
  
  function model = filtering(model, z, u, param, dt)
    n_models = size(model, 1);
    states_len = size(model{1}.x_est, 1);
    for j=1:n_models
      param{2} = model{j}.rho;
      % compute the linearized matrices
      matrices = compute_matrices(model{j}.x_est(1:states_len), param, dt);
      matrices.R = param{4};
      matrices.Q = param{5};
  
      [model{j}] = EKF(model{j}, z, u, param, matrices, dt);
    end
  end

  %   _____ _  _______ 
  %  | ____| |/ /  ___|
  %  |  _| | ' /| |_   
  %  | |___| . \|  _|  
  %  |_____|_|\_\_|    
                     
  function [model] = EKF(model, z, u, param, matrices, dt)
    % y -> measurement of the rotational speed [rad/s]
    
    x_est = model.x_est;
    P_est = model.P_est;
    states_len = size(x_est, 1);
    measure_len = size(z, 1);

    F = matrices.F;
    H = matrices.H;
    Q = matrices.Q;
    R = matrices.R;
    G = matrices.G;
  
    % Prediction
    x_est = dynamic(x_est(1:states_len), u, param, dt);
    P_est = F*P_est*F' + G*Q*G';
    P_prior = P_est;
    states_len = size(x_est, 1); % number of states
  
    % Update
    [z_hat] = measurement_est(x_est); % estimation of the measurement using the sensor's model
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


  %    ____                            _                         _     
  %   / ___|___  _ __ ___  _ __  _   _| |_ ___   _ __ ___   __ _| |_   
  %  | |   / _ \| '_ ` _ \| '_ \| | | | __/ _ \ | '_ ` _ \ / _` | __|  
  %  | |__| (_) | | | | | | |_) | |_| | ||  __/ | | | | | | (_| | |_ _ 
  %   \____\___/|_| |_| |_| .__/ \__,_|\__\___| |_| |_| |_|\__,_|\__(_)
  %                       |_|                                          
  
  function matrices = compute_matrices(x_est, param, dt)
  
    I_eq = param{1}; % Inertia
    rho = param{2}; % Density
    K_opt = param{3}; % Power coefficient
    B_eq = param{6};
    R = param{8}; % Rotor radius
    cp_GE = param{9}; % power coefficient
    rho_cp_R2 = rho*cp_GE*R^2;

    omega = x_est(1);
    V0 = x_est(2);

    F = [1+dt/I_eq*(-1/2*rho_cp_R2*pi*V0^3/omega^2-2*K_opt*omega), (3*dt*pi*rho_cp_R2*V0^2)/(2*I_eq*omega); 0, 1];
    H = [1, 0; 0 1];
    G = [dt/I_eq 0; 0 1];

    matrices = struct('F', F, 'H', H, 'G', G);
  end

  
  %   ____                              _      
  %  |  _ \ _   _ _ __   __ _ _ __ ___ (_) ___ 
  %  | | | | | | | '_ \ / _` | '_ ` _ \| |/ __|
  %  | |_| | |_| | | | | (_| | | | | | | | (__ 
  %  |____/ \__, |_| |_|\__,_|_| |_| |_|_|\___|
  %         |___/                              
  
  function x_new = dynamic(x_est, u, param, dt)
  
    states_len = size(x_est, 1);
    I_eq = param{1}; % Inertia
    rho = param{2}; % Density
    K_opt = param{3}; % Power coefficient
    R = param{8}; % Rotor radius
    cp_GE = param{9}; % power coefficient
    rho_cp_R2 = rho*cp_GE*R^2;

    omega = x_est(1);
    V0 = x_est(2);

    x_new = zeros(states_len,1);
  
    x_new(1) = omega + dt/I_eq*(1/2*pi*rho_cp_R2*V0^3/omega - K_opt*omega^2 ); % [rad/s] rotational speed
    x_new(2) = V0; % [m/s] wind speed 
  end
  


  
  %   ____            _                       _       _   _             
  %  |  _ \ _ __ ___ | |__    _   _ _ __   __| | __ _| |_(_)_ __   __ _ 
  %  | |_) | '__/ _ \| '_ \  | | | | '_ \ / _` |/ _` | __| | '_ \ / _` |
  %  |  __/| | | (_) | |_) | | |_| | |_) | (_| | (_| | |_| | | | | (_| |
  %  |_|   |_|  \___/|_.__/   \__,_| .__/ \__,_|\__,_|\__|_|_| |_|\__, |
  %                                |_|                            |___/ 
  
  function model = probability_updating(model, z, param, dt)
    n_models = size(model, 1);
    states_len = size(model{1}.x_est, 1);
    
    for j=1:n_models
      % compute the linearized matrices
      matrices = compute_matrices(model{j}.x_est(1:states_len), param, dt);
      matrices.R = param{4};
      matrices.Q = param{5}; 

      model{j}.Lambda = zeros(1);
      zeta = zeros(states_len, 1);

      H = matrices.H;
      R = matrices.R;
      G = matrices.G;
      z_hat = model{j}.z_hat;
      P_prior = model{j}.P_prior;
      
      S = zeros(states_len, states_len);
      S = H*P_prior*H' + R;
      zeta = z - z_hat;
%       tmp = 1/sqrt(2*pi*max(det(S), 1e-80))*exp(-1/2*zeta'*inv(S)*zeta);
%       tmp = max(tmp, 1e-80);
      tmp = 1/sqrt(2*pi*det(S))*exp(-1/2*zeta'*inv(S)*zeta);
      model{j}.Lambda = tmp(1);
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
  
  %                                       _ 
  %   _ __ _____   ___ __  _ __ _ __   __| |
  %  | '_ ` _ \ \ / / '_ \| '__| '_ \ / _` |
  %  | | | | | \ V /| | | | |  | | | | (_| |
  %  |_| |_| |_|\_/ |_| |_|_|  |_| |_|\__,_|
                                          
  function [z_hat] = measurement_est(x_est)
      % Measurement modelling the sensor
      
      z_hat = x_est;

  end

  %    ___        _               _   
  %   / _ \ _   _| |_ _ __  _   _| |_ 
  %  | | | | | | | __| '_ \| | | | __|
  %  | |_| | |_| | |_| |_) | |_| | |_ 
  %   \___/ \__,_|\__| .__/ \__,_|\__|
  %                  |_|              
  function [rho, K_opt] = weight_models(model)
    % This function computes the local estimation of the global centroid, so where each agent thinks the network centroid is
    
    n_models = size(model, 1);
  
    rho = 0;
    K_opt = 0;
    for j=1:n_models
      rho = rho + model{j}.mu*model{j}.rho;
      K_opt = K_opt + model{j}.mu*model{j}.K_opt;
    end
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
 
