% This function wants to test the effectivness of the IMM for the estimation of the cp coefficient

close all;
clear; 
clc;
rng(4)

parameters;
n_models = 3; % number of models
states_len = 2;
V0 = 5; % wind speed
model = cell(n_models, 1);
prob_transition = 0.99;
Pi = prob_transition*eye(n_models, n_models) + (1 - prob_transition)/(n_models-1)*(ones(n_models, n_models)-eye(n_models)); % mode transition matrix
t_tot = 1000; % number of time sample to investigate
R = (0.05/3)^2; % measurement noise
Q = [(1e-7/3)^2 0; 0 (V0/100/3)^2]; % process noise
mu = zeros(n_models);
model{1}.cp = [0.1];
model{2}.cp = [0.4];
model{3}.cp = [0.6];
% model{4}.cp = 0.20;
cp_real = 0.4; % intial cp

% Vector of parameters
param{1} = I_eq;
param{2} = rho;
param{3} = rotor.R;
param{5} = R;
param{6} = Q;

% Initialize the states
x_real = [0.415; V0]; % state
P_real = 1000*eye(states_len); % covariance
x_est = x_real;
P_est = P_real;

x_real_store = [];
x_est_store = [];
cp_real_store = [];
model_prob = [];
K_store = [];

% Initialize the models
for j = 1:n_models
  model{j}.x_est = x_real; % state
  model{j}.P_est = 10*ones(states_len); % covariance
  model{j}.mu = 1/n_models; % mode probability
  model{j}.x_est_store = [];
  model{j}.cp_est_store = [];
  model{j}.mu_store = [];
end

for t=1:t_tot
  pitch = 0; % pitch angle
  T_G = 0.5*rho*cp_real*pi*param{3}^2*x_real(2)^3/x_real(1); % torque generated by the generator

  % Propagate the dynamic
  lambda =  param{3}*x_real(1)/x_real(2); % TSR
  cp_real = max(0, interp2(lambda_vector, pitch_vector, lookup_cP, lambda, pitch));
  param{4} = cp_real;
  x_real = dynamic(x_real, T_G, param) + mvnrnd(zeros(states_len, 1), param{6})' ; % propagate the dynamic and add noise
  if t < t_tot/2
    x_real(2) = V0; % wind speed (ADD NOISE HERE)
  else
    x_real(1) = 1.01;
    x_real(2) = 4*V0;
  end
  x_real_store = [x_real_store x_real];
  z = measurement(x_real, param{5}) + mvnrnd([0], R)'; % rotational speed measurement [rad/s]

  % Filtering
  model = filtering(model, z, T_G, param);

  % Output estimation
  model = output_estimation(model, pitch, param, lookup_cP, lambda_vector, pitch_vector);
  model = hystory(model);

  % Mode probability updating
  model = probability_updating(model, z, param);

  % State combination
  [x_est, P_est] = state_combination(model);

  % Filter interaction
  [model, model_prob] = filter_interaction(model, Pi, model_prob);
  [~, idx] = max(model_prob(:, end)); % model with the highest probability

  % Compute the gain
  K = 0.5*model{idx}.cp*pi*param{2}*param{3}^2*x_est(2)^3/x_est(1)^3;
  K_store = [K_store K];

  % Save the hystory of the estimation
  x_est_store = [x_est_store x_est];
  cp_real_store = [cp_real_store cp_real];

end

figure();hold on;grid on;
for i=1:n_models
  yline(model{i}.cp, 'DisplayName', ['Model ', num2str(i)], 'LineWidth', line_width)
end
plot(cp_real_store(1, :), 'DisplayName', 'Real', 'LineWidth', line_width)
title('$c_P$')
xlabel('Time')
ylabel('$c_P$')
legend()

figure();hold on;grid on;
for i=1:n_models
  plot(model{i}.x_est_store(1, :), 'DisplayName', ['Model ', num2str(i)], 'LineWidth', line_width)
end
plot(x_real_store(1, :), '-.', 'DisplayName', 'Real', 'LineWidth', line_width)
plot(x_est_store(1, :), '--', 'DisplayName', 'Est.', 'LineWidth', line_width)
title('$\omega$')
xlabel('Time')
ylabel('$\omega$')
legend()
% ylim([0 1])

figure(); hold on; grid on;
for i=1:n_models
  plot(model{i}.mu_store, 'DisplayName', ['Model ', num2str(i)], 'LineWidth', line_width, 'Color', color(i))
end
title('Model probability')
xlabel('Time')
ylabel('$\mu$')
legend()

figure(); hold on; grid on;
plot(K_store, 'LineWidth', line_width)
title('Gain')
xlabel('Time')
ylabel('$K$')


figure();hold on;grid on;
for i=1:n_models
  plot(model{i}.x_est_store(2, :), 'DisplayName', ['Model ', num2str(i)], 'LineWidth', line_width)
end
plot(x_real_store(2, :), '-.', 'DisplayName', 'Real', 'LineWidth', line_width)
plot(x_est_store(2, :), '--', 'DisplayName', 'Est.', 'LineWidth', line_width)
title('$V_0$')
xlabel('Time')
ylabel('$V_0$')
legend()

rms([x_est_store-x_real_store]')

%   _____ _  _______ 
%  | ____| |/ /  ___|
%  |  _| | ' /| |_   
%  | |___| . \|  _|  
%  |_____|_|\_\_|    
                   
function [model] = EKF(model, z, T_G, param, matrices)
  % y -> measurement of the rotational speed [rad/s]
  
  x_est = model.x_est;
  P_est = model.P_est;

  F = matrices.F;
  H = matrices.H;
  Q = matrices.Q;
  R = matrices.R;
  G = matrices.G;

  % Prediction
  x_est = dynamic(x_est, T_G, param);
  P_est = F*P_est*F' + G*Q*G';
  P_prior = P_est;
  states_len = size(x_est, 1); % number of states

  % Update
  [~, z_hat] = measurement(x_est, param{5}); % estimation of the measurement using the sensor's model
  Innovation = z - z_hat;
  S_Inno = H*P_est*H' + R;
  W = P_est*H'*inv(S_Inno); % kalman gain
  x_est = x_est + W*Innovation; % update state estimate
  P_est = (eye(states_len) - W*H)*P_est; % update covariance matrix


  model.x_est = x_est;
  model.P_est = P_est;
  model.z_hat = z_hat;
  model.P_prior = P_prior;
end

%   ____                              _      
%  |  _ \ _   _ _ __   __ _ _ __ ___ (_) ___ 
%  | | | | | | | '_ \ / _` | '_ ` _ \| |/ __|
%  | |_| | |_| | | | | (_| | | | | | | | (__ 
%  |____/ \__, |_| |_|\__,_|_| |_| |_|_|\___|
%         |___/                              

function x_est = dynamic(x_est, T_G, param)

  I_eq = param{1}; % Inertia
  rho = param{2}; % Density
  R = param{3}; % Radius
  cp = param{4}; % Power coefficient

  x_est(1) = x_est(1) + 1/I_eq*(1/2*rho*pi*R^2*cp*x_est(2)^3*x_est(1)^(-1) - T_G);
  x_est(2) = x_est(2);
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
  rho = param{2}; % Density
  R = param{3}; % Radius
  cp = param{4}; % Power coefficient

  F = [1-1/I_eq*(1/2*rho*pi*R^2*cp*x_est(2)^3*x_est(1)^(-2)) 1/I_eq*(3/2*rho*pi*R^2*cp*x_est(2)^2*x_est(1)^(-1));
       0 1];
  H = [1 0];
  G = [1 1];

  matrices.F = F;
  matrices.H = H;
  matrices.G = G;
end

%   _____ _ _ _            _             
%  |  ___(_) | |_ ___ _ __(_)_ __   __ _ 
%  | |_  | | | __/ _ \ '__| | '_ \ / _` |
%  |  _| | | | ||  __/ |  | | | | | (_| |
%  |_|   |_|_|\__\___|_|  |_|_| |_|\__, |
%                                  |___/ 

function model = filtering(model, z, T_G, param)
  n_models = size(model, 1);
  for j=1:n_models
    param{4} = model{j}.cp;
    % compute the linearized matrices
    matrices = compute_matrices(model{j}.x_est, param);
    matrices.R = param{5};
    matrices.Q = param{6};

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
    matrices.R = param{5};
    matrices.Q = param{6}; 
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
    x_est = x_est + model{j}.mu*model{j}.x_est;
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
                                                                            
function [model, model_prob] = filter_interaction(model, Pi, model_prob)
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
    model{j}.x_tilde = zeros(states_len, 1);
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
  tmp = zeros(n_models, 1);
  for j=1:n_models
    model{j}.mu_store = [model{j}.mu_store model{j}.mu_new];
    tmp(j) = model{j}.mu_new;
  end
  model_prob(:,end+1) = tmp;
end

%    ___        _               _               _     
%   / _ \ _   _| |_ _ __  _   _| |_    ___  ___| |_   
%  | | | | | | | __| '_ \| | | | __|  / _ \/ __| __|  
%  | |_| | |_| | |_| |_) | |_| | |_  |  __/\__ \ |_ _ 
%   \___/ \__,_|\__| .__/ \__,_|\__|  \___||___/\__(_)
%                  |_|                                

function model = output_estimation(model, pitch, param, lookup_cP, lambda_vector, pitch_vector)
  n_models = size(model, 1);
  R = param{3}; % Radius
  for j=1:n_models
    lambda = R*model{j}.x_est(1)/model{j}.x_est(2); % TSR
    model{j}.cp_est = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, pitch);
  end
end

%   _   _           _                   
%  | | | |_   _ ___| |_ ___  _ __ _   _ 
%  | |_| | | | / __| __/ _ \| '__| | | |
%  |  _  | |_| \__ \ || (_) | |  | |_| |
%  |_| |_|\__, |___/\__\___/|_|   \__, |
%         |___/                   |___/ 

function model = hystory(model)
  n_models = size(model, 1);
  for j=1:n_models
    model{j}.x_est_store = [model{j}.x_est_store model{j}.x_est];
    model{j}.cp_est_store = [model{j}.cp_est_store model{j}.cp_est];
  end
end
