model = cell(2,1);
model{1}.K_opt = 10;
model{2}.K_opt = 20;
model{1}.P_est = 1;
model{2}.P_est = 1;

wls_K_opt(model)

function K_opt = wls_K_opt(models)
    % This function computes the local estimation of the global centroid, so where each agent thinks the network centroid is
    
      n_models = size(models, 1);
    
      % build the measurement matrix
      H = kron(ones(n_models,1), 1);
    
      Z = []; % measurements
      C = []; % covariance matrix
      for j=1:n_models
        % rearrange the measurements in a vector
        Z = [Z; models{i}.K_opt];
        
        % build the covariance matrix
        C = blkdiag(C, models{i}.P_est);
      end
      % compute the WLS
      K_opt = inv(H'*inv(C)*H)*H'*inv(C)*Z;

     
end

function K_opt = compute_K_opt(model)
% This function computes the local estimation of the global centroid, so where each agent thinks the network centroid is

  n_models = size(model, 1);

  K_opt = 0;
  for j=1:n_models
    K_opt = K_opt + model{j}.mu*model{j}.K_opt;
  end
      
     
end