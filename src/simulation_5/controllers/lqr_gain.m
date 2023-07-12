function K = lqr_gain(A, B, S, R, T, Sf, states_len, input_len, t)

  % A -> state matrix
  % B -> input matrix
  % minimization function: cost = x_f*Sf*x_f + sum(i=1:N) (x(:,i)-target)*S*(x(:,i)-target)' + u(:,i)*R*u(:,i)'
  % S -> cost state matrix
  % Sf -> final cost state matrix
  % R -> cost input matrix
  % T -> time horizon
  % n -> number of agents
  
  K = zeros(input_len, states_len, T);
  
  % Backward cycle
  P = zeros(states_len, states_len, T+1);
  P(:, :, T+1) = Sf;
  % LQR algorithm
  for j=T:-1:t
    P(:,:,j) = S + A'*P(:,:,j+1)*A - A'*P(:,:,j+1)*B*inv(R + B'*P(:,:,j+1)*B)*B'*P(:,:,j)*A;
  end
  
  % Gain matrix calculation
  for t=1:T
    K(:,:,t) = inv(R + B'*P(:, :, t + 1)*B)*B'*P(:, :, t + 1)*A;
  end
  
  end