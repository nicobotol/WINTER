clear; close all; clc;
parameters;
tuning_controllers;
L = generator.Lq;
Rs = generator.Rs;
p = generator.p;
J = I_eq;
B = B_eq;
Lambda =  generator.Lambda;

A = [ -Rs/L Lambda*p/L 0;
      -3/2*p*Lambda/J -B/J 0;
      0 1 0];
B = [-1/L 0 0]';% 0 0 0; 0 0 0]';
E = [0 1/J 0]';
C = [1 0 0];
D = [0];

K = lqr(A,B,eye(3),1);
S = eye(3);
R = 1;
Sf = eye(3);
dt  = 1e-3;

% stop_time = 3;
% time = dt:dt:stop_time;
% states_len = 3;
% inputs_len = 1;
% K = lqr_gain(A, B, S, R, stop_time/dt, Sf, states_len, inputs_len, dt);
% for i = 1:length(time)
%   K_gain(:,i) = K(:,:, i);
% end
omega_R0 = 0;
stop_time = simulation.stop_time;
Q = 100*[1 0 0;
  0 1 0;
  0 0 1];
R = 0.1;
K_lqr = lqr(A, B, Q, R);
integrator_initial_condition = [0 omega_R0 0]; % [iq, omega, theta] 

open_system('PMSM_generator');                    % open the model
in = Simulink.SimulationInput('PMSM_generator');  % set simulation parameters
out = sim(in, 'ShowProgress','on'); 

% sys1 = ss(A-B*K_lqr,B,C,D);
% step(sys1)