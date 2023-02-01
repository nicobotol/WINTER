%% Compare cP computed with the loookup table and with the formulas

clear;
close all;
clc;

% load the parameters
parameters

% redefine a range for the pitch and the lambda
pitch_r = deg2rad([-3, 12]);
lambda_r = [2, 10];
pitch_i = 6;
lambda_i = 10;
pitch_v = linspace(pitch_r(1), pitch_r(2), pitch_i);
lambda_v = linspace(lambda_r(1), lambda_r(2), lambda_i);

cP = zeros(pitch_i, lambda_i);
for l = 1:lambda_i
  for p = 1:pitch_i
    lambda = lambda_v(l);
    pitch = rad2deg(pitch_v(p));
    lambda1 = 1/(1/(lambda + 0.089) - 0.035/(pitch^3 + 1));
    cP(p, l) = 0.5*(116/lambda1 - 0.4*(pitch - 5))*exp(-16.5/lambda1);
  end
end

LegS = cell(1, pitch_i);
cP_interp = zeros(lambda_i, 1);
figure(1)
hold on
for p = 1:pitch_i
%   plot(lambda_v, cP(p, :));
  cP_interp = interp2(lambda_vector, pitch_vector, cP_store, lambda_v(:), pitch_v(p));
  plot(lambda_v, cP_interp, '--');
  LegS{p} = ['lambda {',num2str(rad2deg(pitch_v(p))),'}'];
%  LegS{2*p - 1} = ['lambda {',num2str(pitch_v(p)),'}'];
%   LegS{2*p} = ['lambda {',num2str(pitch_v(p)),'}'];
end
legend(LegS, 'location', 'best');
xlabel('\theta [rad]')
ylabel('c_P [-]')
grid on