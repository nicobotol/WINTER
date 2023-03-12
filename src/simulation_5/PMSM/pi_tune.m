function [ki_r, kp_r] = pi_tune(G, omega)
% Tune a PI controller for the open loop TF G with crossover frequency
% omega

syms s ki

pol = poles(G, s);                  % pole of the TF
pol_abs = eval(abs(pol));
pol_abs_sort = sort(pol_abs, 'ascend');
tau = 1/pol_abs_sort(end);          % const. of the controller tau = kp/ki
R = (1 + s*tau)*ki/s;               % regulator TF
GH = G*R;                           % series of regulator and system
GH_mag = abs(subs(GH, s, 1j*omega));% magnitude at the crossover frequency
ki = solve(GH_mag == 1, ki);        % integral gain
kp = tau*ki;                        % proportional gain
kp_r = eval(kp);                    % evaluate the gains
ki_r = eval(ki);

fprintf('ki = %f\n', ki);
fprintf('kp = %f\n', kp);

end