function [ki_r, kp_r, kd_r] = pid_tune(G, omega)
% Tune a PI controller for the open loop TF G with crossover frequency
% omega

syms s ki

pol = poles(G, s);                  % pole of the TF
pol_abs = sort(pol,'ascend');
tau = eval(1/abs(max(pol)));        % const. of the controller tau = kp/ki
GH = G*(1 + s*tau)*ki/s;            % series of regulator and system
GH_mag = abs(subs(GH, s, 1j*omega));% magnitude at the crossover frequency
ki = solve(GH_mag == 1, ki);       % integral gain
kp = tau*ki;                        % proportional gain
kp_r = eval(kp);                    % evaluate the gains
ki_r = eval(ki);
tau_d = eval(1/abs(pol_abs(end-1)));    % const. of the derivative tau = kd/ki
kd_r = tau_d*ki_r;
end