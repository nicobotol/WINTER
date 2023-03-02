function [ki, tau_p, tau_d, tau_d1] = pid_tune(G, omega)
% Tune a PI controller for the open loop TF G with crossover frequency
% omega

syms s ki

pol = poles(G, s);                        % pole of the TF
pol_sort = sort(pol,'descend');
tau_p = eval(1/abs(pol_sort(end)));       % 1st zero
tau_d = eval(1/abs(pol_sort(end - 1)));   % 2nd zero
tau_d1 = 1/(15*omega);                    % high freq pole
R = ki/s*(1 + s*tau_p)*(1 + s*tau_d)/(1+s*tau_d1); % regulator

GH = G*R;                                 % series of regulator and system
GH_mag = abs(subs(GH, s, 1j*omega));      % mag at the crossover frequency
ki = solve(GH_mag == 1, ki);              % integral gain
ki = eval(ki);

end