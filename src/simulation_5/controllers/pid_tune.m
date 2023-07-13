function [kp, ki, kd, tau_d1] = pid_tune(G, omega)
% Tune a PI controller for the open loop TF G with crossover frequency
% omega

syms s ki

pol = poles(G, s);                        % pole of the TF
pol_abs = eval(abs(pol));
pol_abs_sort = sort(pol_abs, 'ascend');
% tau_p = 1/pol_abs_sort(1);                % 1st zero 
tau_p=1/10^((log10(pol_abs_sort(1))+log10(pol_abs_sort(end)))/2); %1st zero 
tau_d = 1/pol_abs_sort(end);              % 2nd zero 
tau_d1 = 1/(1e1*omega);                    % high freq pole
R = ki/s*(1 + s*tau_p)*(1 + s*tau_d)/(1 + s*tau_d1); % regulator
GH = G*R;                                  % series of regulator and system
GH_mag = abs(subs(GH, s, 1j*omega));      % mag at the crossover frequency
ki = solve(GH_mag == 1, ki);              % integral gain
ki = 1*eval(ki);
kp = ki*(tau_d + tau_p); % 0.5*ki*(tau_d + tau_p);
kd = ki*tau_d*tau_p;     % 0.5*ki*tau_d*tau_p;

fprintf('ki = %f\n', ki);
fprintf('kp = %f\n', kp);
fprintf('kd = %f\n', kd);


% ku = 22.06;
% % kp=ku;
% Tu=2*(0.508535-509213);
% kp=0.5*ku;
% Ti = 0.83*Tu;
% Td = 0.125*Tu;
% ki = 0.54*ku/Tu;
% kd =0.075*ku*Tu;
% kd=0;ki=0;

end