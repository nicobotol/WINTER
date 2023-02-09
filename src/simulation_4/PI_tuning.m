function [kp_d, ki_d, kd_d] = PI_tuning(G_noR, tA)
% This function tunes a PI controller
% G_noR -> open-loop TF
% tA -> raising time [s]
syms ki tau_i
s = tf('s');                    % define s as complex variable

R = ki*(1 + tau_i*s)/s;         % controller TF
GH = R*G_noR;                   % close loop TF

tau_i_d = min(pole(G_noR));               % designed tau_i = kp/ki
GH = simplify(subs(GH, tau_i, tau_i_d)); % TF after zero-pole cancellation

% CASE 1: assuming Phase margin >= 75°
omega_BP = 5/tA; % [rad/s]
[Pm, ki_d, kp_d] = shaping(GH, omega_BP, tau_i_d);

% CASE 2: case 1 fails and so try to shape the PI iteratively
if Pm < 75*pi/180
  for i=1:i_max 
    omega_BP = 500/(Pm*tA); % new guess on the band pass
    [Pm, ki_d, kp_d] = shaping(GH, omega_BP, tau_i_d);
    if omega_BP - 500/(Pm*tA) < fake_zero % convergence
      break
    else                                  % reduce the phase margin
      Pm = Pm - (75 - 60)*pi/180/i_max;
    end
  end
end

kd_d = 0.0; % derivative gain
end


function [Pm, ki_d, kp_d] = shaping(GH, omega_BP, tau_i_d)
% This function shapes the proportional and integral gain according to the 
% guess on the bandpass frequency
norm_TF = norm(evalfr(GH, 1j*omega_BP));
ki_d = solve(norm_TF == 1, ki); % integral gain

kp_d = ki_d*tau_i_d;

% rewrite the final close loop TF
GH_d = subs(GH, [ki, kp], [ki_d, kp_d]);

% check the phase margin
[~, Pm] = margin(GH_d);
if Pm < pi/3 % unstable 
  disp('The TF is unstable');
else
end

end