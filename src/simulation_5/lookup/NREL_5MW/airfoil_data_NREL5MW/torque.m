%% Compute the reference torque

references = load('DTU_10MW_reference.txt'); 

len = size(reference, 1);
rho = 1.225;          % air density [kg/m^3]
cp = 0.476;           % max power coefficient
omega_rated = 1.005;  % rated rotational speed [rad/s]
lambda = 7.5;         % optimal TSR
R = 89.17;            % rotor radius [m]
K_opt = rho*cp*pi*R^5/(2*lambda^3);

for i=1:len
  if references(i, 1) < 11.4  % below rated
    references(i, 7) = K_opt*(references(i, 3)*pi/30)^2;
  else                        % above rated
    references(i, 7) = references(i, 6)*1e3/(references(i, 3)*pi/30);
  end
end

writematrix(references,'DTU_10MW_reference.txt','Delimiter','tab')