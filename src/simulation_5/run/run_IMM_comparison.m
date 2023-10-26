function IMM = run_IMM_comparison(IMM, i)
  % In this function set the value for enabling or disabling the IMM algoritm. In the simulink model this is done by a gain block. When IMM.enable=1 then the algorithm is enable, while when IMM.enable=0 then the constant gain is applied

  parameters;

  if rem(i, 2) == 1 % enable IMM
    IMM.enable = 1;
    rng(i);
  elseif rem(i, 2) == 0 % disable IMM and use constant gain
    IMM.enable = 0;
    rng(i-1);
  end

  IMM.sigma_omega = IMM.sigma_gain(i)*omega_rated*0.05/3; % fixed as 5% of the nominal value
  IMM.sigma_rho = IMM.sigma_gain(i)*rho*0.05/3; % fixed as 5% of the nominal value
  IMM.sigma_R = IMM.sigma_gain(i)*(4/3); % assuming a  deflection of 4 meters
  IMM.sigma_V0_rated = IMM.sigma_gain(i)*V0_rated*0.15/3; % fixed as 15% of the nominal value 
  IMM.sigma_theta = IMM.sigma_gain(i)*1*pi/180/3; % assuming 1 deg of uncertainty


end