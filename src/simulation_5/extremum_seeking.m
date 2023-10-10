function extremum_seeking()
%% This function computes the optimal power-speed gain

switch mode
case 0 % prevoius step used HSC

  vc = logic_0();
  switch vc
  case 0
    if abs(delta_P/delta_omega) < A_th
      omega_new = omega;
    else
      omega_new = omega + beta*(delta_P/delta_omega);
    end
    mode_new = 0;
  case 1
    omega_new = (P/k_opt)^(1/3);
    mode_new = 1;
  otherwise
    error('Case not defined')
  end

case 1 % previous step used PSF

  vc = logic_1();
  switch vc
  case 0
    omega_new = omega + delta_P/A_th;
    mode_new = 0;
  case 1
    
  otherwise
    error('Case not defined') 
  end

otherwise % error
  error('Case not defined')
end



end