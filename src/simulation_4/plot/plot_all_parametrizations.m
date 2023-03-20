% This function plots the oarametrization between two different signals

addpath('\..')
addpath('lookup\')

% torque as function of rotational speed
plot_parametrization('fig_torque_vs_omega',out_store,'omega_R','T_G', ...
  '$\omega_{R}$ [rpm]','T [MNm]','Generator torque',pi/30,1e6,date_fig)