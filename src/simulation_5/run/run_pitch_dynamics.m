function [blade] = run_pitch_dynamics(blade, i)
% This functions sets different transfer functions to model the pitching
% actuation mechanism

actuator_nominal = tf(blade.omegap^2, [1 2*blade.zetap*...
  blade.omegap blade.omegap^2]); % nominal tf of the pitch actuator

  if i == 1 % with reference blade pitching dynamic
    blade.actuator_dynamic = actuator_nominal;
    

  elseif i == 2 % with faster blade pitch dynamic
%     omega = blade.omegap*0.5; % [rad/s];
%     zeta = 0.5;
% 
%     actuator_fast =  tf(omega^2, [1 2*zeta*omega omega^2]);
%     blade.actuator_dynamic = actuator_fast;
    blade.actuator_dynamic = tf(1, [1 1e-1]);

  elseif i ==3  % with slower blade pitch dynamic
    zeta = blade.zetap;

    actuator_fast = tf(1, [1 1e-1]);
    blade.actuator_dynamic = actuator_fast;
    
    omega = blade.omegap*10;
    actuator_slow =  tf(omega^2, [1 2*blade.zetap*omega omega^2]);
    blade.actuator_dynamic = actuator_slow;

    % plot     
    tf_plot = [actuator_nominal, actuator_fast, actuator_slow];
    title_plot = "Pitching mechanism dynamic";
    legends_plot = ["Sim. 1", "Sim. 2", "Sim. 3"];
    plot_name = "pitching_dynamic_TF";
    bode_plot(tf_plot, legends_plot, title_plot, plot_name)
  end

end

function actuator_nominal = rescale_tf(actuator_nominal, gain)
  pol = actuator_nominal.den{1}; % pole of the actuator's tf
  pol_new = [pol(1) pol(2)/gain pol(3)/gain^2];
  actuator_nominal = tf(pol_new(3), pol_new);
end
