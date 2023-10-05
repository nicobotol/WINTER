
actuator_nominal = tf(blade.omegap^2, [1 2*blade.zetap*blade.omegap blade.omegap^2]); % nominal tf of the pitch actuator
blade.actuator_dynamic = actuator_nominal;

omega = blade.omegap*5;
actuator_fast = tf(omega^2, [1 2*blade.zetap*omega omega^2]);
% actuator_fast = tf(1, 1);

omega = blade.omegap/5;
actuator_slow =  tf(omega^2, [1 2*blade.zetap*omega omega^2]);
blade.actuator_dynamic = actuator_slow;
