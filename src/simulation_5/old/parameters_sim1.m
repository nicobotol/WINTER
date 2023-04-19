%% Parameters for simulation 1

wind.ramp_WS_start = 4;         % wind speed at the start of the ramp [m/s]
wind.ramp_WS_stop = 10;         % wind speed at the stop of the ramp [m/s]
wind.ramp_time_start = [0 0 0];  % time speed at the start of the ramp [s]
wind.ramp_time_stop = [10 20 40];  % time speed at the stop of the ramp [s]

rotor.omega_r = 0.33;       % initial rotational speed [rad/s]

simulation.stop_time = [40 60 80];  % max time to investigaste [s]
simulation.plot_step = simulation.plot_time/simulation.time_step;