function [in, s_time, out_store] = run_generator_step()
% This function tests the step resonse of the generator

parameters;

open_system('simulink\PMSM_generator.slx');     % open the model
in = Simulink.SimulationInput('PMSM_generator');  % set simulation parameters

s_time = simulation.stop_time; % time for the simulation

end