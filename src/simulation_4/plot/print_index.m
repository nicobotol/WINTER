function [print_range] = print_index(stop_time, plot_time, time_step, ...
  WS_len)
% This fucntion provides the indeces of the element to be printed
% stop time -> time from the beggining of the simulation where to stop the
%               printing [s]
% print_time -> total time to be printed [s]
% time_step -> time step of the simulation [s]

print_rage = zeros(WS_len, 2);
print_range(:, 1) = (stop_time - plot_time)/time_step;
print_range(:, 2) = stop_time/time_step;
print_range = ceil(print_range);
