close all
clear 
clc

wind_l = load('8ms.csv');
wind(:,:) = wind_l(wind_l(:,2)>3.3, :);

figure()
plot(wind(:, 1), wind(:, 2))
grid on

fprintf('mean: %3.3f [m/s] \n', mean(wind(:, 2)));
fprintf('std: %3.3f [m/s] \n', std(wind(:, 2)));
fprintf('Theoretical std: %3.3f [m/s] \n', 0.16*(8.105*0.75 + 3.8));
