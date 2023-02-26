close all
clear 
clc

wind = load('16ms.csv');

figure()
plot(wind(:, 1), wind(:, 2))
grid on

std(wind(:, 2))