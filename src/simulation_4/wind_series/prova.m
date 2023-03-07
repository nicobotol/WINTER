clear 
close all
clc
% This function produces a wind series provided:
% V10     -> averaged wind speed blowing in 10 minutes [m/s]
% V10_std -> std of the V10 [m/s]
% f       -> sampling freq. [Hz] (i.e. f=1/delta_t, delta_t sampling time)
% h       -> height above ground level [m]
% T       -> total time [s]

% This function is based on Hansen M.O.L. "Aerodynamics of wind turbines",
% 2015
%%
h = 150;
% l parameter based on the height
if h <= 30 
  l = 20;
else
  l = 600;
end

fs = 50; % sample frequency
V10 = 10;
V10_std = 0.50;
T = 1000; % second to generate

delta_ts = 1/fs;  % sampling time [s]
N = T/delta_ts;   % number of samples

I = V10_std/V10;  % turbulence intensity

t = delta_ts*[1:1:N]; % [s] time  
cos_v = zeros(N, 1);  % vector of cosines
p_sum = 0;            % partial sum
PSD_store = zeros(N/2, 1);
for n = 1:N/2
  f_n = n/T;                                  % frequency [Hz]
  PSD = I^2*V10*l/(1 + 1.5*f_n*l/V10)^(5/3);  % PSD
  PSD_store(n) = PSD;
  phi_n = rand(1)*2*pi;                       % Random phase [rad]
  cos_v = cos(2*pi*f_n.*t - phi_n);           % vector of cosines
  p_sum = p_sum + sqrt(2*PSD/T)*cos_v;        % partial sum
end

sum(PSD_store)/T

u = V10 + p_sum;  % add the mean to the windspeed [m/s]
%%
% Rescale the std
u1 = V10_std*u + (1 - V10_std)*V10;

figure()
plot([1/fs:1/fs:T], u)
hold on
plot([1/fs:1/fs:T], u1)
legend('u', 'u1', 'location','southeast')
%%
fs2 = 10000;
window = hamming(256);  % Window for spectral estimation
noverlap = 128;  % Number of samples to overlap between segments
nfft = 1024;  % Number of points for FFT
[Pxx,f] = periodogram(u,hamming(length(u)),length(u),fs2);
trapz(f, Pxx)
sum(Pxx)*(f(2)-f(1))

figure()
plot([1:1:N/2]/T, PSD_store)
hold on
plot(f, Pxx)
