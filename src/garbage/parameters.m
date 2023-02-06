%% PARAMETERS

%%
V0 = 15; % windspeed (m/s)
omega = 1.2; % rotational speed (rad/s)
P = 10e6; % power where to pitch

%% From here on, nothing should be changed
R = 89.17; % rotor radius (m)
A = R^2 * pi; % rotor area (m^2)
B = 3; % number of blades (#)
rho = 1.225;  % air density (kg/m^3)
V0_cutin = 4; % cut in wind velocity (m/s)
V0_cut_out = 25; % cut out wind velocity (m/s)
P_rated = 10.64e6; % rated power (W)
lambda_opt = 7.857; % optimal lambda
cp_opt = 0.465; % optimal cp
omega_max = 1.01; % maximum rotational speed (rad/s)
velocity_pitch = [11.55:1:20]; % range of velocities where to compute the pitch
velocity_item = size(velocity_pitch, 2);

pitch_item = 50; % number of division for the pitch 
pitch_range = [-deg2rad(20) deg2rad(20)];  % (rad) range in within look for picth angle 

stall_lim = -1*pi/180;
feather_lim = 4*pi/180;

% first guess for a and a_prime
a_guess = 0;
a_prime_guess = 0.1;

% file with airfoils parameter 
filenames = [ "airfoil_data\cylinder", "airfoil_data\FFA-W3-600",  ...
  "airfoil_data\FFA-W3-480", "airfoil_data\FFA-W3-360", ...
  "airfoil_data\FFA-W3-301", "airfoil_data\FFA-W3-241"];
% file with blade parameters
blade_filename = "airfoil_data\bladedat.txt";
% t/c ratio (they are in the same order provided as the file uploaded)
thick_prof = [100 60 48 36 30.1 24.1];

i_max = 300;  % maximum number of accepted iterations
fake_zero = 1e-8; % thershold for exiting the loop

% pitch parameters from DTU reference turbine pag 33
velocity_reference = [4:1:25];
pitch_reference = [2.751 1.966 0.896 0.000 0.000 0.000 0.000 0.000 ... 
  4.502 7.266 9.292 10.958 12.499 13.896 15.200 16.432 17.618  18.758 ...
  19.860 20.927 21.963 22.975 ];

font_size = 25; 
line_size = 2;
