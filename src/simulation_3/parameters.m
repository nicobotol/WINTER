%% PARAMETERS

% general parameters
rho = 1.225;                % air density [kg/m^3]
stop_time = 300;            % max time to investigaste [s]
font_size = 25;             % fontsize for plots
line_width = 2;             % line width for plots
i_max = 300;                % max number of accepted iterations
fake_zero = 1e-8;           % thershold for exiting the loop
a_guess = 0;                % initial guess for the BEM code
a_prime_guess = 0.1;        % initial guess for the BEM code
V0_cut_in = 4;              % cut in wind speed [m/s]
V0_cut_out = 25;            % cut out wind speed [m/s]

% rotor parameters
rotor.R = 89.17;            % rotor radius [m]
rotor.A = rotor.R^2*pi;     % rotor area [m^2]
rotor.blades = 3;           % number of blades [#]
rotor.V0_cutin = 4;         % cut in wind velocity [m/s]
rotor.V0_cutout = 25;       % cut out wind velocity [m/s]
rotor.P_rated = 10.64e6;    % rated power [W]
rotor.I = 325671;           % inertia [kgm^2]
rotor.omega_r = 0;          % initial rotational speed [rad/s]
rotor.B  = 1;               % rotational friction [kgm^2/s] (random placeholder)\

% generator parameters
generator.I = 1e3;          % generator iniertia [kgm^2]
generator.B = 1;            % rotational friction [kgm^2/s] (random placeholder)
generator.vll = 4e3;        % rated line-to-line voltage [V]
generator.is = 1443.4;      % rated stator current [A]
generator.fe = 26.66;       % rated stator frequency [Hz]
generator.p = 320;          % number of poles
generator.Ld = 1.8e-3;      % d-axis stator inductance [H]
generator.Lq = 1.8e-3;      % q-axis stator inductance [H]
generator.Rs = 64;          % stator resistance [ohm]
generator.Lambda = 19.49;   % magnet flux-linkage [Wb]
generator.tau_c = 50e-6;    % inverter time constant [s] (pag. 141 notes 'azionemanti elettrici')
generator.p_ctrl = 10;      % gain for the Ig reference
generator.k_ctrl = 0.01;    % parmter for the Iq refernce

% gearbox_parameters
gearbox.ratio = 1;          % gearbox transmission ratio 

% lambda_opt = 7.857; % optimal lambda
% cp_opt = 0.465; % optimal cp
% omega_max = 1.01; % maximum rotational speed (rad/s)
% velocity_pitch = [11.55:1:20]; % range of velocities where to compute the pitch
% velocity_item = size(velocity_pitch, 2);

% parameters for the file lookup_cp.m
pitch_range = deg2rad([-15 90]);              % range for picth angle [rad]
pitch_item = ceil(rad2deg(diff(pitch_range))); % # of guess pitch 
lambda_range = [3 9];                         % range for the TSR
lambda_item = 20;                             % # of guess pitch 
% distribute lambda and TSR in their ranges
lambda_vector = linspace(lambda_range(1), lambda_range(2), lambda_item); 
pitch_vector = linspace(pitch_range(1), pitch_range(2), pitch_item);

% load mesh for cP, cT, rated values, and pitch angle (computed in 
% lookup_cp.m and lookup_pitch.m)
if isfile('lookup_cP_theta_lambda.mat')
  load('lookup_cP_theta_lambda.mat'); % cP as function of TSR and pitch angle
end

if isfile('rated_values.mat')
  load('rated_values.mat');           % rated values of wind speed and omega
  % parameters for the file lookup_pitch.m
  stall_lim = -4*pi/180;                    % initial stall limit [°]
  feather_lim = 4*pi/180;                   % initial feathering limit [°]
  p_item = 50;                              % # of item to investigate
  V0_rated = rated_values(1);               % rated windspeed [m/s]
  velocity_vector = [0:0.5:30];             % range of wind speed [m/s]
  velocity_vector(end) = V0_rated;          % add the rated windspeed
  velocity_vector = sort(velocity_vector);
  velocity_item = size(velocity_vector, 2); % # of wind speed
  omega_rated = rated_values(2);            % rated wind speed [rad/s]
  lambda_opt = rated_values(4);             % optimum TSR
  cp_max = rated_values(5);                 % maximum cp
else
  disp(['Attention: rated values may not have been computed. ' ...
    'Run lookup_cp.m first']);
end

if isfile('lookup_pitch.mat') 
  load('lookup_pitch.mat');
else
   disp(['Attention: pitch angle values may not have been computed. ' ...
    'Run lookup_pitch.m first']);
end

% airfoil parameters 
filenames = [ "airfoil_data\cylinder", "airfoil_data\FFA-W3-600",  ...
  "airfoil_data\FFA-W3-480", "airfoil_data\FFA-W3-360", ...
  "airfoil_data\FFA-W3-301", "airfoil_data\FFA-W3-241"];
blade_filename = "airfoil_data\bladedat.txt";
thick_prof = [100 60 48 36 30.1 24.1]; % t/c ratio
[aoa_mat, cl_mat, cd_mat] = load_airfoil_data(filenames);
[r_vector, c_vector, beta_vector, thick_vector] = load_blade_data( ...
  blade_filename);
r_item = size(r_vector, 2); % number of cross sections along the blade
r_item_no_tip = r_item - 1; % number of cross section without the tip

% pitch parameters from DTU reference turbine pag 33
velocity_reference = [4:1:25];
pitch_reference = [2.751 1.966 0.896 0.000 0.000 0.000 0.000 0.000 ... 
  4.502 7.266 9.292 10.958 12.499 13.896 15.200 16.432 17.618  18.758 ...
  19.860 20.927 21.963 22.975 ];


