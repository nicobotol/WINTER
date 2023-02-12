%% PARAMETERS

% Add some folders to the path
addpath("aerodynamic_functions")
addpath("aerodynamic_functions\airfoil_data")
addpath("lookup")

%% Parameters for the lookup tables generation
pitch_range = deg2rad([-15 90]);              % range for picth angle [rad]
pitch_item = ceil(rad2deg(diff(pitch_range)));% # of guess pitch 
lambda_range = [0 18];                        % range for the TSR
lambda_item = diff(lambda_range)*3;           % # of guess pitch 
% distribute lambda and TSR in their ranges
lambda_vector = linspace(lambda_range(1), lambda_range(2), lambda_item); 
pitch_vector = linspace(pitch_range(1), pitch_range(2), pitch_item);

% load mesh for cP, cT, rated values, and pitch angle (computed in 
% lookup_cp.m and lookup_pitch.m)
if isfile('lookup\lookup_cP_theta_lambda.mat')
  load('lookup\lookup_cP_theta_lambda.mat'); % cP(TSR, pitch angle)
  load('lookup\lookup_cT_theta_lambda.mat'); % cT(TSR, pitch angle)
end

% parameters for the file lookup_pitch.m
if isfile('lookup\rated_values.mat')
  load('lookup\rated_values.mat');          % rated values of ws and omega
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

if isfile('lookup\lookup_pitch.mat') 
  load('lookup\lookup_pitch.mat');
else
   disp(['Attention: pitch angle values may not have been computed. ' ...
    'Run lookup_pitch.m first']);
end

%% Physical parameters
rho = 1.225;                % air density [kg/m^3]
stop_time = 150;            % max time to investigaste [s]
font_size = 25;             % fontsize for plots
line_width = 2;             % line width for plots
marker_size = 12;           % marker size for plots
i_max = 300;                % max number of accepted iterations
fake_zero = 1e-8;           % thershold for exiting the loop
a_guess = 0;                % initial guess for the BEM code
a_prime_guess = 0.1;        % initial guess for the BEM code
V0_cut_in = 4;              % cut in wind speed [m/s]
V0_cut_out = 25;            % cut out wind speed [m/s]
time_step = 1e-3;           % simulation time step [s]

% Rotor parameters
rotor.R = 89.17;            % rotor radius [m]
rotor.A = rotor.R^2*pi;     % rotor area [m^2]
rotor.blades = 3;           % number of blades [#]
rotor.V0_cutin = 4;         % cut in wind velocity [m/s]
rotor.V0_cutout = 25;       % cut out wind velocity [m/s]
rotor.P_rated = 10.64e6;    % rated power [W]
rotor.mass = 1.3016e5;      % mass [kg]
rotor.I = 1.5617e8;         % inertia wrt rotational axis [kgm^2]
% rotor.omega_r = 1.0457;     % initial rotational speed [rad/s]
rotor.omega_r = 0.2;      % initial rotational speed [rad/s]
rotor.B  = 0;               % rotational friction [kgm^2/s] (random placeholder)

% Gearbox_parameters
gearbox.ratio = 1;          % gearbox transmission ratio 

% Generator parameters
generator.P_rated = rotor.P_rated; % rated power [W]
generator.I = 4800;         % generator iniertia [kgm^2]
generator.B = 0.0;          % rotational friction [kgm^2/s] (random placeholder)
generator.vll = 4e3;        % rated line-to-line voltage [V]
generator.is = 1443.4;      % rated stator current [A]
generator.fe = 26.66;       % rated stator frequency [Hz]
generator.p = 32;           % number of poles
generator.Ld = 1.8e-3;      % d-axis stator inductance [H]
generator.Lq = 1.8e-3;      % q-axis stator inductance [H]
generator.Rs = 64e-3;       % stator resistance [ohm]
generator.Lambda = 19.49;   % magnet flux-linkage [Wb]
generator.tau_c = 835e-6;   % q-axis control time constant [s]
% generator.p_ctrl = 1e3;     % gain for the Ig reference
% generator.k_ctrl = 0.01;    % paramter for the Iq refernce
generator.iq_pm = 50;       % phase margin for the Iq controller [°]
generator.iq_omegaBP = 1e2; % Iq controller crossover frequency [rad/s]
generator.omega_pm = 60;    % phase margin for the speed controller [°]
generator.omega_omegaBP=1e3;% speed controller crossover frequency [rad/s]
generator.K_opt = ...
rho*pi*rotor.R^5*cp_max*gearbox.ratio^3/(2*lambda_opt^3); % ref. torque
                                                          % const. [kgm^2]
generator.omega_LP = 0.2;   % freq. of the II order speed LP filter [rad/s]  
generator.zeta_LP = 0.7;    % damping of the II order speed LP filter [-]

% Blade parameters
blade.Kp = 1e3;             % proportional gain
blade.Ki = 5e2;             % integrative gain
blade.mass = 4.3388e4;      % mass [kg]
blade.I = 5.2056e7;         % inertia wrt the rotor rotational axis [kgm^2]
blade.zetap = 0.7;          % damping ratio of the pitch actuator
blade.omegap = 2*pi;        % undamped nat. freq. of the actuator [rad/s]
blade.pitch_rate =10*pi/180;% maximum pitch rate [rad/s]
% Equivlent inertia and damping, referred to the rotor side of the
% transmission
I_eq = rotor.I + generator.I/gearbox.ratio^2; % equiv. inertia [kgm^2]
B_eq = rotor.B + generator.B/gearbox.ratio^2; % equiv. damping [kgm^2/s]

% Transform the struct of parameters into buses for simulink
rotor_bus_info = Simulink.Bus.createObject(rotor); 
rotor_bus = evalin('base', rotor_bus_info.busName);
generator_bus_info = Simulink.Bus.createObject(generator); 
generator_bus = evalin('base', generator_bus_info.busName);
gearbox_bus_info = Simulink.Bus.createObject(gearbox); 
gearbox_bus = evalin('base', gearbox_bus_info.busName);
blade_bus_info = Simulink.Bus.createObject(blade); 
blade_bus = evalin('base', blade_bus_info.busName);

% Wind time history
wind_speed = load('usim.dat');                    % [m/s] 
sample_time = wind_speed(2,1) - wind_speed(1, 1); % WS sample time [s]

% Pitching startegy (feathering or stall)
pitch_strategy = 0;  % 0     -> feathering
                     % 1     -> stalling
                     % else  -> no pitch control

%% Airfoil parameters 
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

%% Parameters from DTU reference turbine pag 33
% column 1 -> windspeed [m/s]
% column 2 -> pitch angle [°]
% column 3 -> rotaional speed [rpm]
% column 4 -> cP [-]
% column 5 -> cT [-]
reference = load('airfoil_data\DTU_10MW_reference.txt'); 

