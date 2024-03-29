%% PARAMETERS

% Add some folders to the path
addpath("aerodynamic_functions")
addpath("aerodynamic_functions\airfoil_data")
addpath("lookup\")
addpath("run")
addpath("wind_series")
addpath("PMSM")
addpath("simulink\")
addpath("plot\")

%% Parameters for the lookup tables generation
pitch_range = deg2rad([-15 90]);              % range for picth angle [rad]
pitch_item = 2*ceil(rad2deg(diff(pitch_range)));% # of guess pitch 
lambda_range = [0 18];                        % range for the TSR (original)
lambda_item = 2*diff(lambda_range)*3;           % # of guess TSR 
% distribute lambda and TSR in their ranges
lambda_vector = linspace(lambda_range(1), lambda_range(2), lambda_item); 
pitch_vector = linspace(pitch_range(1), pitch_range(2), pitch_item);

% load mesh for cP, cT, rated values, and pitch angle (computed in 
% lookup_cp.m and lookup_pitch.m)
if exist('lookup_cP_theta_lambda.mat', 'file')
  load('lookup_cP_theta_lambda.mat'); % cP(TSR, pitch angle)
  load('lookup_cT_theta_lambda.mat'); % cT(TSR, pitch angle)
  load('lookup_static_values.mat');   % static values for rows:
                                      % 1 -> wind speed [m/s]
                                      % 2 -> rotor rotational speed [rad/s]
                                      % 3 -> gen. rotational speed [rad/s]
                                      % 4 -> rotor torque [Nm]
                                      % 5 -> generator torque [Nm]
                                      % 6 -> rotor power [W]
                                      % 7 -> generator power [W]
end

% parameters for the file lookup_pitch.m
if exist('rated_values.mat', 'file')
  load('rated_values.mat');          % rated values of ws and omega
  stall_lim = -4*pi/180;                    % initial stall limit [°]
  feather_lim = 4*pi/180;                   % initial feathering limit [°]
  p_item = 50;                              % # of item to investigate
  V0_rated = rated_values(1);               % rated windspeed [m/s]
  velocity_spacing = 0.5;                   % spcae between two WS [m/s]
  velocity_vector = [0:velocity_spacing:30];% range of wind speed [m/s]
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

if exist('lookup_pitch.mat', 'file') 
  load('lookup_pitch.mat');
else
   disp(['Attention: pitch angle values may not have been computed. ' ...
    'Run lookup_pitch.m first']);
end

%% Physical parameters
rho = 1.225;                % air density [kg/m^3]
font_size = 25;             % fontsize for plots
line_width = 2;             % line width for plots
marker_size = 12;           % marker size for plots
i_max = 300;                % max number of accepted iterations
fake_zero = 1e-8;           % thershold for exiting the loop
a_guess = 0;                % initial guess for the BEM code
a_prime_guess = 0.1;        % initial guess for the BEM code
V0_cut_in = 4;              % cut in wind speed [m/s]
V0_cut_out = 25;            % cut out wind speed [m/s]

simulation.model = 2;       % choice of the model
                            % 1 -> without power controller
                            % 2 -> with power controller
if simulation.model == 1    % without power controller
  simulation.mdl = 'winter_simulink_without_PC'; 
elseif simulation.model == 2 % with power controller
  simulation.mdl = 'winter_simulink_with_PC'; 
end
simulation.stop_time = [80]; % max time to investigaste [s]
simulation.time_step_H=1e-2;% time step for the mechanical part [s]
simulation.time_step_L=5e-5;% time step for the electrical part [s]
simulation.type = 6;        % 1 -> constant wind speed
                            % 2 -> ramp
                            % 3 -> generated wind series
                            % 4 -> generator step response
                            % 5 -> generated WS and parametrization plot
                            % 6 -> ramp and parametrization plot
simulation.plot_time = [];  % time from the end of the simulation to 
                            % average the response [s]
% simulation.plot_step = simulation.plot_time/simulation.time_step;
simulation.print_figure = 1;% enables or disable plot's autosaving 
                            % 1 -> plot enabled
                            % 0 -> plot disable
simulation.seed = 3;        % seed for the random number generation

% Rotor parameters
rotor.R = 89.17;            % rotor radius [m]
rotor.A = rotor.R^2*pi;     % rotor area [m^2]
rotor.blades = 3;           % number of blades [#]
rotor.V0_cutin = 4;         % cut in wind velocity [m/s]
rotor.V0_cutout = 25;       % cut out wind velocity [m/s]
rotor.P_rated = 10.64e6;    % rated power [W]
rotor.mass = 1.3016e5;      % mass [kg]
rotor.I = 1.5617e8;         % inertia wrt rotational axis [kgm^2]
rotor.omega_R = lambda_opt*10.5/rotor.R;  % initial rotational speed [rad/s]
rotor.B  = 1000;               % rotational friction [kgm^2/s] (random placeholder)
rotor.K_opt = rho*pi*rotor.R^5*cp_max/(2*lambda_opt^3);

% Gearbox_parameters
gearbox.ratio = 1;          % gearbox transmission ratio 

% Generator parameters
generator.P_rated = rotor.P_rated; % rated power [W]
generator.omega_rated = omega_rated/gearbox.ratio; % rated speed gen. side [rad/s]
generator.I = 4800;         % generator iniertia [kgm^2]
generator.B = 0.0;          % rotational friction [kgm^2/s] (random placeholder)
generator.vll = 4e3;        % rated line-to-line voltage [V]
generator.is = 1443.4;      % rated stator current [A]
generator.fe = 26.66;       % rated stator frequency [Hz]
generator.p = 320;          % number of poles
generator.Ld = 1.8e-3;      % d-axis stator inductance [H]
generator.Lq = 1.8e-3;      % q-axis stator inductance [H]
generator.Rs = 64e-3;       % stator resistance [ohm]
generator.Lambda = 19.49;   % magnet flux-linkage [Wb]
generator.tau_c = 500e-6;   % q-axis control time constant [s]
% generator.p_ctrl = 1e3;   % gain for the Ig reference
% generator.k_ctrl = 0.01;    % paramter for the Iq refernce
generator.iq_pm = 70;       % phase margin for the Iq controller [°]
generator.iq_omegaBP = 1.5e3; % Iq controller crossover freq. [rad/s]
generator.omega_pm = 70;  % phase margin for the speed controller [°]
generator.omega_omegaBP=generator.iq_omegaBP/10;% speed controller crossover frequency [rad/s]
generator.K_opt = ...
rho*pi*rotor.R^5*cp_max*gearbox.ratio^3/(2*lambda_opt^3); % ref. torque
                                                          % const. [kgm^2]
generator.design = 0;       % 0 enables manual design of the controller
                            % 1 enables pidtune design of the controller
generator.design_omega = 0; % 0 enables manual design of speed controller
                            % 1 enables pidtune design of the controller                          
generator.bode_plot = 0;    % 1 enables bode plot, 0 disables it
generator.alpha_omega= 2.51;% Speed low pass filter frequency [rad/s]  
generator.power_ctrl_kp=0.5;% power controller gain 0.5
generator.power_ctrl_ki=5.5;% power controller gain
% generator.kp = 7.33e7;
% generator.kd = 0;
% generator.ki = 1.32e7;
generator.omega1_min = 0;
generator.omega2_min = 1.05*generator.omega1_min;
generator.omega1_max = 0.90*generator.omega_rated;
generator.omega2_max = 0.95*generator.omega_rated;
generator.torque_full = generator.K_opt*generator.omega_rated^2; % full load torque [Nm]

% Blade parameters
blade.mass = 4.3388e4;      % mass [kg]
blade.I = 5.2056e7;         % inertia wrt the rotor rotational axis [kgm^2]
blade.zetap = 0.7;          % damping ratio of the pitch actuator
blade.omegap = 2*pi;        % undamped nat. freq. of the actuator [rad/s]
blade.pitch_rate=10*pi/180; % maximum pitch rate [rad/s]
blade.alpha_beta = 2*pi*0.4;% constant for the pitch error filter [rad/s]
blade.kp_schedule = [-59.871 46.281 -7.814 -2.541 1];
blade.ki_schedule = [27.689 -31.926 13.128 -2.405 0.351]*3;
% blade.kp_schedule = 0.4;
% blade.ki_schedule = 0.2;
% blade.kp_tab = [-2, 0,4,6,8,10.5,12,13,14,16,17,18,19,20,21,22,23,24,...
%                       25,26,27;
%                    2.35,2.35,1.5,1.25,1.08,0.98,0.9,0.82,0.78,0.71,0.68,...
%                       0.65,0.6,0.57,0.55,0.5,0.49,0.43,0.43,0.42,0.41];
% blade.ki_tab = [ -2, 0,4,6,8,10.5,12,13,14,16,17,18,19,20,21,22,23,24,...
%                       25,26,27;
%                     0.92, 0.92,0.58,0.49,0.42,0.39,0.36,0.32,0.305,0.29,...
%                       0.28,0.26,0.23,0.22,0.21,0.205,0.195,0.19,0.185,...
%                       0.18,0.17];
blade.theta_f = 0.5*pi/180; % added term in the generator control [rad]
blade.initial_pitch = 0*pi/180;  % initial condition for pitch ctrl 
blade.K1 = 164.13; % Linear coeff. in aero gain scheduling [deg]
blade.K2 = 702.09; % Quadratic coeff. in aero gain scheduling [deg^2]
blade.omega2omega0ratio = 1.3; % Relative speed for double nonlinear gain
blade.pitch_min = 0;        % minimum pitch angle [rad]
blade.kp = 0.592;
blade.kpp = 4e-9;
blade.ki = 0.133;
blade.kip = 4e-9;
blade.kd = 0;


% Wind parameters
wind.mean = [15 10];                % 10 minutes mean wind speed [m/s]]
wind.turbulence = 0.1*wind.mean; % 10 min std (i.e. turbulence) [m/s]
wind.height = 119.0;            % height where to measure the wind [m]
wind.sample_f = 50;             % wind sample frequncy [Hz]
wind.sample_t = 1/wind.sample_f;% wind sample time [s]
wind.ramp_WS_start = 10.5;        % wind speed at the start of the ramp [m/s]
wind.ramp_WS_stop = 14;         % wind speed at the stop of the ramp [m/s]
wind.ramp_time_start = [1]; % time speed at the start of the ramp [s]
wind.ramp_time_stop = [simulation.stop_time];  % time speed at the stop of the ramp [s]

switch simulation.type
  case {1, 3, 5}
    wind.WS_len = length(wind.mean);  % number of separated WSs to test
  case {2, 6}
    wind.WS_len = length(wind.ramp_time_start);
  case 4
    wind.WS_len = 1; 
end

% Struct where to save the simulations results
out_store = cell(1, wind.WS_len);

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

% Pitching startegy (feathering or stall)
pitch_strategy = 0;  % 0     -> feathering
                     % 1     -> stalling
                     % else  -> no pitch control

%% Airfoil parameters 
filenames = [ "airfoil_data\cyli1", "airfoil_data\cyli1", "airfoil_data\DU40_A17",... 
  "airfoil_data\DU35_A17", "airfoil_data\DU30_A17", ...
  "airfoil_data\DU25_A17", "airfoil_data\DU21_A17","airfoil_data\NA64_A17"];
blade_filename = "airfoil_data\bladedat.txt";
thick_prof = [100 100 40.5 35.09 30 25 21 18 ]; % t/c ratio
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

%% Plotting options
% colors
colors_vect = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; ...
               [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; ...
               [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; ...
               [0.6350 0.0780 0.1840]];
% pathe where to save the images
path_images = ['C:\Users\Niccolò\Documents\UNIVERSITA\TESI_MAGISTRALE\' ...
  '\report\images'];

% parameters for the generator step response plot
step_t_start = 0.4;
step_t_stop = 0.65;
step_y_min = 0.9;
step_y_max = 1.06;

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  font_size)
set(0,'DefaultLegendFontSize', font_size)