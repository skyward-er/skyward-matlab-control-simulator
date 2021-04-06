%{

CONFIG - This script sets up all the parameters for the simulation
All the parameters are stored in the "settings" structure.

Author: Francesco Colombi
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: francesco.colombi@skywarder.eu
Release date: 16/04/2016

%}

%% SIMULATION SETTINGS
settings.electronics        =   false;                                     % Switch on when testing with Hardware in the loop HIL
settings.ascentOnly         =   true;                                      % Switch on to simulate only the ascent phase untill the apogee
settings.ballisticFligth    =   true;                                     % Switch on to simulate the balistic fligth without any parachute
settings.control            =   true;                                      % Switch on to simulate the control
settings.dataNoise          =   true;                                      % Switch on to simulate the data acquisiton from sensors
settings.launchWindow       =   false;                                     % Switch off this to avoid pausing the launch till you press the launch button
settings.Kalman             =   true;                                      % Switch on to run the kalman algorithm
settings.Ada                =   true;                                      % Switch on to run the apogee detection algorithm

%% LAUNCH SETUP
% launchpad for Pont De Sor
settings.z0            =   109;                                            %[m] Launchpad Altitude
lpin                   =   1.150;                                          %[m] Distance from base of second pin
settings.lrampa        =   5.9 - lpin;                                     %[m] LaunchPad route (total available route)
settings.lat0          =   39.201778;                                      %[°] Launchpad latitude
settings.lon0          =  -8.138368;                                       %[°] Launchpad longitude
settings.launchDate    =  [2021, 10, 15];                                  %[YYYY, mm, dd] date of launch

settings.g0 = gravitywgs84(settings.z0, settings.lat0);                    % Gravity costant at launch latitude and altitude

% launchpad directions
settings.OMEGA         =   84*pi/180;                                      %[rad] Minimum Elevation Angle, user input in degrees (ex. 80)
settings.PHI           =   0*pi/180;                                       %[rad] Minimum Azimuth Angle from North Direction, user input in degrees (ex. 90)


%% ENGINE DETAILS
% load motors data 
DATA_PATH   =    '../data/';
filename    =   strcat(DATA_PATH,'Motors.mat');
Motors      =   load(filename);
Motors      =   [Motors.Cesaroni, Motors.Aerotech];


% name = 'M2000Rbis';
name = 'M2020';
% name = 'M1890';
% name = 'M1800';

n_name = [Motors.MotorName] == name;
settings.motor.exp_time     =    Motors(n_name).t;
settings.motor.exp_thrust   =    0.97*Motors(n_name).T;
settings.motor.exp_m        =    Motors(n_name).m;
settings.mp                 =    Motors(n_name).mp;                        % [kg]   Propellant Mass                                                
settings.tb                 =    Motors(n_name).t(end) ;                   % [s]    Burning time
mm                          =    Motors(n_name).mm;                        % [kg]   Total Mass of the Motor 
settings.ms                 =    17.873 + mm - settings.mp;                % [kg]   Structural Mass
settings.m0                 =    settings.ms + settings.mp;                % [kg]   Total Mass
settings.mnc                =    0.400;                                    % [kg]   Nosecone Mass

clear ('Motors','name')


%% GEOMETRY DETAILS
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.15;                                                         % [m]      Caliber (Fuselage Diameter)
settings.S = pi*settings.C^2/4;                                            % [m^2]    Cross-sectional Surface

%% MASS GEOMERTY DETAILS
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% inertias for full configuration (with all the propellant embarqued) obtained with CAD's
settings.Ixxf = 0.08;                                                      % [kg*m^2] Inertia to x-axis
settings.Iyyf = 13.05;                                                     % [kg*m^2] Inertia to y-axis
settings.Izzf = 13.05;                                                     % [kg*m^2] Inertia to z-axis

% inertias for empty configuration (all the propellant consumed) obtained with CAD's
settings.Ixxe = 0.07;                                                      % [kg*m^2] Inertia to x-axis
settings.Iyye = 10.06;                                                     % [kg*m^2] Inertia to y-axis
settings.Izze = 10.06;                                                     % [kg*m^2] Inertia to z-axis


%% AERODYNAMICS DETAILS
% These coefficients are obtained using MISSILE DATCOM
% (after parsing with the tool datcom_parser.py)
% The files are stored in the ../data folder with the following rule:
% rocket_name_full.mat | rocket_name_empty.mat
% e.g. R2a_full.mat    | R2a_empty.mat
% Relative Path of the data files (default: ../data/). Remember the trailing slash!!

% Coeffs is a 4D matrix given by Datcom that contains the aerodynamics
% coefficient computed for the input parameters (AoA,Betas,Altitudes,Machs)
% Note: All the parameters (AoA,Betas,Altitudes,Machs) must be the same for
% empty and full configuration

% DATA_PATH = '../data/';

% Coefficients in full configuration
filename_full       =       strcat(DATA_PATH,'full.mat');
CoeffsF             =       load(filename_full,'Coeffs');
settings.CoeffsF    =       CoeffsF.Coeffs;
clear('CoeffsF');

% Coefficients in empty configuration
filename_empty      =       strcat(DATA_PATH,'empty.mat');
CoeffsE             =       load(filename_empty,'Coeffs');
settings.CoeffsE    =       CoeffsE.Coeffs;
clear('CoeffsE');

s = load(filename_empty,'State');
settings.Alphas     =       s.State.Alphas';
settings.Betas      =       s.State.Betas';
settings.Altitudes  =       s.State.Altitudes';
settings.Machs      =       s.State.Machs';
settings.Controls   =       s.State.hprot';
clear('s');

%% CONTROL AND SENSOR FREQUENCIES
settings.frequencies.controlFrequency           =   10;                    % [hz] control action frequency 
settings.frequencies.accelerometerFrequency     =   100;                   % [hz] control action frequency 
settings.frequencies.gyroFrequency              =   100;                   % [hz] control action frequency 
settings.frequencies.magnetometerFrequency      =   100;                   % [hz] control action frequency 
settings.frequencies.gpsFrequency               =   10;                    % [hz] control action frequency 
settings.frequencies.barometerFrequency         =   20;                    % [hz] control action frequency 

%% KALMAN TUNING PARAMETERS
settings.kalman.dt_k          =   0.01;                                    % [s]        kalman time step
settings.kalman.sigma_baro    =   10;                                       % [m/2]   estimated barometer variance    
settings.kalman.sigma_mag     =   0.1;                                     % [mgauss^2] estimated magnetometer variance    
settings.kalman.sigma_GPS     =   5;                                       % [mg^2]     estimated GPS variance
settings.kalman.sigma_w       =   10*(1000*pi/180)^2;                      % [mdps^2]   estimated gyroscope variance;
settings.kalman.sigma_beta    =   1e-2;                                    % [mdps^2]   estimated gyroscope bias variance;

settings.kalman.v_thr         =   2.5;                                     % Velocity threshold for the detected apogee
settings.kalman.count_thr     =   5;                                       % If the apogee is detected count_thr time, the algorithm will return the apogee event
settings.kalman.counter       =   0;

settings.kalman.t_kalman      =   -1;                                      % Apogee detection timestamp
settings.kalman.flag_apo      =   false;                                   % True when the apogee is detected

% Process noise covariance matrix for the linear dynamics
settings.kalman.QLinear       =   1*...
                                 [1       0       0        0        0       0;
                                  0       1       0        0        0       0;
                                  0       0       1        0        0       0;
                                  0       0       0        0.1      0       0;
                                  0       0       0        0        0.1     0;
                                  0       0       0        0        0       0.1];

% Process noise covariance matrix for the quaternion dynamics
settings.kalman.Qq              =   [(settings.kalman.sigma_w^2*settings.kalman.dt_k+(1/3)*settings.kalman.sigma_beta^2*settings.kalman.dt_k^3)*eye(3)          0.5*settings.kalman.sigma_beta^2*settings.kalman.dt_k^2*eye(3);
                                      0.5*settings.kalman.sigma_beta^2*settings.kalman.dt_k^2*eye(3)                                              settings.kalman.sigma_beta^2*settings.kalman.dt_k*eye(3)];
%% ADA TUNING PARAMETER

settings.ada.Q           =   [1     0       0;                             % Process noise covariance matrix
                              0     1       0;
                              0     0       1;];
settings.ada.R           =   1;                                            % Measurement noise covariance matrix
settings.ada.P0          =   [  10    0      0;                            % Initial condition fo the 
                                0      10     0;                            % state covariance matrix 
                                0      0      10;];
[settings.ada.temp_ref, ~,...
 settings.ada.p_ref, ~]  =   atmosisa(0);                                  % Reference temperature in kelvin and pressure in Pa 

settings.ada.v0          =   0;                                            % Vertical velocity initial condition
settings.ada.a0          =   -100;                                         % Acceleration velocity initial condition
settings.ada.x0          =  [settings.ada.p_ref, settings.ada.v0, settings.ada.a0];         
                                                                           % Ada initial condition

settings.ada.v_thr       =   2.5;                                          % Velocity threshold for the detected apogee
settings.ada.count_thr   =   5;                                            % If the apogee is detected count_thr time, the algorithm will return the apogee event
settings.ada.counter     =   0;

settings.ada.t_ada       =   -1;                                           % Apogee detection timestamp
settings.ada.flag_apo    =   false;                                        % True when the apogee is detected


%% CONTROL SETTINGS 

settings.Mach_control    =   0.7;                                          % Mach of activation of aerobrakes 
settings.brakesWidth     =   0.088;                                        % [m] aerobrakes width (the fixed in-plane length)
settings.brakesHeigth    =   0.0387;                                       % [m] max aerobrakes heigth (the control action)
settings.Atot            =   settings.brakesWidth*settings.brakesHeigth*3; % [m^2] total area of aerobrakes (100% out)

%% PARACHUTES DETAILS
% parachute 1
settings.para(1).S      =    0.4;                                          % [m^2]   Surface
settings.para(1).mass   =    0.2;                                          % [kg]   Parachute Mass
settings.para(1).CD     =    0.78;                                         % [/] Parachute Drag Coefficient
settings.para(1).CL     =    0;                                            % [/] Parachute Lift Coefficient
settings.para(1).delay  =    1;                                            % [s] drogue opening delay
settings.para(1).z_cut  =    450;                                          % [m] Final altitude of the parachute

% parachute 2
settings.para(2).S      =    10.5;                                         % [m^2]   Surface
settings.para(2).mass   =    1.5;                                          % [kg]   Parachute Mass
settings.para(2).CD     =    0.7;                                          % [/] Parachute Drag Coefficient
settings.para(2).CL     =    0;                                            % [/] Parachute Lift Coefficient
settings.para(2).z_cut  =    0;                                            % [m] Final altitude of the parachute
                               
%% INTEGRATION OPTIONS
settings.ode.final_time =  2000;                                           % [s] Final integration time

% create an option structure for the integrations:

% - AbsTol is the threshold below which the value of the solution becomes unimportant
% - RelTol is the tolerance betweeen two consecutive values
% - Events is the event function that defines when the integration must be
% - stopped (it has to be created)
% - InitialStep is the highest value tried by the solver

settings.ode.optionsasc1 = odeset('Events',@event_mach,'InitialStep',1);   % ODE options for not control phase 

%% WIND DETAILS
% select which model you want to use:

%%%%% Input wind
settings.wind.input = true;
% Wind is generated for every altitude interpolating with the coefficient defined below

settings.wind.input_ground  = 7;                                           % wind magnitude at the ground [m/s]
settings.wind.input_alt     = [0 100 600 750 900 1500 2500 3000 3500];     % altitude vector [m]
settings.wind.input_mult    = [0 0 5 5 10 10 15 15 15];                  % percentage of increasing magnitude at each altitude
settings.wind.input_azimut  = [0 0 0 0 0 0 0 0 0];                         % wind azimut angle at each altitude (toward wind incoming direction) [deg]

settings.wind.input_uncertainty = [1, 1];
% settings.wind.input_uncertainty = [a,b];      wind uncertanties:
% - a, wind magnitude percentage uncertanty: magn = magn *(1 +- a)
% - b, wind direction band uncertanty: dir = dir 1 +- b

%%%%% Random wind model
% if both the model above are false

% Wind is generated randomly from the minimum to the maximum parameters which defines the wind.
% Setting the same values for min and max will fix the parameters of the wind.
settings.wind.MagMin    =   1;                                             % [m/s] Minimum Magnitude
settings.wind.MagMax    =   1;                                             % [m/s] Maximum Magnitude
settings.wind.ElMin     =   0*pi/180;                                      % [rad] Minimum Elevation, user input in degrees (ex. 0)
settings.wind.ElMax     =   0*pi/180;                                      % [rad] Maximum Elevation, user input in degrees (ex. 0) (Max == 90 Deg)
settings.wind.AzMin     =  (360)*pi/180;                                   % [rad] Minimum Azimuth, user input in degrees (ex. 90)
settings.wind.AzMax     =  (360)*pi/180;                                   % [rad] Maximum Azimuth, user input in degrees (ex. 90)

% NOTE: wind azimuth angle indications (wind directed towards):
% 0 deg (use 360 instead of 0)  -> North
% 90 deg                        -> East
% 180 deg                       -> South
% 270 deg                       -> West

%% PLOT DETAILS
settings.plots   =   true;