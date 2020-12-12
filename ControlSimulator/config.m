%{

CONFIG - This script sets up all the parameters for the simulation
All the parameters are stored in the "settings" structure.

Author: Francesco Colombi
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: francesco.colombi@skywarder.eu
Release date: 16/04/2016

%}

%% SIMULATION SETTINGS
settings.electronics = false;
settings.ascentOnly = false;
settings.ballisticFligth = true;
settings.dataNoise = false;

%% LAUNCH SETUP
% launchpad
settings.z0 = 109;                                                          %[m] Launchpad Altitude
settings.lrampa = 4.9;                                                      %[m] LaunchPad route (distance from ground of the first hook)
settings.lat0 = 44.519272;                                                  % Launchpad latitude
settings.lon0 = 11.642333;                                                  % Launchpad longitude
settings.launchDate = [2021, 10, 15];                                       % [YYYY, mm, dd] date of launch

% launchpad directions
% for a single run the maximum and the minimum value of the following
% angles must be the same.
settings.OMEGA = 84*pi/180;                                                 %[rad] Minimum Elevation Angle, user input in degrees (ex. 80)
settings.PHI = 0*pi/180;                                                  %[rad] Minimum Azimuth Angle from North Direction, user input in degrees (ex. 90)

%% ENGINE DETAILS
% load motors data 
DATA_PATH = '../data/';
filename_full = strcat(DATA_PATH,'MotorsList.mat');
motors = load(filename_full,'MotorsByName');
motors = motors.MotorsByName;

name = 'M2020'; % USE THIS
%name = 'M1890';  % UNFEASIBLE, TOO STRONG TO DECELERATE
%name = 'M1800'; % UNFEASIBLE, TOO STRONG TO DECELERATE
settings.motor.exp_time = motors.(name).t;
settings.motor.exp_thrust = motors.(name).T;
settings.mp = motors.(name).mp;                                             % [kg]   Propellant Mass                                                
settings.tb = motors.(name).t(end) ;                                        % [s]    Burning time
settings.mfr = settings.mp/settings.tb;                                     % [kg/s] Mass Flow Rate
settings.ms = 21.05;                                                        % [kg]   Total Mass
settings.m0 = settings.ms + settings.mp;                                    % [kg]   Structural Mass
settings.mnc = 0.400;                                                       % [kg]   Nosecone Mass

clear ('motors','name')

%% GEOMETRY DETAILS
% This parameters should be the same parameters set up in MISSILE DATCOM
% simulation.

settings.C = 0.15;                                                          % [m]      Caliber (Fuselage Diameter)
settings.S = pi*settings.C^2/4;                                             % [m^2]    Cross-sectional Surface
L = 3;                                                                      % [m]      Rocket length

%% MASS GEOMERTY DETAILS
% x-axis: along the fuselage
% y-axis: right wing
% z-axis: downward

% inertias for full configuration (with all the propellant embarqued) obtained with CAD's
settings.Ixxf = 0.0540;                     % [kg*m^2] Inertia to x-axis
settings.Iyyf = 13.7274;                    % [kg*m^2] Inertia to y-axis
settings.Izzf = 13.7302;                    % [kg*m^2] Inertia to z-axis

% inertias for empty configuration (all the propellant consumed) obtained with CAD's
settings.Ixxe = 0.0498;                     % [kg*m^2] Inertia to x-axis
settings.Iyye = 11.5612;                    % [kg*m^2] Inertia to y-axis
settings.Izze = 11.5640;                    % [kg*m^2] Inertia to z-axis

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
filename_full = strcat(DATA_PATH,'full.mat');
CoeffsF = load(filename_full,'Coeffs');
settings.CoeffsF = CoeffsF.Coeffs;
clear('CoeffsF');

% Coefficients in empty configuration
filename_empty = strcat(DATA_PATH,'empty.mat');
CoeffsE = load(filename_empty,'Coeffs');
settings.CoeffsE = CoeffsE.Coeffs;
clear('CoeffsE');

s = load(filename_empty,'State');
settings.Alphas = s.State.Alphas';
settings.Betas = s.State.Betas';
settings.Altitudes = s.State.Altitudes';
settings.Machs = s.State.Machs';
settings.Controls = s.State.hprot';
clear('s');

%% CONTROL AND SENSOR FREQUENCIES
settings.frequencies.controlFrequency = 10;                                 % [hz] control action frequency 
settings.frequencies.accelerometerFrequency = 100;                          % [hz] control action frequency 
settings.frequencies.gyroFrequency = 100;                                   % [hz] control action frequency 
settings.frequencies.magnetometerFrequency = 100;                           % [hz] control action frequency 
settings.frequencies.gpsFrequency = 10;                                     % [hz] control action frequency 
settings.frequencies.barometerFrequency = 20;                               % [hz] control action frequency 

%% CONTROL SETTINGS 

settings.Mach_control = 0.7;                                                % Mach of activation of aerobrakes 
settings.brakesWidth = 0.088;                                               % [m] aerobrakes width (the fixed in-plane length)
settings.brakesHeigth = 0.0387;                                               % [m] max aerobrakes heigth (the control action)
settings.Atot = settings.brakesWidth*settings.brakesHeigth*3;               % [m^2] total area of aerobrakes (100% out)

%% PARACHUTES DETAILS
% parachute 1
settings.para(1).S = 1.55;                                                  % [m^2]   Surface
settings.para(1).mass = 0.4;                                                % [kg]   Parachute Mass
settings.para(1).CD = 0.75;                                                 % [/] Parachute Drag Coefficient
settings.para(1).CL = 0;                                                    % [/] Parachute Lift Coefficient
settings.para(1).z_cut = 300;                                               % [m] Final altitude of the parachute

% parachute 2
settings.para(2).S = 10.5;                                                  % [m^2]   Surface
settings.para(2).mass = 0.8;                                                % [kg]   Parachute Mass
settings.para(2).CD = 0.7;                                                  % [/] Parachute Drag Coefficient
settings.para(2).CL = 0;                                                    % [/] Parachute Lift Coefficient
settings.para(2).z_cut = 0;                                                 % [m] Final altitude of the parachute
                               
%% INTEGRATION OPTIONS
settings.ode.final_time =  2000;                                            % [s] Final integration time

% create an option structure for the integrations:

% - AbsTol is the threshold below which the value of the solution becomes unimportant
% - RelTol is the tolerance betweeen two consecutive values
% - Events is the event function that defines when the integration must be
% - stopped (it has to be created)
% - InitialStep is the highest value tried by the solver

settings.ode.optionsasc1 = odeset('Events',@event_mach,'InitialStep',1);    % ODE options for not control phase 

%% WIND DETAILS
% select which model you want to use:

%%%%% Input wind
settings.wind.input = false;
% Wind is generated for every altitude interpolating with the coefficient defined below

settings.wind.input_ground = 7;                                             % wind magnitude at the ground [m/s]
settings.wind.input_alt = [0 100 600 750 900 1500 2500];                    % altitude vector [m]
settings.wind.input_mult = [0 0 10 15 20 30 30];                            % percentage of increasing magnitude at each altitude
settings.wind.input_azimut = [30 30 30 30 30 30 30];                        % wind azimut angle at each altitude (toward wind incoming direction) [deg]

settings.wind.input_uncertainty = [1, 1];
% settings.wind.input_uncertainty = [a,b];      wind uncertanties:
% - a, wind magnitude percentage uncertanty: magn = magn *(1 +- a)
% - b, wind direction band uncertanty: dir = dir 1 +- b

%%%%% Random wind model
% if both the model above are false

% Wind is generated randomly from the minimum to the maximum parameters which defines the wind.
% Setting the same values for min and max will fix the parameters of the wind.
settings.wind.MagMin = 1;                                                   % [m/s] Minimum Magnitude
settings.wind.MagMax = 1;                                                   % [m/s] Maximum Magnitude
settings.wind.ElMin = 0*pi/180;                                             % [rad] Minimum Elevation, user input in degrees (ex. 0)
settings.wind.ElMax = 0*pi/180;                                             % [rad] Maximum Elevation, user input in degrees (ex. 0) (Max == 90 Deg)
settings.wind.AzMin = (360)*pi/180;                                         % [rad] Minimum Azimuth, user input in degrees (ex. 90)
settings.wind.AzMax = (360)*pi/180;                                         % [rad] Maximum Azimuth, user input in degrees (ex. 90)

% NOTE: wind azimuth angle indications (wind directed towards):
% 0 deg (use 360 instead of 0)  -> North
% 90 deg                        -> East
% 180 deg                       -> South
% 270 deg                       -> West

%% PLOT DETAILS
settings.plots = true;
