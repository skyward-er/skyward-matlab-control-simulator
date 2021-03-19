function contSettings = controlConfig
%{

CONTROLCONFIG - This script sets up all the parameters for the control
All the parameters are stored in the "contSetting" structure.

 %}

% Load coefficients for Cd
data                    =     load('coeffs.mat');
contSettings.coeff_Cd   =     data.coeffs;

% Load the trajectories
struct_trajectories              =	load('Trajectories');
contSettings.data_trajectories   =  struct_trajectories.trajectories_saving;

% PI controler tune parameter
contSettings.Kp_1    =   50;                                               % using Fdrag nel pid --> da migliorare (magari si può ottenere variabile controllo più smooth)
contSettings.Ki_1    =   10;                                               % using Fdrag nel pid
contSettings.Kp_2    =   50;                                               % using u nel pid --> da migliorare (magari si può ottenere variabile controllo più smooth)
contSettings.Ki_2    =   37;                                               % using u nel pid
contSettings.Kp_3    =   50;                                               % using alfa_degree nel pid --> ancora da tunare
contSettings.Ki_3    =   10;                                               % using alfa_degree nel pid

% Internal parameter of controler
contSettings.I                   =   0;
contSettings.alpha_degree_prec   =   0;
contSettings.iteration_flag      =   1;
contSettings.saturation          =   false;

contSettings.m  = 22;
contSettings.g  = 9.81;
contSettings.D  = 0.15; 
contSettings.S0 = (pi*contSettings.D^2)/4; 

contSettings.a = -9.43386/1000;
contSettings.b = 19.86779/1000;

contSettings.rate_limiter_max =    60/0.2; % datasheet: 60deg/0.13s --> increased for robustness
contSettings.rate_limiter_min =   -60/0.2;
contSettings.filter_coeff     =    0.9;

% Possible range of values for the control variable
contSettings.delta_S_available = (0.0:0.001/2:0.01)'; 
% Select the PID algorithm
contSettings.flagPID          =    1;                                   % 1: Fdrag;  2: u;  3: alfa_degree;