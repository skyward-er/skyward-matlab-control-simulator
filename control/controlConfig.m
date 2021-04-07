function contSettings = controlConfig

% Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@kywarder.eu
% email: alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

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

% Just for plotting the setpoints of the chosen trajectory without spikes
contSettings.starting_index = 0;

% Control time sample
contSettings.sample_time         =  0.1;

% PI controler tune parameter
contSettings.Kp_1    =   20;                                               % using Fdrag nel pid --> da migliorare (magari si può ottenere variabile controllo più smooth)
contSettings.Ki_1    =   1;                                               % using Fdrag nel pid
contSettings.Kp_2    =   50;                                               % using u nel pid --> da migliorare (magari si può ottenere variabile controllo più smooth)
contSettings.Ki_2    =   40;                                               % using u nel pid
contSettings.Kp_3    =   20;                                               % using alfa_degree nel pid --> ancora da tunare
contSettings.Ki_3    =   20;                                               % using alfa_degree nel pid

% Select the PID algorithm
contSettings.flagPID           =    1;                                     % 1: control_PID (Fdrag);  2: control_LIN (u);  3: control_Servo (alfa_degree);

% Internal parameter of controler
contSettings.I                   =   0;
contSettings.alpha_degree_prec   =   0;
contSettings.iteration_flag      =   1;
contSettings.saturation          =   false;

% Physical parameters of the rocket
contSettings.m  = 22;
contSettings.g  = 9.81;
contSettings.D  = 0.15; 
contSettings.S0 = (pi*contSettings.D^2)/4; 

% Parameters for the function get extension from angle
contSettings.a  = -9.43386/1000;                                            
contSettings.b  = 19.86779/1000;                                           

contSettings.rate_limiter      =    60/0.2;                                % datasheet: 60deg/0.13s --> increased for robustness
contSettings.filter_coeff      =    0.9;

% Possible range of values for the control variable
contSettings.delta_S_available = (0.0:0.001/2:0.01)'; 

