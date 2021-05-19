%{

HIL CONFIG - This script sets up all the parameters for the 
hardware-in-the-loop simulation.
All the parameters are stored in the "hil_settings" structure.

Author: Luca Conterio
Skyward Experimental Rocketry | ELC Dept | electronics@skywarder.eu
email: luca.conterio@skywarder.eu
Release date: 10/03/2021

%}

%% HIL SIMULATION SETTINGS

% Windows, e.g. "COM6"
% Linux, e.g. "/dev/ttyACM0"
hil_settings.serial_port = "COM6";

% Windows, e.g. 256000
% Linux, e.g. 115200
% See serialib for supported baudrates
hil_settings.baudrate = 115200;