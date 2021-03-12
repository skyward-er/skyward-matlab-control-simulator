
%{

READ CONTROL OUTPUT FROM SERIAL - use the serial to get the control
                                  output from the microcontroller

Author: Emilio Corigliano
Skyward Experimental Rocketry | ELC Dept | electronics@skywarder.eu
email: emilio.corigliano@skywarder.eu
Release date: 10/03/2021

%}

function [alpha_degree] = readControlOutputFromSerial()

alpha_degree = serialbridge("Read", 1);

end

