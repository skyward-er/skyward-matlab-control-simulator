%{

SEND DATA OVER SERIAL - use the serial to send flight phases flags and 
                        sensors and velocities data to the microcontroller

Author: Emilio Corigliano
Skyward Experimental Rocketry | ELC Dept | electronics@skywarder.eu
email: emilio.corigliano@skywarder.eu
Release date: 10/03/2021

%}

function [] = sendDataOverSerial(data, flags)

data.barometer = rmfield(data.barometer, 'temperature');

data.flags.flagFligth = cast(flags(1), "double");
data.flags.flagAscent = cast(flags(2), "double");
data.flags.flagBurning = cast(flags(3), "double");
data.flags.flagAeroBrakes = cast(flags(4), "double");
data.flags.flagPara1 = cast(flags(5), "double");
data.flags.flagPara2 = cast(flags(6), "double");

serialbridge("Write", structToSingles(data));

end