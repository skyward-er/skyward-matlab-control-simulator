%{

SEND DATA OVER SERIAL - use the serial to send flight phases flags and 
                        sensors and velocities data to the microcontroller

Author: Emilio Corigliano
Skyward Experimental Rocketry | ELC Dept | electronics@skywarder.eu
email: emilio.corigliano@skywarder.eu
Release date: 10/03/2021

%}

function [] = sendDataOverSerial(sensorData, flags)

sensorData.flags.flagFligth = flags(1);
sensorData.flags.flagAscent = flags(2);
sensorData.flags.flagBurning = flags(3);
sensorData.flags.flagAeroBrakes = flags(4);
sensorData.flags.flagPara1 = flags(5);
sensorData.flags.flagPara2 = flags(6);

serialbridge("Write", structToSingles(sensorData));

end