%{

SEND DATA OVER SERIAL - use the serial to send flight phases flags and 
                        sensors and velocities data to the microcontroller

Author: Emilio Corigliano
Skyward Experimental Rocketry | ELC Dept | electronics@skywarder.eu
email: emilio.corigliano@skywarder.eu
Release date: 10/03/2021

%}

function [] = sendDataOverSerial(data, flags)
if(isfield(data, 'h_baro'))
    data = rmfield(data, 'h_baro');
end

if(isfield(data.barometer, 'temperature'))
    data.barometer = rmfield(data.barometer, 'temperature');
end

dataToBeSent.accelerometer = data.accelerometer;

dataToBeSent.gyro = data.gyro;

dataToBeSent.magnetometer = data.magnetometer;

dataToBeSent.gps.positionMeasures = data.gps.positionMeasures;
dataToBeSent.gps.velocityMeasures = data.gps.velocityMeasures;

dataToBeSent.barometer = data.barometer;

dataToBeSent.kalman.z = data.kalman.z;
dataToBeSent.kalman.vz = data.kalman.vz;
dataToBeSent.kalman.vMod = data.kalman.vMod;



dataToBeSent.flags.flagFligth = cast(flags(1), "double");
dataToBeSent.flags.flagAscent = cast(flags(2), "double");
dataToBeSent.flags.flagBurning = cast(flags(3), "double");
dataToBeSent.flags.flagAeroBrakes = cast(flags(4), "double");
dataToBeSent.flags.flagPara1 = cast(flags(5), "double");
dataToBeSent.flags.flagPara2 = cast(flags(6), "double");

serialbridge("Write", structToSingles(dataToBeSent));

end