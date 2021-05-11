function SensorData = manageCalibration(magField, settings)
freq = settings.frequencies;
Na   =  freq.accelerometerFrequency/freq.controlFrequency;
Nm 	 =  freq.magnetometerFrequency/freq.controlFrequency;
Ngps =  freq.gpsFrequency/freq.controlFrequency;
Nbar =  freq.barometerFrequency/freq.controlFrequency;

ome = settings.OMEGA;
phi = settings.PHI;
q   = eul2quat([phi,ome,0],'ZYX');
q   = [q(2:4), q(1)];
A   = cdmrot(q);
for i =  1:Na
    SensorData.accelerometer.measures(i,:) = (A*[0; 0; +9.81])';
    SensorData.gyro.measures(i,:) = [0, 0, 0];
end
    SensorData.accelerometer.time = zeros(Na,1);
    SensorData.gyro.time = zeros(Na,1);
    
for i =  1:Nm
    SensorData.magnetometer.measures(i,:) = (A * magField)';
end
    SensorData.magnetometer.time = zeros(Nm,1);
    
for i = 1:Ngps
    SensorData.gps.velocityMeasures(i,:)  = [0,0,0];
    SensorData.gps.positionMeasures(i,:)  = [0,0,settings.z0];
end
    SensorData.gps.time = zeros(Ngps,1);
    
    [T,~,P,~] = atmosisa(settings.z0);
for i = 1:Nbar
    SensorData.barometer.measures(i) = P;
    SensorData.barometer.temperature(i) = T;
end
    SensorData.barometer.time = zeros(Nbar,1);
    
end


function A = cdmrot(q)
A       = [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2,               2*(q(1)*q(2) + q(3)*q(4)),                 2*(q(1)*q(3) - q(2)*q(4));
                 2*(q(1)*q(2) - q(3)*q(4)),      -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,                2*(q(2)*q(3) + q(1)*q(4)) ;
                 2*(q(1)*q(3) + q(2)*q(4)),               2*(q(2)*q(3) - q(1)*q(4)),       -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2];

end