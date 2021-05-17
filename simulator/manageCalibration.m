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
    SensorData.accelerometer.time = zeros(1,Na);
    SensorData.gyro.time = zeros(1,Na);
    
for i =  1:Nm
    SensorData.magnetometer.measures(i,:) = (A * magField)';
end
    SensorData.magnetometer.time = zeros(1,Nm);
    
for i = 1:Ngps
    SensorData.gps.velocityMeasures(i,:)  = [0,0,0];
    SensorData.gps.positionMeasures(i,:)  = [0,0,settings.z0];
end
    SensorData.gps.time = zeros(1,Ngps);
    
    [T,~,P,~] = atmosisa(settings.z0);
for i = 1:Nbar
    SensorData.barometer.measures(i,1) = P;
    SensorData.barometer.temperature(i,1) = T;
end
    SensorData.barometer.time = zeros(1,Nbar);
    
end


function A = cdmrot(q)
A       = [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2,               2*(q(1)*q(2) + q(3)*q(4)),                 2*(q(1)*q(3) - q(2)*q(4));
                 2*(q(1)*q(2) - q(3)*q(4)),      -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,                2*(q(2)*q(3) + q(1)*q(4)) ;
                 2*(q(1)*q(3) + q(2)*q(4)),               2*(q(2)*q(3) - q(1)*q(4)),       -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2];

end