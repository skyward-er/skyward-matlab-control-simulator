function [sp, c] = acquisition_Sys(sensorData, s, c)
%{
Routine to simulate the data acquisition from the sensors, that use the
class sensors in: "skyward-matlab-control-simulator\sensors"

INPUT: 
    - sensorData: STRUCT WITH ALL THE MEASURED DATA AND TIMESTAMPS

    - s:          STRUCT WITH ALL THE SENSOR OBJECTS 

    - c:          STRUCT WITH THE ASSEMBLED TOTAL MEASUREMENT VECTORS       

OUTPUT:
    - sp          PROCESSED MEASUREMENTS
    
    - c           STRUCT WITH THE ASSEMBLED TOTAL MEASUREMENT VECTORS       
%}

% Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept
% email: alessandro.delduca@skywarder.eu
% Revision date: 18/03/2021

%% Baro Acquisition loop
        sp.barometer.measures   = zeros(1,length(sensorData.barometer.time));
        sp.h_baro               = zeros(1,length(sensorData.barometer.time));
        sp.barometer.time       = sensorData.barometer.time;
        
        for ii=1:length(sensorData.barometer.time)
                sp.barometer.measures(ii)        =      s.MS580301BA01.sens(sensorData.barometer.measures(ii)/100,...
                                                            sensorData.barometer.temperature(ii) - 273.15);  
                sp.barometer.measures(ii)        =      sp.barometer.measures(ii)*100;   
                sp.h_baro(ii)    =     -atmospalt(sp.barometer.measures(ii),'None');
        end 
        c.pn_tot(c.np_old:c.np_old + size(sp.barometer.measures,2) - 1,1)    = sp.barometer.measures(1:end);
        c.hb_tot(c.np_old:c.np_old + size(sp.barometer.measures,2) - 1,1)    = sp.h_baro(1:end);
        c.time_baro(c.np_old:c.np_old + size(sp.barometer.measures,2) - 1)   =  sp.barometer.time ;
        c.np_old = c.np_old + size(sp.barometer.measures,2);      
      
%% IMU Acquisition loop
        sp.accelerometer.measures   = zeros(length(sensorData.accelerometer.time),3);
        sp.gyro.measures            = zeros(length(sensorData.gyro.time),3);
        sp.magnetometer.measures    = zeros(length(sensorData.magnetometer.time),3);  
        sp.accelerometer.time       = sensorData.accelerometer.time;
        sp.magnetometer.time        = sensorData.magnetometer.time;
        for ii=1:length(sensorData.accelerometer.time)
                [sp.accelerometer.measures(ii,1),sp.accelerometer.measures(ii,2),sp.accelerometer.measures(ii,3)] =      ...
                                                 s.ACCEL_LSM9DS1.sens(...
                                                 sensorData.accelerometer.measures(ii,1)*1000/9.81,...
                                                 sensorData.accelerometer.measures(ii,2)*1000/9.81,...
                                                 sensorData.accelerometer.measures(ii,3)*1000/9.81,...
                                                 14.8500);  
                 [sp.gyro.measures(ii,1),sp.gyro.measures(ii,2),sp.gyro.measures(ii,3)]   =      ...
                                                 s.GYRO_LSM9DS1.sens( ...
                                                 sensorData.gyro.measures(ii,1)*1000*360/2/pi,...
                                                 sensorData.gyro.measures(ii,2)*1000*360/2/pi,...
                                                 sensorData.gyro.measures(ii,3)*1000*360/2/pi,...
                                                 14.8500);
                 [sp.magnetometer.time(ii,1),sp.magnetometer.time(ii,2),sp.magnetometer.time(ii,3)]      =      ...
                                                 s.MAGN_LSM9DS1.sens( ...
                                                 sensorData.magnetometer.measures(ii,1)*0.01,...
                                                 sensorData.magnetometer.measures(ii,2)*0.01,...
                                                 sensorData.magnetometer.measures(ii,3)*0.01,...
                                                 14.8500);   
                 sp.accelerometer.measures(ii,:) = sp.accelerometer.measures(ii,:)*9.81/1000;
                 sp.gyro.measures(ii,:)  = sp.gyro.measures(ii,:)*2*pi/360/1000;                          
        end 
        c.accel_tot(c.na_old:c.na_old + size(sp.accelerometer.measures,1) - 1,:) = sp.accelerometer.measures(1:end,:) ;
        c.gyro_tot(c.na_old:c.na_old + size(sp.gyro.measures,1) - 1,:)           = sp.gyro.measures(1:end,:) ;
        c.mag_tot(c.na_old:c.na_old + size(sp.magnetometer.measures,1) - 1,:)    = sp.magnetometer.measures(1:end,:) ;
        c.time_imu(c.na_old:c.na_old + size(sp.accelerometer.measures,1) - 1)    =  sp.accelerometer.time;
        c.na_old = c.na_old + size(sp.accelerometer.measures,1);
        

%% GPS Acquisition loop
        sp.gps.positionMeasures     = zeros(length(sensorData.gps.time),3);
        sp.gps.velocityMeasures     = zeros(length(sensorData.gps.time),3);
        sp.gps.time   = sensorData.gps.time;
        for ii=1:length(sensorData.gps.time)
            [sp.gps.positionMeasures(ii,1),sp.gps.positionMeasures(ii,2),sp.gps.positionMeasures(ii,3)]   =            ...
                                                 s.GPS_NEOM9N.sens( ...
                                                 sensorData.gps.positionMeasures(ii,1),...
                                                 sensorData.gps.positionMeasures(ii,2),...
                                               - sensorData.gps.positionMeasures(ii,3),...
                                                 14.8500);  
            [sp.gps.velocityMeasures(ii,1),sp.gps.velocityMeasures(ii,2),sp.gps.velocityMeasures(ii,3)] =           ...
                                                 s.GPS_NEOM9N.sens( ...
                                                 sensorData.gps.velocityMeasures(ii,1),...
                                                 sensorData.gps.velocityMeasures(ii,2),...
                                               - sensorData.gps.velocityMeasures(ii,3),...
                                                 14.8500);  
        end
        c.gps_tot(c.ngps_old:c.ngps_old + size(sp.gps.positionMeasures,1) - 1,:)   =  sp.gps.positionMeasures(1:end,:) ;
        c.gpsv_tot(c.ngps_old:c.ngps_old + size(sp.gps.velocityMeasures,1) - 1,:)  =  sp.gps.velocityMeasures(1:end,:) ;
        c.time_gps(c.ngps_old:c.ngps_old + size(sp.gps.velocityMeasures,1) - 1)    =  sp.gps.time;
        c.ngps_old = c.ngps_old + size(sp.gps.positionMeasures,1);
end