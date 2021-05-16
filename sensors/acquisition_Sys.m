function [sensorData, tot] = acquisition_Sys(sensorData, s, tot)
%{
Routine to simulate the data acquisition from the sensors, that use the
class sensors in: "skyward-matlab-control-simulator\sensors"

INPUT: 
    - sensorData: STRUCT WITH ALL THE MEASURED DATA AND TIMESTAMPS

    - s:          STRUCT WITH ALL THE SENSOR OBJECTS 

    - c:          STRUCT WITH THE ASSEMBLED TOTAL MEASUREMENT VECTORS       

OUTPUT:
    - sensorData          PROCESSED MEASUREMENTS
    
    - c           STRUCT WITH THE ASSEMBLED TOTAL MEASUREMENT VECTORS       
%}

% Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept
% email: alessandro.delduca@skywarder.eu
% Revision date: 18/03/2021

%% Baro Acquisition loop

        for ii=1:length(sensorData.barometer.time)
                sensorData.barometer.measures(ii,:)        =      s.MS580301BA01.sens(sensorData.barometer.measures(ii)/100,...
                                                                  sensorData.barometer.temperature(ii) - 273.15);  
                sensorData.barometer.measures(ii,:)        =      sensorData.barometer.measures(ii)*100;   
                sensorData.h_baro(ii,:)                    =     -atmospalt(sensorData.barometer.measures(ii),'None');
        end 
        tot.pn_tot(tot.np_old:tot.np_old + size(sensorData.barometer.measures,1) - 1,1)    = sensorData.barometer.measures(1:end);
        tot.hb_tot(tot.np_old:tot.np_old + size(sensorData.barometer.measures,1) - 1,1)    = sensorData.h_baro(1:end);
        tot.time_baro(tot.np_old:tot.np_old + size(sensorData.barometer.measures,1) - 1)   = sensorData.barometer.time ;
        tot.np_old = tot.np_old + size(sensorData.barometer.measures,2);      
      
%% IMU Acquisition loop

        for ii=1:length(sensorData.accelerometer.time)
                [sensorData.accelerometer.measures(ii,1),sensorData.accelerometer.measures(ii,2),sensorData.accelerometer.measures(ii,3)] =      ...
                                                 s.ACCEL_LSM9DS1.sens(...
                                                 sensorData.accelerometer.measures(ii,1)*1000/9.81,...
                                                 sensorData.accelerometer.measures(ii,2)*1000/9.81,...
                                                 sensorData.accelerometer.measures(ii,3)*1000/9.81,...
                                                 14.8500);  
                 [sensorData.gyro.measures(ii,1),sensorData.gyro.measures(ii,2),sensorData.gyro.measures(ii,3)]   =      ...
                                                 s.GYRO_LSM9DS1.sens( ...
                                                 sensorData.gyro.measures(ii,1)*1000*360/2/pi,...
                                                 sensorData.gyro.measures(ii,2)*1000*360/2/pi,...
                                                 sensorData.gyro.measures(ii,3)*1000*360/2/pi,...
                                                 14.8500);
                 [sensorData.magnetometer.measures(ii,1),sensorData.magnetometer.measures(ii,2),sensorData.magnetometer.measures(ii,3)]      =      ...
                                                 s.MAGN_LSM9DS1.sens( ...
                                                 sensorData.magnetometer.measures(ii,1)*0.01,...
                                                 sensorData.magnetometer.measures(ii,2)*0.01,...
                                                 sensorData.magnetometer.measures(ii,3)*0.01,...
                                                 14.8500);   
                 sensorData.accelerometer.measures(ii,:) = sensorData.accelerometer.measures(ii,:)*9.81/1000;
                 sensorData.gyro.measures(ii,:)  = sensorData.gyro.measures(ii,:)*2*pi/360/1000;
                 sensorData.magnetometer.measures(ii,:) =  sensorData.magnetometer.measures(ii,:)/0.01;
        end 
        tot.accel_tot(tot.na_old:tot.na_old + size(sensorData.accelerometer.measures,1) - 1,:) = sensorData.accelerometer.measures(1:end,:) ;
        tot.gyro_tot(tot.na_old:tot.na_old + size(sensorData.gyro.measures,1) - 1,:)           = sensorData.gyro.measures(1:end,:) ;
        tot.mag_tot(tot.na_old:tot.na_old + size(sensorData.magnetometer.measures,1) - 1,:)    = sensorData.magnetometer.measures(1:end,:) ;
        tot.time_imu(tot.na_old:tot.na_old + size(sensorData.accelerometer.measures,1) - 1)    =  sensorData.accelerometer.time;
        tot.na_old = tot.na_old + size(sensorData.accelerometer.measures,1);
        

%% GPS Acquisition loop

        for ii=1:length(sensorData.gps.time)
            [sensorData.gps.positionMeasures(ii,1),sensorData.gps.positionMeasures(ii,2),sensorData.gps.positionMeasures(ii,3)]   =            ...
                                                 s.GPS_NEOM9N.sens( ...
                                                 sensorData.gps.positionMeasures(ii,1),...
                                                 sensorData.gps.positionMeasures(ii,2),...
                                               - sensorData.gps.positionMeasures(ii,3),...
                                                 14.8500);  
            [sensorData.gps.velocityMeasures(ii,1),sensorData.gps.velocityMeasures(ii,2),sensorData.gps.velocityMeasures(ii,3)] =           ...
                                                 s.GPS_NEOM9N.sens( ...
                                                 sensorData.gps.velocityMeasures(ii,1),...
                                                 sensorData.gps.velocityMeasures(ii,2),...
                                               - sensorData.gps.velocityMeasures(ii,3),...
                                                 14.8500);  
        end
        tot.gps_tot(tot.ngps_old:tot.ngps_old + size(sensorData.gps.positionMeasures,1) - 1,:)   =  sensorData.gps.positionMeasures(1:end,:) ;
        tot.gpsv_tot(tot.ngps_old:tot.ngps_old + size(sensorData.gps.velocityMeasures,1) - 1,:)  =  sensorData.gps.velocityMeasures(1:end,:) ;
        tot.time_gps(tot.ngps_old:tot.ngps_old + size(sensorData.gps.velocityMeasures,1) - 1)    =  sensorData.gps.time;
        tot.ngps_old = tot.ngps_old + size(sensorData.gps.positionMeasures,1);
end