function [s, sensorCounter, sensorProcessed]  = initSensors
% Initialize all sensors
    
% Author: Jan Hammelman
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@skywarder.eu
% email: jan.hammelmann@skywarder.eu,alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

% initial barometer sensor MS580301BA01
ep_p_0     =  csvread('skyward-matlab-control-simulator/sensors/data/MS580301BA01/ep_p_0.csv');
ep_p_25    =  csvread('skyward-matlab-control-simulator/sensors/data/MS580301BA01/ep_p_25.csv');
ep_p_85    =  csvread('skyward-matlab-control-simulator/sensors/data/MS580301BA01/ep_p_85.csv');
ep_p_neg40 =  csvread('skyward-matlab-control-simulator/sensors/data/MS580301BA01/ep_p_-40.csv');

p_table = [ep_p_0(:,1);ep_p_25(:,1);ep_p_85(:,1);ep_p_neg40(:,1)];         % presure in mbar
ep      = [ep_p_0(:,2);ep_p_25(:,2);ep_p_85(:,2);ep_p_neg40(:,2)];         % error presure in mbar
T       = [0*ones(size(ep_p_0(:,1)));25*ones(size(ep_p_25(:,1)));85*ones(size(ep_p_85(:,1)));-40*ones(size(ep_p_neg40(:,1)));];
ep_data = [p_table,T,ep];

s.MS580301BA01 = Sensor(); % presure in mbar, temp should be in C°
s.MS580301BA01.maxMeasurementRange  =   1100;                   % 1100, 1300 in mbar
s.MS580301BA01.minMeasurementRange  =   300;                    % 300, 10 in mbar
s.MS580301BA01.resolution           =   0.012;                  % 0.012, 0.018, 0.027, 0.042, 0.065 in mbar
s.MS580301BA01.noiseVariance        =   4;                      % guess in mbar
s.MS580301BA01.error2dOffset        =   ep_data;                % [p in mbar, T in celsius, ep in mbar]

% initial accelerometer sensor from LSM9DS1
s.ACCEL_LSM9DS1 = Sensor3D(); % acceleration in mg
s.ACCEL_LSM9DS1.maxMeasurementRange =   16000;                  % 2000, 4000, 8000, 16000 in mg
s.ACCEL_LSM9DS1.minMeasurementRange =  -16000;                  % -2000, -4000, -8000, -16000 in mg
s.ACCEL_LSM9DS1.resolution          =   0.732;                  % 0.061, 0.122, 0.244, 0.732 in mg 
s.ACCEL_LSM9DS1.noiseVariance       =   4;                      % guess in mg original was 4
s.ACCEL_LSM9DS1.offsetX             =   0;                      % +-90 in mg
s.ACCEL_LSM9DS1.offsetY             =   0;                      % +-90 in mg
s.ACCEL_LSM9DS1.offsetZ             =   0;                      % +-90 in mg
s.ACCEL_LSM9DS1.walkDiffusionCoef   =   1;                      % guess
s.ACCEL_LSM9DS1.dt                  =   0.01;                   % sampling time

% initial gyroscope sensor from LSM9DS1
s.GYRO_LSM9DS1 = Sensor3D(); % angular rate in mdps
s.GYRO_LSM9DS1.maxMeasurementRange  =   245e3;                  % 245e3, 500e3, 2000e3 in mdps
s.GYRO_LSM9DS1.minMeasurementRange  =  -245e3;                  % -245e3, -500e3, -2000e3 in mdps
s.GYRO_LSM9DS1.resolution           =   8.75;                   % 8.75, 17.5, 70 in mdps
s.GYRO_LSM9DS1.noiseVariance        =   100;                    % guess in mdps    100 was original
s.GYRO_LSM9DS1.offsetX              =   0;                      % +-30e3 in mdps
s.GYRO_LSM9DS1.offsetY              =   0;                      % +-30e3 in mdps
s.GYRO_LSM9DS1.offsetZ              =   0;                      % +-30e3 in mdps
s.GYRO_LSM9DS1.walkDiffusionCoef    =   1;                      % guess
s.GYRO_LSM9DS1.dt                   =   0.01;                   % sampling time
s.GYRO_LSM9DS1.transMatrix          =   diag([1 1 1]);          % axis transformation

% initial megnetometer sensor from LSM9DS1
s.MAGN_LSM9DS1=Sensor3D(); % magnetic field in mgauss
s.MAGN_LSM9DS1.maxMeasurementRange  =   16000;                  % 4000, 8000, 12000, 16000 in mgauss
s.MAGN_LSM9DS1.minMeasurementRange  =  -16000;                  % -4000, -8000, -12000, -16000 in mgauss
s.MAGN_LSM9DS1.resolution           =   0.58;                   % 0.14, 0.29, 0.43, 0.58 in mgauss
s.MAGN_LSM9DS1.noiseVariance        =   4;                      % guess in mgauss    original guess 2
s.MAGN_LSM9DS1.offsetX              =   0;                      % +-1000 in mgauss
s.MAGN_LSM9DS1.offsetY              =   0;                      % +-1000 in mgauss
s.MAGN_LSM9DS1.offsetZ              =   0;                      % +-1000 in mgauss
s.MAGN_LSM9DS1.walkDiffusionCoef    =   1;                      % guess
s.MAGN_LSM9DS1.dt                   =   0.01;                   % sampling time
s.MAGN_LSM9DS1.transMatrix          =   diag([1 1 1]);          % axis transformation

% initial GPS sensor from NEO-M9N
s.GPS_NEOM9N = Sensor3D();                                      % lon, in degree lat in deree, alt in m
s.GPS_NEOM9N.noiseVariance          =   2;                      % in m
s.GPS_NEOM9N.transMatrix            =   diag([1 1 1]);          % axis transformation

% initial megnetometer sensor from IIS2MDC: TODO
s.MAGN_IIS2MDC = Sensor3D();                                    % magnetic field in mgauss, temp should be in C°-25C°
s.MAGN_IIS2MDC.maxMeasurementRange  =   49152;                  % in mgauss
s.MAGN_IIS2MDC.minMeasurementRange  =   -49152;                 % in mgauss
s.MAGN_IIS2MDC.resolution           =   1.5;                    % in mgauss
s.MAGN_IIS2MDC.tempOffset           =   0;                      % +-0.3 in mgauss
s.MAGN_IIS2MDC.noiseVariance        =   90;                     % guess in mgauss
s.MAGN_IIS2MDC.offsetX              =   0;                      % +-60 in mgauss
s.MAGN_IIS2MDC.offsetY              =   0;                      % +-60 in mgauss
s.MAGN_IIS2MDC.offsetZ              =   0;                      % +-60 in mgauss
s.MAGN_IIS2MDC.walkDiffusionCoef    =   1;                      % guess
s.MAGN_IIS2MDC.dt                   =   0.01;                   % sampling time
s.MAGN_IIS2MDC.transMatrix          =   diag([1 1 1]);          % axis transformation
s.MAGN_IIS2MDC.transMatrix          =   diag([1 1 1]);          % axis transformation

sensorCounter.np_old        =   1;
sensorCounter.na_old        =   1;
sensorCounter.ngps_old      =   1;
sensorCounter.n_est_old     =   1;
sensorCounter.n_ada_old     =   1;

sensorProcessed.pn_tot      =   0;
sensorProcessed.hb_tot      =   0;
sensorProcessed.accel_tot   =   [0, 0, 0];
sensorProcessed.gyro_tot    =   [0, 0, 0];
sensorProcessed.mag_tot     =   [0, 0, 0];
sensorProcessed.gps_tot     =   [0, 0, 0];
sensorProcessed.gpsv_tot    =   [0, 0, 0];