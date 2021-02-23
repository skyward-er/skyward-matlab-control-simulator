% Initialize all sensors


% initial barometer sensor MS580301BA01
ep_p_0=csvread('skyward-matlab-control-simulator/sensors/data/MS580301BA01/ep_p_0.csv');
ep_p_25=csvread('skyward-matlab-control-simulator/sensors/data/MS580301BA01/ep_p_25.csv');
ep_p_85=csvread('skyward-matlab-control-simulator/sensors/data/MS580301BA01/ep_p_85.csv');
ep_p_neg40=csvread('skyward-matlab-control-simulator/sensors/data/MS580301BA01/ep_p_-40.csv');

p_table=[ep_p_0(:,1);ep_p_25(:,1);ep_p_85(:,1);ep_p_neg40(:,1)]; % presure in mbar
ep=[ep_p_0(:,2);ep_p_25(:,2);ep_p_85(:,2);ep_p_neg40(:,2)]; % error presure in mbar
T=[0*ones(size(ep_p_0(:,1)));25*ones(size(ep_p_25(:,1)));85*ones(size(ep_p_85(:,1)));-40*ones(size(ep_p_neg40(:,1)));];
ep_data=[p_table,T,ep];

MS580301BA01=Sensor(); % presure in mbar, temp should be in C°
MS580301BA01.maxMeasurementRange=1100; % 1100, 1300 in mbar
MS580301BA01.minMeasurementRange=300; % 300, 10 in mbar
MS580301BA01.resolution=0.065; % 0.012, 0.018, 0.027, 0.042, 0.065 in mbar
MS580301BA01.noiseVariance=4; % guess in mbar
MS580301BA01.error2dOffset=ep_data; % [p in mbar, T in celsius, ep in mbar]

% initial accelerometer sensor from LSM9DS1
ACCEL_LSM9DS1=Sensor3D(); % acceleration in mg
ACCEL_LSM9DS1.maxMeasurementRange=16000; % 2000, 4000, 8000, 16000 in mg
ACCEL_LSM9DS1.minMeasurementRange=-16000; % -2000, -4000, -8000, -16000 in mg
ACCEL_LSM9DS1.resolution=0.732; % 0.061, 0.122, 0.244, 0.732 in mg 
ACCEL_LSM9DS1.noiseVariance=4; % guess in mg
ACCEL_LSM9DS1.offsetX=0; % +-90 in mg
ACCEL_LSM9DS1.offsetY=0; % +-90 in mg
ACCEL_LSM9DS1.offsetZ=0; % +-90 in mg
ACCEL_LSM9DS1.walkDiffusionCoef=1; % guess
ACCEL_LSM9DS1.dt=0.01; % sampling time

% initial gyroscope sensor from LSM9DS1
GYRO_LSM9DS1=Sensor3D(); % angular rate in mdps
GYRO_LSM9DS1.maxMeasurementRange=2000e3; % 245e3, 500e3, 2000e3 in mdps
GYRO_LSM9DS1.minMeasurementRange=-2000e3; % -245e3, -500e3, -2000e3 in mdps
GYRO_LSM9DS1.resolution=70; % 8.75, 17.5, 70 in mdps
GYRO_LSM9DS1.noiseVariance=100; % guess in mdps
GYRO_LSM9DS1.offsetX=0; % +-30e3 in mdps
GYRO_LSM9DS1.offsetY=0; % +-30e3 in mdps
GYRO_LSM9DS1.offsetZ=0; % +-30e3 in mdps
GYRO_LSM9DS1.walkDiffusionCoef=1; % guess
GYRO_LSM9DS1.dt=0.01; % sampling time
GYRO_LSM9DS1.transMatrix=diag([1 1 1]); % axis transformation

% initial megnetometer sensor from LSM9DS1
MAGN_LSM9DS1=Sensor3D(); % magnetic field in mgauss
MAGN_LSM9DS1.maxMeasurementRange=16000; % 4000, 8000, 12000, 16000 in mgauss
MAGN_LSM9DS1.minMeasurementRange=-16000; % -4000, -8000, -12000, -16000 in mgauss
MAGN_LSM9DS1.resolution=0.58; % 0.14, 0.29, 0.43, 0.58 in mgauss
MAGN_LSM9DS1.noiseVariance=2; % guess in mgauss
MAGN_LSM9DS1.offsetX=0; % +-1000 in mgauss
MAGN_LSM9DS1.offsetY=0; % +-1000 in mgauss
MAGN_LSM9DS1.offsetZ=0; % +-1000 in mgauss
MAGN_LSM9DS1.walkDiffusionCoef=1; % guess
MAGN_LSM9DS1.dt=0.01; % sampling time
MAGN_LSM9DS1.transMatrix=diag([1 1 1]); % axis transformation

% initial GPS sensor from NEO-M9N
GPS_NEOM9N=Sensor3D(); % lon, in degree lat in deree, alt in m
GPS_NEOM9N.noiseVariance=4; % in m
GPS_NEOM9N.transMatrix=diag([1 1 1]); % axis transformation

% initial megnetometer sensor from LSM9DS1
MAGN_IIS2MDC=Sensor3D(); % magnetic field in mgauss, temp should be in C°-25C°
MAGN_IIS2MDC.maxMeasurementRange=49152; % in mgauss
MAGN_IIS2MDC.minMeasurementRange=-49152; % in mgauss
MAGN_IIS2MDC.resolution=1.5; % in mgauss
MAGN_IIS2MDC.tempOffset=0; % +-0.3 in mgauss
MAGN_IIS2MDC.noiseVariance=9; % guess in mgauss
MAGN_IIS2MDC.offsetX=0; % +-60 in mgauss
MAGN_IIS2MDC.offsetY=0; % +-60 in mgauss
MAGN_IIS2MDC.offsetZ=0; % +-60 in mgauss
MAGN_IIS2MDC.walkDiffusionCoef=1; % guess
MAGN_IIS2MDC.dt=0.01; % sampling time
MAGN_IIS2MDC.transMatrix=diag([1 1 1]); % axis transformation

