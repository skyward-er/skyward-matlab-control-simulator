function [Yf, Tf, cpuTimes, flagMatr, dataBallisticFlight] = std_run_control(settings)
%{

STD_RUN_BALLISTIC - This function runs a standard ballistic (non-stochastic) simulation

INTPUTS:
            - settings, rocket data structure;

OUTPUTS:
            - Tf, Total integration time vector;
            - Yf, Total State Matrix;

Author: Ruben Di Battista
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu

Author: Francesco Colombi
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: francesco.colombi@skywarder.eu

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept
email: adriano.filippo.inno@skywarder.eu
Revision date: 09/10/2019

%}

if not(settings.ballisticFligth) && settings.ascentOnly
   error('To simulate a landing with the parachutes, settings.ascentOnly must be false') 
end

%% STARTING CONDITIONS
% Attitude
Q0 = angle2quat(settings.PHI, settings.OMEGA, 0*pi/180, 'ZYX')';

% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
Y0 = [X0; V0; W0; Q0; settings.Ixxf; settings.Iyyf; settings.Izzf];

%% WIND GENERATION
if settings.wind.input   % will be computed inside the integrations
    uw = 0; vw = 0; ww = 0;
else
    [uw,vw,ww,~] = wind_const_generator(settings.wind.AzMin, settings.wind.AzMax,...
        settings.wind.ElMin, settings.wind.ElMax, settings.wind.MagMin, settings.wind.MagMax);
    
    if ww ~= 0
        warning('Pay attention using vertical wind, there might be computational errors')
    end
    
end

if settings.wind.input && all(settings.wind.input_uncertainty ~= 0)
    signn = randi([1, 4]); % 4 sign cases
    unc = settings.wind.input_uncertainty;
    
    switch signn
        case 1
            %                       unc = unc;
        case 2
            unc(1) = - unc(1);
        case 3
            unc(2) = - unc(2);
        case 4
            unc = - unc;
    end
    
    uncert = rand(1,2).*unc;
else
    uncert = [0,0];
end

%% MAGNETIC FIELD MODEL
dy = decyear(settings.launchDate);
hmax = 6000;
[XYZ0] = wrldmagm(0, settings.lat0, settings.lon0, dy, '2020');
[XYZh] = wrldmagm(hmax, settings.lat0, settings.lon0, dy, '2020');

magneticFieldApprox = @(zSlm) XYZ0 + (XYZh-XYZ0)./hmax.*zSlm;
%% SENSORS DEFINITION
addpath('../sensors');
addpath('../sensors/data/MS580301BA01');
initSensors;
%% INTEGRATION
% setting initial condition before control phase
dt = 1/settings.frequencies.controlFrequency;
t0 = 0;
t1 = t0 + dt;
vz = 1;
z = 1;
flagStopIntegration = true;
nmax = 10000;
mach = 0;
x = 0;
flagMatr = false(nmax, 6);
flagAscent = false;
Yf_tot = zeros(nmax, 16);
Tf_tot = zeros(nmax, 1);
p_tot  = zeros(nmax, 1);
C = zeros(nmax, 1);
n_old = 1;
np_old = 1;
na_old = 1;
ngps_old = 1;
cpuTimes = zeros(nmax,1);
iTimes = 0;
g = 9.81;

while flagStopIntegration || n_old < nmax
    tic 
    iTimes = iTimes + 1;
    
    lastFlagAscent = flagAscent;

    if t0 <= settings.tb
        flagBurning = true;
    else
        flagBurning = false;
    end
    
    if flagAscent && not(flagBurning) && mach <=0.7
        flagAeroBrakes = true;
    else
        flagAeroBrakes = false;
    end
    
    if z < 0
        flagFligth = false;
    else
        flagFligth = true;
    end
    
    if vz >= 0
        flagAscent = true;
    else
        flagAscent = false;
    end
    
    if not(flagAscent) 
        if z >= settings.para(1).z_cut
            flagPara1 = true;
            flagPara2 = false;
        else
            flagPara1 = false;
            flagPara2 = true;
        end
    else
        flagPara1 = false;
        flagPara2 = false;
    end
    
    % dynamics
    if settings.ballisticFligth
        [Tf, Yf] = ode45(@ascent, [t0, t1], Y0, [], settings, x, uw, vw, ww, uncert);
    else
        if flagAscent
            [Tf, Yf] = ode45(@ascent, [t0, t1], Y0, [], settings, x, uw, vw, ww, uncert);
        else
            if flagPara1 
                para = 1; 
            end
            if flagPara2 
                para = 2; 
            end
            
            Y0 = Y0(1:6);
            [Tf, Yd] = ode45(@descentParachute, [t0, t1], Y0, [], settings, uw, vw, ww, para, uncert);
            [nd, ~] = size(Yd);
            Yf = [Yd, zeros(nd, 7), settings.Ixxe*ones(nd, 1), ...
                settings.Iyye*ones(nd, 1), settings.Iyye*ones(nd, 1)];
        end
    end

    [sensorData] = manageSignalFrequencies(magneticFieldApprox, flagAscent, settings, Yf, Tf, x, uw, vw, ww, uncert);
    [~, ~, p, ~]  = atmosisa(-Yf(:,3)) ; 
    
    if settings.dataNoise
        
        for ii=1:length(sensorData.barometer.time)
                pn(ii,1) = MS580301BA01.sens(sensorData.barometer.measures(ii)/100,...
                                           sensorData.barometer.temperature(ii) - 273.15);  
        end 
        pn_tot(np_old:np_old + length(pn) - 1,1) = pn(1:end,1)';
        np_old = np_old + length(pn);
        
        for ii=1:length(sensorData.accelerometer.time)
                [accelx(ii,1),accely(ii,1),accelz(ii,1)] = ACCEL_LSM9DS1.sens(...
                                                 sensorData.accelerometer.measures(ii,1)*1000/g,...
                                                 sensorData.accelerometer.measures(ii,2)*1000/g,...
                                                 sensorData.accelerometer.measures(ii,3)*1000/g,...
                                                 14.8500);  
                 [gyrox(ii,1),gyroy(ii,1),gyroz(ii,1)] = GYRO_LSM9DS1.sens(...
                                                 sensorData.gyro.measures(ii,1)*1000*360/2/pi,...
                                                 sensorData.gyro.measures(ii,2)*1000*360/2/pi,...
                                                 sensorData.gyro.measures(ii,3)*1000*360/2/pi,...
                                                 14.8500);  
                 [magx(ii,1),magy(ii,1),magz(ii,1)] = MAGN_LSM9DS1.sens(...
                                                 sensorData.magnetometer.measure(ii,1)*0.01,...
                                                 sensorData.magnetometer.measure(ii,2)*0.01,...
                                                 sensorData.magnetometer.measure(ii,3)*0.01,...
                                                 14.8500);  
                                            
        end 
        accel_tot(na_old:na_old + length(accelx) - 1,:) =[accelx, accely, accelz] ;
        gyro_tot(na_old:na_old + length(gyrox) - 1,:) =[gyrox, gyroy, gyroz] ;
        mag_tot(na_old:na_old + length(magx) - 1,:) =[magx, magy, magz] ;
        na_old = na_old + length(accelx);
        
        for ii=1:length(sensorData.gps.time)
            [gpsx(ii,1),gpsy(ii,1),gpsz(ii,1)] = GPS_NEOM9N.sens(...
                                                 sensorData.gps.positionMeasures(ii,1),...
                                                 sensorData.gps.positionMeasures(ii,2),...
                                                 sensorData.gps.positionMeasures(ii,3),...
                                                 14.8500);  
            [gpsvx(ii,1),gpsvy(ii,1),gpsvz(ii,1)] = GPS_NEOM9N.sens(...
                                                 sensorData.gps.velocityMeasures(ii,1),...
                                                 sensorData.gps.velocityMeasures(ii,2),...
                                                 sensorData.gps.velocityMeasures(ii,3),...
                                                 14.8500);  
        end
        gps_tot(ngps_old:ngps_old + length(gpsx) - 1,:) =[gpsx, gpsy, gpsz] ;
        gpsv_tot(ngps_old:ngps_old + length(gpsvx) - 1,:) =[gpsvx, gpsvy, gpsvz] ;
        ngps_old = ngps_old + length(gpsx);

    end
  

    
    %%%%%%% kalmann filter %%%%%%%%
    % kalman(p, acc_body, ang_vel, q, [u, v, w]ned, [x, y, z]ned )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flagAeroBrakes
         alpha_degree = controlAlgorithm(z, vz, vx);
         x = get_extension_from_angle(alpha_degree);
    else 
        x = 0;
    end    

    % vertical velocity and position
    if flagAscent || (not(flagAscent) && settings.ballisticFligth)
        Q = Yf(end, 10:13);
        vels = quatrotate(quatconj(Q), Yf(end, 4:6));
        vz = - vels(3);
        vx = vels(1); % Needed for the control algorithm. Ask if it is right
    else
        vz = -Yf(end, 6);
        vx = Yf(end, 4);  % Needed for the control algorithm. Ask if it is right
    end
    z = -Yf(end, 3);

    if lastFlagAscent && not(flagAscent)
        Y0 = [Yf(end, 1:3), vels, Yf(end, 7:end)];
    else
        Y0 = Yf(end, :);
    end
    
    % atmosphere
    [~, a, ~, ~] = atmosisa(z);        % pressure and temperature at each sample time
    normV = norm(Yf(end, 4:6));
    mach = normV/a;
    
    % time update
    t0 = t0 + dt;
    t1 = t1 + dt;
    
    % assemble total state
    [n, ~] = size(Yf);
    Yf_tot(n_old:n_old+n-1, :) = Yf(1:end, :);
    Tf_tot(n_old:n_old+n-1,1) = Tf(1:end, 1);
    p_tot(n_old:n_old+n-1,1) = p(1:end, 1);
    C(n_old:n_old+n-1) = x;
    
    n_old = n_old + n -1;
   
    cpuTimes(iTimes) = toc;
    
     if settings.ascentOnly
         flagStopIntegration = flagAscent;
     else
         flagStopIntegration = flagFligth;
     end        
    
     flagMatr(n_old:n_old+n-1, :) = repmat([flagFligth, flagAscent, flagBurning, flagAeroBrakes, flagPara1, flagPara2], n, 1);
end
cpuTimes = cpuTimes(1:iTimes);

%% ASSEMBLE TOTAL FLIGHT STATE
Yf = Yf_tot(1:n_old, :);
Tf = Tf_tot(1:n_old, :);
p = p_tot(1:n_old, :);
fbaro = settings.frequencies.barometerFrequency;
faccel = settings.frequencies.accelerometerFrequency;
tp = Tf(1):1/fbaro:Tf(end);
ta = Tf(1):1/faccel:Tf(end);
flagMatr = flagMatr(1:n_old, :);

%% RETRIVE PARAMETERS FROM THE ODE
if not(settings.electronics)
    dataBallisticFlight = RecallOdeFcn(@ascent, Tf(flagMatr(:, 2)), Yf(flagMatr(:, 2), :), settings, C, uw, vw, ww, uncert);
end

figure 
subplot(2,1,1)
plot (tp,pn_tot'), grid on;
xlabel('time [s]'), ylabel('|P| [mBar]');
subplot(2,1,2)
plot (Tf,p_tot'/100), grid on;

figure 
subplot(3,1,1)
plot (ta,accel_tot(:,1)/1000'), grid on;
subplot(3,1,2)
plot (ta,accel_tot(:,2)/1000'), grid on;
subplot(3,1,3)
plot (ta,accel_tot(:,3)/1000'), grid on;
xlabel('time [s]'), ylabel('|Acc| [mg]');

figure 
subplot(3,1,1)
plot (ta,gyro_tot(:,1)/1000'), grid on;
subplot(3,1,2)
plot (ta,gyro_tot(:,2)/1000'), grid on;
subplot(3,1,3)
plot (ta,gyro_tot(:,3)/1000'), grid on;
xlabel('time [s]'), ylabel('|Ang vel| [°/s]');

figure 
subplot(3,1,1)
plot (ta,mag_tot(:,1)/1000'), grid on;
subplot(3,1,2)
plot (ta,mag_tot(:,2)/1000'), grid on;
subplot(3,1,3)
plot (ta,mag_tot(:,3)/1000'), grid on;
xlabel('time [s]'), ylabel('|Mag field| [gauss]');
end