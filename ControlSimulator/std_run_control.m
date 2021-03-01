function [Yf, Tf, t_ada, cpuTimes, flagMatr, dataBallisticFlight] = std_run_control(settings)
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

%% KALMAN PATH 
addpath('../kalman');

%% MAGNETIC FIELD MODEL
hmax = 6000;
% %Use this lines if your MATLAB version is up to 2020
dy = decyear(settings.launchDate);
[XYZ0,H0,D0,I0,F0] = wrldmagm(0, settings.lat0, settings.lon0, dy, '2020');
[XYZh] = wrldmagm(hmax, settings.lat0, settings.lon0, dy, '2020');
% fprintf('Horizontal intensity in nT: %g [nT] \n', H0)
% fprintf('Declination in degrees: %g [°]  \n', D0)
% fprintf('Inclination in degrees: %g [°]  \n', I0)
% fprintf('Total intensity in nT:  %g [nT] \n', F0)
    
% %Use this next line if your MATLAB version is previous to 2020
% load('magn_field.mat');

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
n_old     = 1;
np_old    = 1;
na_old    = 1;
ngps_old  = 1;
n_est_old = 1;
n_ada_old = 1;
cpuTimes = zeros(nmax,1);
iTimes = 0;
g = 9.81;
flagADA = false;

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

    
    %rotation of velocity to inertial axis
    A = zeros(3,3,size(Yf,1));
    v_NED = zeros(size(Yf,1),3);
     for ii=1:size(Yf,1)
         
         q               = [Yf(ii,11:13),Yf(ii,10)];        %composition of the quaternion in the form [q,q0]

         v_b             = Yf(ii,4:6)';                     %Velocity in body axis (column vector)
         
         q_mat           =  [0      -q(3)      q(2);        %Matrix needed for the
                             q(3)    0        -q(1);        %definition of rotation
                            -q(2)    q(1)      0   ;];      %matrix
                   
         A(:,:,ii)       =   (q(4)^2-q(1:3)*q(1:3)')*eye(3) + 2*q(1:3)'*q(1:3) - 2*q(4)*q_mat;               
                                                            %Rotation matrix to
                                                            %body axis from inertial
                                                            
         v_NED(ii,:)     = (A(:,:,ii)'*v_b)';               %Matrix with the velocity components in NED;
                                                            %each column is a velocity component
     end
   
     
    [sensorData] = manageSignalFrequencies(magneticFieldApprox, flagAscent, settings, Yf, Tf, x, uw, vw, ww, uncert);
    [~, ~, p, ~]  = atmosisa(-Yf(:,3)) ; 
 
 
    if settings.dataNoise

        % Baro Acquisition loop
        pn      = zeros(1,length(sensorData.barometer.time));
        h_baro  = zeros(1,length(sensorData.barometer.time));

        for ii=1:length(sensorData.barometer.time)
                pn(ii)        = MS580301BA01.sens(sensorData.barometer.measures(ii)/100,...
                                                  sensorData.barometer.temperature(ii) - 273.15);  
                h_baro(ii)    = -atmospalt(pn(ii)*100,'None');
        end 
        pn_tot(np_old:np_old + size(pn,2) - 1,1) = pn(1:end);
        hb_tot(np_old:np_old + size(pn,2) - 1,1) = h_baro(1:end);
        np_old = np_old + size(pn,2);
        
      
        % IMU Acquisition loop
        accel   = zeros(length(sensorData.accelerometer.time),3);
        gyro    = zeros(length(sensorData.gyro.time),3);
        mag     = zeros(length(sensorData.magnetometer.time),3);  
        
        for ii=1:length(sensorData.accelerometer.time)
                [accel(ii,1),accel(ii,2),accel(ii,3)] = ACCEL_LSM9DS1.sens(...
                                                 sensorData.accelerometer.measures(ii,1)*1000/g,...
                                                 sensorData.accelerometer.measures(ii,2)*1000/g,...
                                                 sensorData.accelerometer.measures(ii,3)*1000/g,...
                                                 14.8500);  
                 [gyro(ii,1),gyro(ii,2),gyro(ii,3)]   = GYRO_LSM9DS1.sens(...
                                                 sensorData.gyro.measures(ii,1)*1000*360/2/pi,...
                                                 sensorData.gyro.measures(ii,2)*1000*360/2/pi,...
                                                 sensorData.gyro.measures(ii,3)*1000*360/2/pi,...
                                                 14.8500);
                 [mag(ii,1),mag(ii,2),mag(ii,3)]      = MAGN_LSM9DS1.sens(...
                                                 sensorData.magnetometer.measures(ii,1)*0.01,...
                                                 sensorData.magnetometer.measures(ii,2)*0.01,...
                                                 sensorData.magnetometer.measures(ii,3)*0.01,...
                                                 14.8500);   
                 accel(ii,:) = accel(ii,:)*g/1000;
                 gyro(ii,:)  = gyro(ii,:)*2*pi/360/1000;
                                            
        end 
        accel_tot(na_old:na_old + size(accel,1) - 1,:) = accel(1:end,:) ;
        gyro_tot(na_old:na_old + size(gyro,1) - 1,:)   = gyro(1:end,:) ;
        mag_tot(na_old:na_old + size(mag,1) - 1,:)     = mag(1:end,:) ;
        na_old = na_old + size(accel,1);
        

        % GPS Acquisition loop
        gps     = zeros(length(sensorData.gps.time),3);
        gpsv    = zeros(length(sensorData.gps.time),3);
        
        for ii=1:length(sensorData.gps.time)
            [gps(ii,1),gps(ii,2),gps(ii,3)]   =  GPS_NEOM9N.sens(...
                                                 sensorData.gps.positionMeasures(ii,1),...
                                                 sensorData.gps.positionMeasures(ii,2),...
                                               - sensorData.gps.positionMeasures(ii,3),...
                                                 14.8500);  
            [gpsv(ii,1),gpsv(ii,2),gpsv(ii,3)] = GPS_NEOM9N.sens(...
                                                 sensorData.gps.velocityMeasures(ii,1),...
                                                 sensorData.gps.velocityMeasures(ii,2),...
                                               - sensorData.gps.velocityMeasures(ii,3),...
                                                 14.8500);  
        end
        gps_tot(ngps_old:ngps_old + size(gps,1) - 1,:)   =  gps(1:end,:) ;
        gpsv_tot(ngps_old:ngps_old + size(gpsv,1) - 1,:) =  gpsv(1:end,:) ;
        ngps_old = ngps_old + size(gps,1);

    end
  
    if iTimes==1
        x_prev    =  [X0; V0; Q0(2:4); Q0(1);0;0;0];
        P_prev    =   0.01*eye(12);
        ada_prev  =   settings.x0_ada;
        Pada_prev =   settings.P0_ada;
    else
        x_prev    =   x_est_tot(end,:);
        P_prev    =   P_c(:,:,end);
        ada_prev  =   x_ada_tot(end,:);
        Pada_prev =   P_ada(:,:,end);
    end
    
    %%%%%%%%%%%%% ADA %%%%%%%%%%%%%
    if flagADA == false
    [x_ada, P_ada, flagADA, t_ada]   =  run_ADA(ada_prev, Pada_prev, - h_baro, sensorData.barometer.time, settings.Q_ada, settings.N_ada);
    end
     x_ada_tot(n_ada_old:n_ada_old + size(x_ada(:,1),1)-1,:)  = x_ada(1:end,:);
     t_ada_tot(n_ada_old:n_ada_old + size(x_ada(:,1),1)-1)    = sensorData.barometer.time;              
     n_ada_old = n_ada_old + size(x_ada(1,:)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% kalmann filter %%%%%%%%
    n_satellite = 4;
    flagGPS_fix = true;
    [x_c,P_c]   =  run_kalman(x_prev,P_prev,...
                              sensorData.accelerometer.time, accel,...
                              gyro,...
                              sensorData.barometer.time, h_baro, settings.sigma_baro,...
                              sensorData.magnetometer.time, mag,settings.sigma_mag, XYZ0*0.01,...
                              sensorData.gps.time, gps,gpsv, settings.sigma_GPS,...
                              n_satellite,flagGPS_fix,settings.QLinear,settings.Qq);
     x_est_tot(n_est_old:n_est_old + size(x_c(:,1),1)-1,:)  = x_c(1:end,:);
     t_est_tot(n_est_old:n_est_old + size(x_c(:,1),1)-1) = sensorData.accelerometer.time;              
     n_est_old = n_est_old + size(x_c(1,:)); 
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
    
    v_NED_tot(n_old:n_old+n-1,:) = v_NED;
    
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

i_apo = find(Tf<24.8);
i_apo = max(i_apo);
i_apo_est = find(t_est_tot<Tf(i_apo));
i_apo_est=max(i_apo_est);

flagMatr = flagMatr(1:n_old, :);
%% RETRIVE PARAMETERS FROM THE ODE
if not(settings.electronics)
    dataBallisticFlight = RecallOdeFcn(@ascent, Tf(flagMatr(:, 2)), Yf(flagMatr(:, 2), :), settings, C, uw, vw, ww, uncert);
end

if true && not(settings.electronics)
%% FIGURE: Barometer reads 
fbaro = settings.frequencies.barometerFrequency;
tp = Tf(1):1/fbaro:Tf(end);
figure 
subplot(2,1,1);plot(tp,pn_tot',Tf,p_tot'/100);grid on;xlabel('time [s]');ylabel('|P| [mBar]');
legend('Pressure','Ground-truth','location','southeast');
title('Barometer pressure reads');
subplot(2,1,2);plot(tp,-hb_tot',Tf,-Yf(:,3));grid on;xlabel('time [s]');ylabel('|h| [m]');
legend('Altitude','Ground-truth','location','northeast');
title('Barometer altitude reads');
figure 
subplot(3,1,1);plot(t_ada_tot, x_ada_tot(:,1));grid on;xlabel('time [s]');ylabel('|P| [mBar]');
title('ADA pressure estimation');
subplot(3,1,2);plot(t_ada_tot,x_ada_tot(:,2));grid on;xlabel('time [s]');ylabel('|P_dot| [mBar/s]');
title('ADA velocity estimation');
subplot(3,1,3);plot(t_ada_tot,x_ada_tot(:,3));grid on;xlabel('time [s]');ylabel('|P_dot^2| [mBar/s^2]');
title('ADA acceleration estimation');
%% FIGURE: Accelerometer reads
faccel = settings.frequencies.accelerometerFrequency;
ta = Tf(1):1/faccel:Tf(end);
figure 
subplot(3,1,1);plot(ta,accel_tot(:,1)/g');grid on;xlabel('time [s]');ylabel('|ax| [g]');
title('Accelerometer reads along x');
subplot(3,1,2);plot(ta,accel_tot(:,2)/g');grid on;xlabel('time [s]');ylabel('|ay| [g]');
title('Accelerometer reads along y');
subplot(3,1,3);plot(ta,accel_tot(:,3)/g');grid on;xlabel('time [s]');ylabel('|az| [g]');
title('Accelerometer reads along z');
%% FIGURE: Gyroscope reads
figure 
subplot(3,1,1);plot(ta,gyro_tot(:,1)/1000');grid on;xlabel('time [s]');ylabel('|wx| [°/s]');
title('Gyroscope reads along x');
subplot(3,1,2);plot(ta,gyro_tot(:,2)/1000');grid on;xlabel('time [s]');ylabel('|wy| [°/s]');
title('Gyroscope reads along y');
subplot(3,1,3);plot(ta,gyro_tot(:,3)/1000');grid on;xlabel('time [s]');ylabel('|wz| [°/s]');
title('Gyroscope reads along z');
%% FIGURE: Magnetometer reads
figure 
subplot(3,1,1);plot(ta,mag_tot(:,1)');grid on;xlabel('time [s]');ylabel('|mx| [Gauss]');
title('Magnetometer reads along x');
subplot(3,1,2);plot(ta,mag_tot(:,2)');grid on;xlabel('time [s]');ylabel('|my| [Gauss]');
title('Magnetometer reads along y');
subplot(3,1,3);plot(ta,mag_tot(:,3)');grid on;xlabel('time [s]');ylabel('|mz| [Gauss]');
title('Magnetometer reads along z');
%% FIGURE: Gps reads
fgps = settings.frequencies.gpsFrequency;
tgps = Tf(1):1/fgps:Tf(end);
figure 
subplot(3,1,1);plot(tgps, gps_tot(:,1)');grid on;xlabel('time [s]');ylabel('|Pn| [m]');
title('GPS position  North');
subplot(3,1,2);plot(tgps, gps_tot(:,2)');grid on;xlabel('time [s]');ylabel('|Pe| [m]');
title('GPS position  East');
subplot(3,1,3);plot(tgps,-gps_tot(:,3)');grid on;xlabel('time [s]');ylabel('|Pu| [m]');
title('GPS position Upward');
figure 
subplot(3,1,1);plot(tgps,gpsv_tot(:,1)');grid on;xlabel('time [s]');ylabel('|Velocity N| [m/s]');
subplot(3,1,2);plot(tgps,gpsv_tot(:,2)');grid on;xlabel('time [s]');ylabel('|Velocity E| [m/s]');
subplot(3,1,3);plot(tgps,gpsv_tot(:,3)');grid on;xlabel('time [s]');ylabel('|Velocity D| [m/s]');
title('GPS velocity reads');
subplot(3,1,1);plot(tgps, gpsv_tot(:,1)');grid on;xlabel('time [s]');ylabel('|Vn| [m/s]');
title('GPS velocity  North');
subplot(3,1,2);plot(tgps, gpsv_tot(:,2)');grid on;xlabel('time [s]');ylabel('|Ve| [m/s]');
title('GPS velocity  East');
subplot(3,1,3);plot(tgps,-gpsv_tot(:,3)');grid on;xlabel('time [s]');ylabel('|Vu| [m/s]');
title('GPS velocity Upward');
%% FIGURE: Estimated position vs ground-truth
figure
subplot(3,1,1);plot(t_est_tot(1:i_apo_est), x_est_tot(1:i_apo_est,1),Tf(1:i_apo), Yf(1:i_apo,1));grid on;xlabel('time [s]');ylabel('|Pn| [m]');
legend('North','Ground-truth','location','best');
title('Estimated North position vs ground-truth');
subplot(3,1,2);plot(t_est_tot(1:i_apo_est), x_est_tot(1:i_apo_est,2),Tf(1:i_apo), Yf(1:i_apo,2));grid on;xlabel('time [s]');ylabel('|Pe| [m]');
legend('East','Ground-truth','location','best');
title('Estimated East position vs ground-truth');
subplot(3,1,3);plot(t_est_tot(1:i_apo_est), -x_est_tot(1:i_apo_est,3),Tf(1:i_apo), -Yf(1:i_apo,3));grid on;xlabel('time [s]');ylabel('|Pu| [m]');
legend('Upward','Ground-truth','location','best');
title('Estimated Upward position vs ground-truth');
%% FIGURE: Estimated velocities vs ground-truth
figure
subplot(3,1,1);plot(t_est_tot(1:i_apo_est), x_est_tot(1:i_apo_est,4),Tf(1:i_apo), v_NED_tot(1:i_apo,1));grid on;xlabel('time [s]');ylabel('|Vn| [m/s]');
legend('North','Ground-truth','location','best');
title('Estimated North velocity vs ground-truth');
subplot(3,1,2);plot(t_est_tot(1:i_apo_est), x_est_tot(1:i_apo_est,5),Tf(1:i_apo), v_NED_tot(1:i_apo,2));grid on;xlabel('time [s]');ylabel('|Ve| [m/s]');
legend('East','Ground-truth','location','best');
title('Estimated East velocity vs ground-truth');
subplot(3,1,3);plot(t_est_tot(1:i_apo_est+1),-x_est_tot(1:i_apo_est+1,6),Tf(1:i_apo),-v_NED_tot(1:i_apo,3));grid on;xlabel('time [s]');ylabel('|Vu| [m/s]');
legend('Upward','Ground-truth','location','best');
title('Estimated Upward velocity vs ground-truth');
%% FIGURE: Estimated quaternions vs ground-truth
figure
subplot(4,1,1);plot(t_est_tot(1:i_apo_est),x_est_tot(1:i_apo_est,10),Tf(1:i_apo),Yf(1:i_apo,10));grid on;ylabel('|q0| [-]');
legend('Estimated q0','Ground-truth','location','northeast');
title('Estimated q0 vs ground-truth');
subplot(4,1,2);plot(t_est_tot(1:i_apo_est),x_est_tot(1:i_apo_est,7),Tf(1:i_apo),Yf(1:i_apo,11));grid on;ylabel('|q1| [-]');
legend('Estimated q1','Ground-truth','location','northeast');
title('Estimated q1 vs ground-truth');
subplot(4,1,3);plot(t_est_tot(1:i_apo_est),x_est_tot(1:i_apo_est,8),Tf(1:i_apo),Yf(1:i_apo,12));grid on;ylabel('|q2| [-]');
legend('Estimated q2','Ground-truth','location','northeast');
title('Estimated q2 vs ground-truth');
subplot(4,1,4);plot(t_est_tot(1:i_apo_est),x_est_tot(1:i_apo_est,9),Tf(1:i_apo),Yf(1:i_apo,13));grid on;ylabel('|q3| [-]');
legend('Estimated q3','Ground-truth','location','northeast');
title('Estimated q3 vs ground-truth');
end
end


