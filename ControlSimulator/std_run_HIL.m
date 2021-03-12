function [Yf, Tf, cpuTimes, flagMatr] = std_run_HIL(settings)
%{

STD_RUN_HIL - This function runs a hardware-in-the-loop simulation

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

Author: Emilio Corigliano
Skyward Experimental Rocketry | ELC Dept | electronics@skywarder.eu
email: emilio.corigliano@skywarder.eu
Release date: 10/03/2021

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
initialCond = [X0; V0; W0; Q0; settings.Ixxf; settings.Iyyf; settings.Izzf];
Y0 = initialCond;

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

%% SENSORS DEFINITION
addpath('../sensors');
addpath('../sensors/data/MS580301BA01');
initSensors;

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

%% INTEGRATION
% setting initial condition before control phase
dt = 1/settings.frequencies.controlFrequency;
t0 = 0;
t1 = t0 + dt;
vz = 1;
z = 1;
nmax = 10000;
mach = 0;
x = 0;
Yf_tot = zeros(nmax, 16);
Tf_tot = zeros(nmax, 1);
C = zeros(nmax, 1);
n_old = 1;
cpuTimes = zeros(nmax,1);
iTimes = 0;
g = 9.81;

index_plot = 1;
plot_control_variable = zeros(nmax,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flagStopIntegration = true;
flagMatr = false(nmax, 6);
flagFligth = false;
flagAscent = false;
flagBurning = false;
flagAeroBrakes = false;
flagPara1 = false;
flagPara2 = false;

if settings.launchWindow
    launchWindow;
    pause(0.01);
    launchFlag = false;
    lastLaunchflag = true;
else
    launchFlag = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while flagStopIntegration || n_old < nmax
    tic 
    iTimes = iTimes + 1;
    
    lastFlagAscent = flagAscent;

    if settings.launchWindow
        if not(lastLaunchflag) && launchFlag
            tLaunch = t0;
        end
    else 
        tLaunch = 0;
    end
    
    if launchFlag && (t0 - tLaunch) <= settings.tb
        flagBurning = true;
    else
        flagBurning = false;
    end
    
    if flagAscent && not(flagBurning) && mach <=0.7
        flagAeroBrakes = true;
    else
        flagAeroBrakes = false;
    end
    
    if z < 0 || not(launchFlag)
        flagFligth = false;
    else
        flagFligth = true;
    end
    
    if vz >= 0 && launchFlag
        flagAscent = true;
    else
        flagAscent = false;
    end
    
    if not(flagAscent) && launchFlag
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
    % dynamics
    if flagFligth
        if settings.ballisticFligth
            [Tf, Yf] = ode45(@ascent, [t0, t1], Y0, [], settings, x, uw, vw, ww, uncert, tLaunch);
        else
            if flagAscent
                [Tf, Yf] = ode45(@ascent, [t0, t1], Y0, [], settings, x, uw, vw, ww, uncert, tLaunch);
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
    else
        Tf = [t0, t1];
        Yf = [initialCond'; initialCond'];
    end
   
    [sensorData] = manageSignalFrequencies(magneticFieldApprox, flagAscent, settings, Yf, Tf, x, uw, vw, ww, uncert);
    %[~, ~, p, ~]  = atmosisa(-Yf(:,3)) ; 
    
    sensorData.barometer.measures(1)
 
    if settings.dataNoise
        % QUA  MI SERVE CHE IL NOISE VENGA APPLICATO AI DATI DENTRO A
        % SENSOR DATA ALTRIMENTI NON POSSO MANDARLI SU SERIALE
    end
    
    %%%%%%%%%%%
    % TEMPORARY SOLUTION UNTIL WE DON'T HAVE THE OBSW KALMAN
    if flagAeroBrakes
         sensorData.kalman.z    = z;
         sensorData.kalman.vz   = vz;
         sensorData.kalman.vMod = normV;
    else 
        sensorData.kalman.z    = 0;
        sensorData.kalman.vz   = 0;
        sensorData.kalman.vMod = 0;
    end 
    %%%%%%%%%%
    
    flagsArray = [flagFligth, flagAscent, flagBurning, flagAeroBrakes, flagPara1, flagPara2];
    
    sendDataOverSerial(sensorData, flagsArray);
    %alpha_degree = readControlOutputFromSerial();
    alpha_degree = 25;
         
    if flagAeroBrakes
         x = get_extension_from_angle(alpha_degree);
    else 
        x = 0;
    end 


    % vertical velocity and position
    if flagAscent || (not(flagAscent) && settings.ballisticFligth)
        Q = Yf(end, 10:13);
        vels = quatrotate(quatconj(Q), Yf(end, 4:6)); 
        vz = - vels(3);   % down
        vx = vels(2);   % north
        vy = vels(1);   % east
    else
        vz = -Yf(end, 6);  
        vx = Yf(end, 5); 
        vy = Yf(end, 4); 
    end
    z = -Yf(end, 3);
    xxx = Yf(end, 2);
    yyy = Yf(end, 1);
    

    if lastFlagAscent && not(flagAscent)
        Y0 = [Yf(end, 1:3), vels, Yf(end, 7:end)];
    else
        Y0 = Yf(end, :);
    end
    
    % atmosphere
    [~, a, ~, ~] = atmosisa(z); % pressure and temperature at each sample time
    normV = norm(Yf(end, 4:6));
    mach = normV / a;
    
    % time update
    t0 = t0 + dt;
    t1 = t1 + dt;
    
    % assemble total state
    [n, ~] = size(Yf);
    Yf_tot(n_old:n_old+n-1, :) = Yf(1:end, :);
    Tf_tot(n_old:n_old+n-1,1) = Tf(1:end, 1);
    C(n_old:n_old+1) = x;
    
    n_old = n_old + n -1;
   
    cpuTimes(iTimes) = toc;
    
    if settings.launchWindow
        lastLaunchflag = launchFlag;
        pause(1e-6);
        if exist('launchFlag.txt','file') == 2
            launchFlag = true;
        end
    end
    
     if settings.ascentOnly
         flagStopIntegration = flagAscent;
     else
         flagStopIntegration = flagFligth;
     end        
        
     flagsArray = [flagFligth, flagAscent, flagBurning, flagAeroBrakes, flagPara1, flagPara2];
     flagMatr(n_old:n_old+n-1, :) = repmat(flagsArray, n, 1);
end

if settings.launchWindow
    fclose('all');
    delete('launchFlag.txt');
end

cpuTimes = cpuTimes(1:iTimes);

%% ASSEMBLE TOTAL FLIGHT STATE
Yf = Yf_tot(1:n_old, :);
Tf = Tf_tot(1:n_old, :);
flagMatr = flagMatr(1:n_old, :);

%% RETRIVE PARAMETERS FROM THE ODE

%if not(settings.electronics)
%    dataBallisticFlight = RecallOdeFcn(@ascent, Tf(flagMatr(:, 2)), Yf(flagMatr(:, 2), :), settings, C, uw, vw, ww, uncert);
%end

end


