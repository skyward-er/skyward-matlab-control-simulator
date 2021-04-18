function [Yf, Tf, t_ada, t_kalman, cpuTimes, flagMatr, dataBallisticFlight] = std_run_control(settings)
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

Author: Alessandro Del Duca
Skyward Experimental Rocketry | ELC-SCS Dept
email: alessandro.delduca@skywarder.eu
Revision date: 18/03/2021

%}

if not(settings.ballisticFligth) && settings.ascentOnly
   error('To simulate a landing with the parachutes, settings.ascentOnly must be false') 
end

%% STARTING CONDITIONS
% Attitude
Q0 = angle2quat(settings.PHI, settings.OMEGA, 0*pi/180, 'ZYX')';           % Attitude initial condition

% State
X0 = [0; 0; 0;];                                                           % Position initial condition
V0 = [0; 0; 0;];                                                           % Velocity initial condition  
W0 = [0; 0; 0;];                                                           % Angular speed initial condition

initialCond = [X0; V0; W0; Q0; settings.Ixxf; settings.Iyyf; settings.Izzf];    
Y0 = initialCond;

%% WIND GENERATION
if settings.wind.input                                                     % will be computed inside the integrations
    uw = 0; vw = 0; ww = 0;
else
    [uw,vw,ww,~] = wind_const_generator(settings.wind.AzMin, settings.wind.AzMax,...
        settings.wind.ElMin, settings.wind.ElMax, settings.wind.MagMin, settings.wind.MagMax);
    
    if ww ~= 0
        warning('Pay attention using vertical wind, there might be computational errors')
    end
    
end

if settings.wind.input && all(settings.wind.input_uncertainty ~= 0)
    signn = randi([1, 4]);                                                 % 4 sign cases
    unc = settings.wind.input_uncertainty;
    
    switch signn
        case 1
            % unc = unc;
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

%% KALMAN INIT 
addpath('../kalman');

%% SENSORS INIT
addpath('../sensors');
addpath('../sensors/data/MS580301BA01');
addpath('../simulationData');
[s, c] = initSensors;

%% CONTROL INIT
addpath('../control');
addpath('../control/Cd_rho_computation');
csett  = controlConfig;
%% MAGNETIC FIELD MODEL
hmax   =   6000;

%Use this lines if your MATLAB version is up to 2020
dy     =    decyear(settings.launchDate);
XYZ0   =    wrldmagm(0, settings.lat0, settings.lon0, dy, '2020');        % World magnetic map at h = 0
XYZh   =    wrldmagm(hmax, settings.lat0, settings.lon0, dy, '2020');     % World magnetic map at h = 6000

% %Use this next line if your MATLAB version is previous to 2020
% load('magn_field.mat');

magneticFieldApprox = @(zSlm) XYZ0 + (XYZh-XYZ0)./hmax.*zSlm;              % Magnetic field linear interpolation


%% INTEGRATION
% setting initial condition before control phase
dt          =       1/settings.frequencies.controlFrequency;               % Time step of the controller
t0          =       0;
t1          =       t0 + dt;
vz          =       1;
z           =       1;
nmax        =       10000;
mach        =       0;
x           =       0;
n_old       =       1;
Yf_tot      =       zeros(nmax, 16);
Tf_tot      =       zeros(nmax, 1);
C           =       zeros(nmax, 1);
cpuTimes    =       zeros(nmax,1);
iTimes      =       0;
c.ctr_start =      -1;
i           =       1;
settings.kalman.pn_prec  =       settings.ada.p_ref;
%% Flag initializations
flagStopIntegration     =   true;
flagAscent              =   false;
flagMatr                =   false(nmax, 6);

if settings.launchWindow
    launchWindow;
    pause(0.01);
    launchFlag = false;
    lastLaunchflag = true;
else
    launchFlag = true;
end

fprintf('START:\n\n\n');
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
    
    if vz(end) >= 0 && launchFlag
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
    [~, ~, p, ~] = atmosisa(-Yf(:,3)) ; 
 
 
    if settings.dataNoise
        [sp, c] = acquisition_Sys(sensorData, s, c);
    end
  
    if iTimes==1 && settings.Ada
        ada_prev  =   settings.ada.x0;
        Pada_prev =   settings.ada.P0;
    elseif iTimes ~= 1 && settings.Ada
        ada_prev  =   xp_ada_tot(end,:);
        Pada_prev =   P_ada(:,:,end);
    end
   
    if iTimes==1 && settings.Kalman
        x_prev    =  [X0; V0; Q0(2:4); Q0(1);0;0;0];
        P_prev    =   0.01*eye(12);
    elseif iTimes ~= 1 && settings.Kalman
        x_prev    =   x_est_tot(end,:);
        P_prev    =   P_c(:,:,end);
    end
%% ADA 
    if settings.Ada && settings.dataNoise
    [xp_ada, xv_ada, P_ada, settings.ada]   =  run_ADA(ada_prev, Pada_prev,                ...
                                                       sp.pn, sensorData.barometer.time,   ...
                                                       settings.ada);
                                                   
     xp_ada_tot(c.n_ada_old:c.n_ada_old + size(xp_ada(:,1),1) -1,:)  = xp_ada(1:end,:);
     xv_ada_tot(c.n_ada_old:c.n_ada_old + size(xv_ada(:,1),1)-1,:)  = xv_ada(1:end,:);
     t_ada_tot(c.n_ada_old:c.n_ada_old + size(xp_ada(:,1),1)-1)     = sensorData.barometer.time;              
     c.n_ada_old = c.n_ada_old + size(xp_ada,1); 
     
     end
%% Navigation system
    if settings.Kalman && settings.dataNoise

    [x_c,P_c,settings.kalman]   =  run_kalman(x_prev, P_prev, sp, xv_ada, settings.kalman, XYZ0*0.01);
    
     x_est_tot(c.n_est_old:c.n_est_old + size(x_c(:,1),1)-1,:)  = x_c(1:end,:);
     t_est_tot(c.n_est_old:c.n_est_old + size(x_c(:,1),1)-1)    = sensorData.accelerometer.time;              
     c.n_est_old = c.n_est_old + size(x_c,1); 
     
    end
%% Control algorithm
    if flagAeroBrakes && settings.Kalman && settings.control
         zc    =    exp_mean(-x_c(:,3),0.8);
         vzc   =    exp_mean(-x_c(:,6),0.8);
         vc    =    exp_mean(sqrt(x_c(:,4).^2+x_c(:,5).^2+x_c(:,6).^2),0.8);
         if c.ctr_start == -1
            c.ctr_start = 0.1*(n - 1);
         end
         %% selection of controler type
         switch csett.flagPID 
             case 1
             [alpha_degree, vz_setpoint, z_setpoint, pid, U_linear, Cdd, delta_S, csett] = control_PID    (zc, vzc, vc, csett);
             case 2
             [alpha_degree, vz_setpoint, z_setpoint, pid, U_linear, Cdd, delta_S, csett] = control_Lin    (zc, vzc, vc, csett);
             case 3
             [alpha_degree, vz_setpoint, z_setpoint, csett]                              = control_Servo  (zc, vzc,  csett);
         end
         x = extension_From_Angle(alpha_degree);
         i = i + 1; 
    elseif flagAeroBrakes && ~settings.Kalman && settings.control
         if c.ctr_start == -1
            c.ctr_start = 0.1*(n - 1);
         end
        switch csett.flagPID 
             case 1
             [alpha_degree, vz_setpoint, z_setpoint, pid,U_linear, Cdd, delta_S, csett] =   control_PID     (z, vz, sqrt(vxxx^2 + vyyy^2 + vz^2),  csett);
             case 2
             [alpha_degree, vz_setpoint, z_setpoint, pid,U_linear, Cdd, delta_S, csett] =   control_Lin     (z, vz, sqrt(vxxx^2 + vyyy^2 + vz^2),  csett);
             case 3
             [alpha_degree, vz_setpoint, z_setpoint, csett]                             =   control_Servo   (z, vz(end),  csett);
        end
    
         x = extension_From_Angle(alpha_degree);
         i = i + 1; 
    else 
        x = 0;
    end    
    
    if settings.control == true  && flagAeroBrakes == 1    
         % Save the values to plot them
         c.vz_tot(i)    =  vz;
         c.z_tot(i)     =  z;
         c.vz_setpoint_tot(i)  =  vz_setpoint;
         c.z_setpoint_tot(i)   =  z_setpoint;
         c.alpha_degree_tot(i) =  alpha_degree;
         if csett.flagPID ~= 3
             c.Cd_tot(i)    =  Cdd;
             c.pid_tot(i)   =  pid;
             c.U_lin_tot(i) =  U_linear;
             c.dS_tot(i)    =  delta_S;
         end
    end

    % vertical velocity and position
    if flagAscent || (not(flagAscent) && settings.ballisticFligth)
        Q    =   Yf(end, 10:13);
        vels =   quatrotate(quatconj(Q), Yf(:, 4:6)); 
        vz   = - vels(end,3);   % down
        vxxx =   vels(end,2);   % north
        vyyy =   vels(end,1);   % east
    else
        vz   = - Yf(end, 6);  
%         vx = Yf(end, 5); 
%         vy = Yf(end, 4); 
    end
    
    v_ned = quatrotate(quatconj(Yf(:, 10:13)), Yf(:, 4:6)); 

    z    = -Yf(end, 3);
    xxx  =  Yf(end, 2); 
    yyy  =  Yf(end, 1); 
    

    if lastFlagAscent && not(flagAscent)
        Y0 = [Yf(end, 1:3), vels(end,:), Yf(end, 7:end)];
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
    Yf_tot(n_old:n_old+n-1, :)   =  Yf(1:end, :);
    Tf_tot(n_old:n_old+n-1, 1)   =  Tf(1:end, 1);
    c.Yf_tot(n_old:n_old+n-1, :) =  Yf(1:end, :);
    c.Tf_tot(n_old:n_old+n-1, 1) =  Tf(1:end, 1);
    c.p_tot(n_old:n_old+n-1, 1)  =  p(1:end, 1);
    C(n_old:n_old+n-1) = x;
    c.v_ned_tot(n_old:n_old+n-1,:) = v_ned;  
    
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
    
     flagMatr(n_old:n_old+n-1, :) = repmat([flagFligth, flagAscent, flagBurning, flagAeroBrakes, flagPara1, flagPara2], n, 1);
end

if settings.launchWindow
    fclose('all');
    delete('launchFlag.txt');
end

cpuTimes = cpuTimes(1:iTimes);

%% ASSEMBLE TOTAL FLIGHT STATE
Yf = Yf_tot(1:n_old, :);
Tf = Tf_tot(1:n_old, :);

t_ada    = settings.ada.t_ada;
t_kalman = settings.kalman.t_kalman;
i_apo = find(Tf < 24.8);
i_apo = max(i_apo);
if settings.Kalman
i_apo_est = find(t_est_tot < Tf(i_apo));
i_apo_est = max(i_apo_est);
end
flagMatr = flagMatr(1:n_old, :);

%% SAVE THE VARIABLES FOR PLOT PURPOSE
% kalman state plot
if settings.Kalman
    c.x_est_tot    =  x_est_tot;
    c.t_est_tot    =  t_est_tot;
    c.i_apo        =  i_apo;
    c.i_apo_est    =  i_apo_est; 
end

% ada state for plot
if settings.Ada
    c.xp_ada_tot   =  xp_ada_tot;
    c.xv_ada_tot   =  xv_ada_tot;  
    c.t_ada_tot    =  t_ada_tot;
end

% control
if settings.control
    c.flagPID      =  csett.flagPID;
end

c.plot_ada     =  settings.Ada && false; 
c.plot_sensors =  settings.dataNoise && false; 
c.plot_kalman  =  settings.Kalman && false;
c.plot_control =  settings.control && true;

%% RETRIVE PARAMETERS FROM THE ODE

if not(settings.electronics)
    dataBallisticFlight = RecallOdeFcn(@ascent, Tf(flagMatr(:, 2)), Yf(flagMatr(:, 2), :), settings, C, uw, vw, ww, uncert, tLaunch);
end
if ~settings.electronics 
    plot_all(c, csett)
end
if true
save('../simulationData/Ground_truth.mat','sensorData');
if settings.dataNoise 
save('../simulationData/Sensors.mat','c');
end
end

end


