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
theta0 = [0 0 0]';
Y0 = [X0; V0; W0; Q0; settings.m0; settings.Ixxf; settings.Iyyf; settings.Izzf; theta0];

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
Yf_tot = zeros(nmax, 20);
Tf_tot = zeros(nmax, 1);
C = zeros(nmax, 1);
n_old = 1;
cpuTimes = zeros(nmax,1);
iTimes = 0;



%%%%%%%%%%%%%%%%%%%%% VARIABLES NEEDED FOR CONTROL %%%%%%%%%%%%%%%%%%%%%%%%

% Define global variables
global coeff_Cd data_1 data_2 data_3 data_4 data_5 data_6 data_7 data_8 data_9 data_10 data_11

% Load coefficients for Cd
data = load('coeffs.mat');
coeff_Cd = data.coeffs;

% Load the trajectories
data_1 = load('Trajectory_1.mat');
data_2 = load('Trajectory_2.mat');
data_3 = load('Trajectory_3.mat');
data_4 = load('Trajectory_4.mat');
data_5 = load('Trajectory_5.mat');
data_6 = load('Trajectory_6.mat');
data_7 = load('Trajectory_7.mat');
data_8 = load('Trajectory_8.mat');
data_9 = load('Trajectory_9.mat');
data_10 = load('Trajectory_10.mat');
data_11 = load('Trajectory_11.mat');

% Define global variables
global Kp Ki I alpha_degree_prec iteration_flag chosen_trajectory saturation
Kp = 40; % 40
Ki = 8; % 8
I = 0;
alpha_degree_prec = 0;
iteration_flag = 1;
saturation = false;

index_plot = 1; % To plot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
            Yf = [Yd, zeros(nd, 7), settings.m0*ones(nd, 1), settings.Ixxe*ones(nd, 1), ...
                settings.Iyye*ones(nd, 1), settings.Iyye*ones(nd, 1), zeros(nd, 3)];
        end
    end

    [sensorData] = manageSignalFrequencies(magneticFieldApprox, flagAscent, settings, Yf, Tf, x, uw, vw, ww, uncert);
    
    if settings.dataNoise
        Yf = acquisitionSystem(Yf);    
    end
  
    %%%%%%% kalmann filter %%%%%%%%
    % kalman(p, acc_body, ang_vel, q, [u, v, w]ned, [x, y, z]ned )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flagAeroBrakes
         alpha_degree = controlAlgorithm(z, vz, normV, dt);
         x = get_extension_from_angle(alpha_degree);
         
         % Save the values of the control variable to plot it
         plot_control_variable(index_plot) = alpha_degree;
         index_plot = index_plot + 1;
    else 
        x = 0;
    end    

    % vertical velocity and position
    if flagAscent || (not(flagAscent) && settings.ballisticFligth)
        Q = Yf(end, 10:13);
        vels = quatrotate(quatconj(Q), Yf(end, 4:6));
        vz = - vels(3);
    else
        vz = -Yf(end, 6);
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
    Tf_tot(n_old:n_old+n-1) = Tf(1:end, 1);
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
flagMatr = flagMatr(1:n_old, :);

%% RETRIVE PARAMETERS FROM THE ODE (commentato senò non stampava grafici)
% if not(settings.electronics)
%     dataBallisticFlight = RecallOdeFcn(@ascent, Tf(flagMatr(:, 2)), Yf(flagMatr(:, 2), :), settings, C, uw, vw, ww, uncert);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Obtain the control variable
time = 0:dt:(length(plot_control_variable)-1)*dt;
% length(Ya)
% size(Ya)
% size(Tf)
% length(Yf)
% plot_control_variable = [zeros(length(Ya),1);
%                          control_variable(:,1)];
% size(plot_control_variable)                  
                     
% Obtain the altitude
plot_z = -Yf(:,3);

% Obtain the vertical velocity
nStates = length(Yf);
plot_Vz = zeros(nStates, 1);
for index = 1:nStates
    Q = Yf(index,10:13);
    vels = quatrotate(quatconj(Q), Yf(index,4:6));
    plot_Vz(index) = - vels(3);
end

% Control variable: aerobrake_surface_total
figure('Name','Total aerobrake surface after burning phase','NumberTitle','off');
plot(time, plot_control_variable), grid on;
% plot(Tf, plot_control_variable), grid on;
axis([0,20, 0,60])
xlabel('time [s]'), ylabel('A [m^2]');

% Altitude
figure('Name','Time, Altitude','NumberTitle','off');
plot(Tf, plot_z), grid on;
xlabel('time [s]'), ylabel('z [m]');


% Vertical Velocity
figure('Name','Time, Vertical Velocity','NumberTitle','off');
plot(Tf, plot_Vz), grid on;
xlabel('time [s]'), ylabel('Vz [m/s]');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

