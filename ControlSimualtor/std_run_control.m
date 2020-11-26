function [Yf, Tf, data_flight] = std_run_control(settings)
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

if settings.wind.model && settings.wind.input
    error('Both wind model and input wind are true, select just one of them')
end

if settings.wind.HourMin ~= settings.wind.HourMax || settings.wind.HourMin ~= settings.wind.HourMax
    error('In standard simulations with the wind model the day and the hour of launch must be unique, check config.m')
end

%% STARTING CONDITIONS
% Attitude

if settings.upwind
    settings.PHI = mod(Azw + pi, 2*pi);
end

% Attitude
Q0 = angle2quat(settings.PHI, settings.OMEGA, 0*pi/180, 'ZYX')';

% State
X0 = [0 0 0]';
V0 = [0 0 0]';
W0 = [0 0 0]';
theta0 = [0 0 0]';
Y0a = [X0; V0; W0; Q0; settings.m0; settings.Ixxf; settings.Iyyf; settings.Izzf; theta0];

%% WIND GENERATION

if settings.wind.model || settings.wind.input   % will be computed inside the integrations
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

tf = settings.ode.final_time;

%% BURNING ASCENT
c = 0;
[Ta, Ya] = ode113(@ascent, [0, tf], Y0a, settings.ode.optionsasc1, settings, c, uw, vw, ww, uncert);

%% CONTROL PHASE

% setting initial condition before control phase
dt = settings.freq;
t0 = Ta(end);
t1 = t0 + dt;
vz = 1;
Y0 = Ya(end,:);
% [~, ~, p0, ~] = atmosisa(Ya(end,3));
nmax = 10000;
Yc_tot = zeros(nmax, 20);
Tc_tot = zeros(nmax, 1);
C = zeros(nmax, 1);
n_old = 1;


alphaDefault=settings.all_alpha(1);
alpha=alphaDefault;
alphaOld=alphaDefault;

z=-Y0(3);
zLast=z;

start1state=230; % velocity when control should start
start2state=100; % velocity when control should end and alpha keeps last control value

% control phase dynamics integration
while vz > -10 || n_old < nmax
    
    %% velocity and altitude
    z=-Y0(3);
    Vz=(z-zLast)/dt; % Caluclated vertical velocity
    Vu=Y0(4);
    zLast=z;

    %% Control
    if Vu<start1state && Vu>start2state % start control when velocity under 233 m/s
        %% Controller
        % get index of actual hight z
        if z>=settings.hightInterval(end)
            iz = length(settings.hightInterval);
        elseif z<=settings.hightInterval(1)
            iz = 1;
        else
            iz = interp1(settings.hightInterval,1:length(settings.hightInterval),z,'nearest');
        end

        % linear interpolate between the values for all_alpha at the actual hight
        % index
        if Vz>=settings.all_Vz(iz,end)
            alpha = settings.all_alpha(end);
        elseif Vz<=settings.all_Vz(iz,1) 
            alpha = settings.all_alpha(1);
        else
            alpha = interp1(settings.all_Vz(iz,:),settings.all_alpha,Vz);
        end
    % do not change alpha at end of the flight (not controllable)
    elseif Vu<=start2state
        alpha=alphaOld;
    else
        alpha = settings.all_alpha(1);
    end
    
    %% Limit changing speed of alpha
    dAlphaMax = (pi/180*60)/0.13*dt; % 60 degree/0.13 sec from servo
    dAlphaMax = dAlphaMax*0.5; % multiplied by a value < 1
    difAlphas=alpha-alphaOld;
    difAlphas=max(difAlphas,-dAlphaMax);
    difAlphas=min(difAlphas,dAlphaMax);
    alpha=alphaOld+difAlphas;
    alphaOld=alpha;
    
    
    % PLOT SERVO CONTROL ANGLE
    figure(10);
    plot(t0, alpha,'*'),grid on;
    hold on;
    xlabel('Time [s]'); ylabel('Control Angle [rad]')
    
    
    
    %% Calculate area of airbreaks
    % alpha varies from 0 to 0.89 rad
    alpha=max(alpha,0);
    alpha=min(alpha,0.89);
    
    % calculate air brake area
    Aall=(-9.43386*alpha^2 + 19.86779*alpha)*10^(-3); % area of all airbrakes (3)
    A=Aall/3; % area of one airbrake
    
%     [At] = controllo(Y0,t0);           % total aerobrakes wet Area
%     A = At/3;                         % single aerobrake wet Area
%      A = settings.Atot/3;                % waiting for the control
    
    c = A/settings.brakesWidth;        % approximated aerobrakes heigth --> control variable of the simulator
    
    % dynamics
    [Tc,Yc] = ode45(@ascent, [t0, t1], Y0, [], settings, c, uw, vw, ww, uncert);
    
    % evaluate the condition for cycle condition 
    Q = Yc(end,10:13);
    vels = quatrotate(quatconj(Q),Yc(end,4:6));
    vz = - vels(3);
    
    % update ode every cycle 
    t0 = t0 + dt;
    t1 = t1 + dt;
    Y0 = Yc(end,:);
%     [T, ~, p0, ~] = atmosisa(Yc(end, 3));        % pressure and temperature at each sample time  
            
    % assemble total state
    [n, ~] = size(Yc);
    Yc_tot(n_old:n_old+n-1,:) = Yc(1:end,:);
    Tc_tot(n_old:n_old+n-1) = Tc(1:end,1);
    C(n_old:n_old+n-1) = c;
    
    n_old = n_old + n -1;
   
end


%% ASSEMBLE TOTAL FLIGHT STATE

Yc_tot = Yc_tot(1:n_old,:);
Tc_tot = Tc_tot(1:n_old,:);

Yf = [Ya; Yc_tot(2:end,:)];
Tf = [Ta; Tc_tot(2:end,:)];

C = [zeros(length(Ta) - 1, 1); C];

%% RETRIVE PARAMETERS FROM THE ODE
data_flight = RecallOdeFcn(@ascent, Tf, Yf, settings, C, uw, vw, ww, uncert);


