function [Yc,Tc] = std_run_control(settings)
%{ 

STD_RUN_BALLISTIC - This function runs a standard ballistic (non-stochastic) simulation

INTPUTS: 
            - settings, rocket data structure;

OUTPUTS:
            - Tf, Total integration time vector; 
            - Yf, Total State Matrix;
            - Ta, Ascent Integration time vector; 
            - Ya, Ascent State Matrix;
            - bound_value, Usefull values for the plots.

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

[Ta, Ya] = ode113(@ascent, [0, tf], Y0a, settings.ode.optionsasc1, settings, uw, vw, ww, uncert);

%% CONTROL PHASE 
dt = settings.freq;
t0 = Ta(end);
t1 = t0 + dt;
vz = 1;
Y0 = Ya(end,:);
% pressione

while vz > -10 
   
% controllo

%[A] = controllo(Y0,t0);   % area totale aerofreno esposto
 
% dynamics 
[Tc,Yc] = ode45(@ascent, [t0, t1], Y0, [], settings, uw, vw, ww, uncert);

Q = Yc(end,10:13);
vels = quatrotate(quatconj(Q),Yc(end,4:6));
vz = - vels(3);
% aggiungere pressione 

t0 = t0 + dt;
t1 = t1 + dt;
Y0 = Yc(end,:);


end


