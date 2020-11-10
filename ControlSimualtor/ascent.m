function [dY, parout] = ascent(t, Y, settings, c, uw, vw, ww, uncert, Hour, Day)
%{ 

ASCENT - ode function of the 6DOF Rigid Rocket Model

INPUTS:      
            - t, integration time;
            - Y, state vector, [ x y z | u v w | p q r | q0 q1 q2 q3 | m | Ixx Iyy Izz | thetax thetay thetaz]:

                                * (x y z), NED{north, east, down} horizontal frame; 
                                * (u v w), body frame velocities;
                                * (p q r), body frame angular rates;
                                * (thetax thetay thetaz), body angles;
                                * m , total mass;
                                * (Ixx Iyy Izz), Inertias;
                                * (q0 q1 q2 q3), attitude unit quaternion.
 

            - settings, rocket data structure;
            - uw, wind component along x;
            - vw, wind component along y;
            - ww, wind component along z;
            - uncert, wind uncertanties;
            - Hour, hour of the day of the needed simulation;
            - Day, day of the month of the needed simulation;
            - OMEGA, launchpad azimuth angle;

OUTPUTS:    
            - dY, state derivatives;
            - parout, interesting fligth quantities structure (aerodyn coefficients, forces and so on..).


NOTE: To get the NED velocities the body-frame must be multiplied for the
conjugated of the current attitude quaternion
E.G.  quatrotate(quatconj(Y(:,10:13)),Y(:,4:6))


Author: Ruben Di Battista
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu
April 2014; Last revision: 31.XII.2014

Author: Francesco Colombi
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: francesco.colombi@skywarder.eu
Release date: 16/04/2016

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
email: adriano.filippo.inno@skywarder.eu
Release date: 13/01/2018

%}

% recalling the state
% x = Y(1);
% y = Y(2);
z = Y(3);
u = Y(4);
v = Y(5);
w = Y(6);
p = Y(7);
q = Y(8);
r = Y(9);
q0 = Y(10);
q1 = Y(11);
q2 = Y(12);
q3 = Y(13);
m = Y(14);
Ixx = Y(15);
Iyy = Y(16);
Izz = Y(17);

%% QUATERION ATTITUDE

Q = [q0 q1 q2 q3];
Q_conj = [q0 -q1 -q2 -q3];
normQ = norm(Q);

if abs(normQ-1) > 0.1
    Q = Q/normQ;
end


%% ADDING WIND (supposed to be added in NED axes);


if settings.wind.model
   
    if settings.stoch.N > 1
        [uw,vw,ww] = wind_matlab_generator(settings,z,t,Hour,Day);
    else
        [uw,vw,ww] = wind_matlab_generator(settings,z,t);
    end
    
elseif settings.wind.input

    [uw,vw,ww] = wind_input_generator(settings,z,uncert);    
end

wind = quatrotate(Q, [uw vw ww]);

% Relative velocities (plus wind);
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

% Body to Inertial velocities
Vels = quatrotate(Q_conj,[u v w]);
V_norm = norm([ur vr wr]);

%% ATMOSPHERE DATA

if -z < 0     % z is directed as the gravity vector
    z = 0;
end

[~, a, P, rho] = atmosisa(-z+settings.z0);
M = V_norm/a;
M_value = M;

%% CONSTANTS

S = settings.S;              % [m^2] cross surface
C = settings.C;              % [m]   caliber
CoeffsE = settings.CoeffsE;  % Empty Rocket Coefficients
CoeffsF = settings.CoeffsF;  % Full Rocket Coefficients
g = 9.80655;                 % [N/kg] module of gravitational field at zero
tb = settings.tb;            % [s]     Burning Time
mfr = settings.mfr;          % [kg/s]  Mass Flow Rate

OMEGA = settings.OMEGA;      % [rad] Elevation Angle in the launch pad

% inertias for full configuration (with all the propellant embarqued) obtained with CAD's
Ixxf = settings.Ixxf;        % [kg*m^2] Inertia to x-axis
Iyyf = settings.Iyyf;        % [kg*m^2] Inertia to y-axis
Izzf = settings.Izzf;        % [kg*m^2] Inertia to z-axis

% inertias for empty configuration (all the propellant consumed) obtained with CAD's
Ixxe = settings.Ixxe;        % [kg*m^2] Inertia to x-axis
Iyye = settings.Iyye;        % [kg*m^2] Inertia to y-axis
Izze = settings.Izze;        % [kg*m^2] Inertia to z-axis


%% TIME-DEPENDENTS VARIABLES

dI = 1/tb*([Ixxf Iyyf Izzf]'-[Ixxe Iyye Izze]');

if t<tb
    mdot = -mfr;
    Ixxdot = -dI(1);
    Iyydot = -dI(2);
    Izzdot = -dI(3);
    T = interp1(settings.motor.exp_time, settings.motor.exp_thrust, t);
    
else             % for t >= tb the fligth condition is the empty one(no interpolation needed)
    mdot = 0;
    Ixxdot = 0;
    Iyydot = 0;
    Izzdot = 0;
    T = 0;
end

%% AERODYNAMICS ANGLES

if not(ur < 1e-9 || V_norm < 1e-9)
    alpha = atan(wr/ur);
    beta = atan(vr/ur);             % beta = asin(vr/V_norm); is the classical notation, Datcom uses this one though. 
else
    alpha = 0;
    beta = 0;
end

alpha_value = alpha;
beta_value = beta;

%% DATCOM COEFFICIENTS

A_datcom = settings.Alphas*pi/180;
B_datcom = settings.Betas*pi/180;
H_datcom = settings.Altitudes;
M_datcom = settings.Machs;
C_datcom = settings.Controls;
%
[CA, angle0] = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CA, CoeffsE.CA, alpha, M, beta, -z, c, t);
CYB = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CYB, CoeffsE.CYB, alpha, M, beta, -z, c, t);
CY0 = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CY, CoeffsE.CY, alpha, M, beta, -z, c, t);
CNA = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CNA, CoeffsE.CNA, alpha, M, beta, -z, c, t);
CN0 = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CN, CoeffsE.CN, alpha, M, beta, -z, c, t);
Cl = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CLL, CoeffsE.CLL, alpha, M, beta, -z, c, t);
Clp = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CLLP, CoeffsE.CLLP, alpha, M, beta, -z, c, t);
Cma = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CMA, CoeffsE.CMA, alpha, M, beta, -z, c, t);
Cm0 = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CM, CoeffsE.CM, alpha, M, beta, -z, c, t);
Cmad = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CMAD, CoeffsE.CMAD, alpha, M, beta, -z, c, t);
Cmq = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CMQ, CoeffsE.CMQ, alpha, M, beta, -z, c, t);
Cnb = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CLNB, CoeffsE.CLNB, alpha, M, beta, -z, c, t);
Cn0 = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CLN, CoeffsE.CLN, alpha, M, beta, -z, c, t);
Cnr = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CLNR, CoeffsE.CLNR, alpha, M, beta, -z, c, t);
Cnp = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.CLNP, CoeffsE.CLNP, alpha, M, beta, -z, c, t);
XCP_value = InterpAero(settings, A_datcom, M_datcom, B_datcom, H_datcom, C_datcom, CoeffsF.X_C_P, CoeffsE.X_C_P, alpha, M, beta, -z, c, t);

alpha0 = angle0(1); beta0 = angle0(2);

CN = (CN0 + CNA*(alpha-alpha0));
CY = (CY0 + CYB*(beta-beta0));
Cm = (Cm0 + Cma*(alpha-alpha0));
Cn = (Cn0 + Cnb*(beta-beta0));

if -z < settings.lrampa*sin(OMEGA)      % No torque on the Launch
    
    Fg = m*g*sin(OMEGA);                % [N] force due to the gravity
    X = 0.5*rho*V_norm^2*S*CA;
    F = -Fg +T -X;
    du = F/m;
    
    dv = 0;
    dw = 0;
    dp = 0;
    dq = 0;
    dr = 0;
    
    alpha_value = NaN;
    beta_value = NaN;
    Y = 0;
    Z = 0;
    XCP_value = NaN;
    
    
    if T < Fg                           % No velocity untill T = Fg
        du = 0;
    end
    
else
    
    %% FORCES
    % first computed in the body-frame reference system
    
    qdyn = 0.5*rho*V_norm^2;        %[Pa] dynamics pressure
    qdynL_V = 0.5*rho*V_norm*S*C; 

    X = qdyn*S*CA;              %[N] x-body component of the aerodynamics force
    Y = qdyn*S*CY;            %[N] y-body component of the aerodynamics force
    Z = qdyn*S*CN;           %[N] z-body component of the aerodynamics force
    Fg = quatrotate(Q,[0 0 m*g])';  %[N] force due to the gravity in body frame
    
    F = Fg +[-X+T,+Y,-Z]';          %[N] total forces vector
    
    %% STATE DERIVATIVES
    
    % velocity
    du = F(1)/m-q*w+r*v;
    dv = F(2)/m-r*u+p*w;
    dw = F(3)/m-p*v+q*u;
    
    % Rotation
    dp = (Iyy-Izz)/Ixx*q*r + qdynL_V/Ixx*(V_norm*Cl+Clp*p*C/2)-Ixxdot*p/Ixx;
    dq = (Izz-Ixx)/Iyy*p*r + qdynL_V/Iyy*(V_norm*Cm + (Cmad+Cmq)*q*C/2)...
        -Iyydot*q/Iyy;
    dr = (Ixx-Iyy)/Izz*p*q + qdynL_V/Izz*(V_norm*Cn + (Cnr*r+Cnp*p)*C/2)...
        -Izzdot*r/Izz;
    
end
% Quaternion
OM = 1/2* [ 0 -p -q -r  ;
            p  0  r -q  ;
            q -r  0  p  ;
            r  q -p  0 ];

dQQ = OM*Q';

%% FINAL DERIVATIVE STATE ASSEMBLING

dY(1:3) = Vels;
dY(4) = du;
dY(5) = dv;
dY(6) = dw;
dY(7) = dp;
dY(8) = dq;
dY(9) = dr;
dY(10:13) = dQQ;
dY(14) = mdot;
dY(15) = Ixxdot;
dY(16) = Iyydot;
dY(17) = Izzdot;
dY(18:20) = [p q r];
dY = dY';

%% SAVING QUANTITIES FOR PLOTS 

parout.integration.t = t;

parout.interp.M = M_value;
parout.interp.alpha = alpha_value;
parout.interp.beta = beta_value;
parout.interp.alt = -z;

parout.wind.NED_wind = [uw, vw, ww];
parout.wind.body_wind = wind;

parout.velocities=Vels;

parout.forces.AeroDyn_Forces = [X, Y, Z];
parout.forces.T = T;

parout.air.rho = rho;
parout.air.P = P;

parout.accelerations.body_acc = [du, dv, dw];
parout.accelerations.ang_acc = [dp, dq, dr];

parout.coeff.CA = CA;
parout.coeff.CYB = CYB;
parout.coeff.CNA = CNA;
parout.coeff.Cl = Cl;
parout.coeff.Clp = Clp;
parout.coeff.Cma = Cma;
parout.coeff.Cmad = Cmad;
parout.coeff.Cmq = Cmq;
parout.coeff.Cnb = Cnb;
parout.coeff.Cnr = Cnr;
parout.coeff.Cnp = Cnp;
parout.coeff.XCP = XCP_value;



