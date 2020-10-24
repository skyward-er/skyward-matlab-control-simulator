function [dY] = ascent(t, Y, settings, uw, vw, ww, uncert, Hour, Day, OMEGA)
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
x = Y(1);
y = Y(2);
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


[lat, lon, ~] = ned2geodetic(x, y, 0, settings.lat0, settings.lon0, 0, wgs84Ellipsoid);     % geographic coordinates

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

%% INTERPOLATION AT THE BOUNDARIES

if M > M_datcom(end)
    
    M = M_datcom(end);
    
end

if M < M_datcom(1)
    
    M = M_datcom(1);
    
end

if alpha > A_datcom(end)
    
    alpha = A_datcom(end);
    
elseif alpha < A_datcom(1)
    
    alpha = A_datcom(1);
    
end

if beta > B_datcom(end)
    
    beta = B_datcom(end);
    
elseif beta < B_datcom(1)
    
    beta = B_datcom(1);
end

if -z > H_datcom(end)
    
    z = -H_datcom(end);
    
elseif -z < H_datcom(1)
    
    z = -H_datcom(1);
    
end

%% CHOSING THE FULL CONDITION VALUE
% interpolation of the coefficients with the value in the nearest condition of the Coeffs matrix

CAf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CA,alpha,M,beta,-z);
CYBf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CYB,alpha,M,beta,-z);
CNAf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CNA,alpha,M,beta,-z);
Clf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CLL,alpha,M,beta,-z);
Clpf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CLLP,alpha,M,beta,-z);
Cmaf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CMA,alpha,M,beta,-z);
Cmadf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CMAD,alpha,M,beta,-z);
Cmqf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CMQ,alpha,M,beta,-z);
Cnbf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CLNB,alpha,M,beta,-z);
Cnrf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CLNR,alpha,M,beta,-z);
Cnpf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.CLNP,alpha,M,beta,-z);
XCPf = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsF.X_C_P,alpha,M,beta,-z);

%% CHOSING THE EMPTY CONDITION VALUE
% interpolation of the coefficients with the value in the nearest condition of the Coeffs matrix

if t >= tb && M <= settings.Mach_control

    c = 1;
    
else
    c = 1;
end
    
CAe = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CA(:, :, :, :, c) ,alpha,M,beta,-z);
CYBe = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CYB(:, :, :, :, c) ,alpha,M,beta,-z);
CNAe = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CNA(:, :, :, :, c) ,alpha,M,beta,-z);
Cle = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CLL(:, :, :, :, c) ,alpha,M,beta,-z);
Clpe = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CLLP(:, :, :, :, c) ,alpha,M,beta,-z);
Cmae = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CMA(:, :, :, :, c) ,alpha,M,beta,-z);
Cmade = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CMAD(:, :, :, :, c) ,alpha,M,beta,-z);
Cmqe = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CMQ(:, :, :, :, c) ,alpha,M,beta,-z);
Cnbe = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CLNB(:, :, :, :, c) ,alpha,M,beta,-z);
Cnre = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CLNR(:, :, :, :, c) ,alpha,M,beta,-z);
Cnpe = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.CLNP(:, :, :, :, c) ,alpha,M,beta,-z);
XCPe = interp4_easy(A_datcom,M_datcom,B_datcom,H_datcom,CoeffsE.X_C_P(:, :, :, :, c) ,alpha,M,beta,-z);

%% LINEAR INTERPOLATION BETWEEN THE TWO CONDITIONS
% Computing the value of the aerodynamics coefficients at a certain time
% Needed only for t<tb because for t>=tb the condition is the empty one

if t < tb
    CA = t/tb*(CAe-CAf)+CAf;
    CYB = t/tb*(CYBe-CYBf)+CYBf;
    CNA = t/tb*(CNAe-CNAf)+CNAf;
    Cl = t/tb*(Cle-Clf)+Clf;
    Clp = t/tb*(Clpe-Clpf)+Clpf;
    Cma = t/tb*(Cmae-Cmaf)+Cmaf;
    Cmad = t/tb*(Cmade-Cmadf)+Cmadf;
    Cmq = t/tb*(Cmqe-Cmqf)+Cmqf;
    Cnb = t/tb*(Cnbe-Cnbf)+Cnbf;
    Cnr = t/tb*(Cnre-Cnrf)+Cnrf;
    Cnp = t/tb*(Cnpe-Cnpf)+Cnpf;
    XCP_value = t/tb*(XCPe-XCPf)+XCPf;
else
    CA = CAe;
    CYB = CYBe;
    CNA = CNAe;
    Cl = Cle;
    Clp = Clpe;
    Cma = Cmae;
    Cmad =Cmade;
    Cmq = Cmqe;
    Cnb = Cnbe;
    Cnr = Cnre;
    Cnp = Cnpe;
    XCP_value = XCPe;
end

if -z < settings.lrampa*sin(OMEGA)      % No torque on the Launch
    
    Fg = m*g*sin(OMEGA);                % [N] force due to the gravity
    X = 1.4*0.5*rho*V_norm^2*S*CA;
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
    
    if c == 1
        X = 1.4*qdyn*S*CA;              %[N] x-body component of the aerodynamics force
    else
        X = 1.2*qdyn*S*CA;              %[N] x-body component of the aerodynamics force
    end
        
    Y = qdyn*S*CYB*beta;            %[N] y-body component of the aerodynamics force
    Z = qdyn*S*CNA*alpha;           %[N] z-body component of the aerodynamics force
    Fg = quatrotate(Q,[0 0 m*g])';  %[N] force due to the gravity in body frame
    
    F = Fg +[-X+T,+Y,-Z]';          %[N] total forces vector
    
    %% STATE DERIVATIVES
    
    % velocity
    du = F(1)/m-q*w+r*v;
    dv = F(2)/m-r*u+p*w;
    dw = F(3)/m-p*v+q*u;
    
    % Rotation
    dp = (Iyy-Izz)/Ixx*q*r + qdynL_V/Ixx*(V_norm*Cl+Clp*p*C/2)-Ixxdot*p/Ixx;
    dq = (Izz-Ixx)/Iyy*p*r + qdynL_V/Iyy*(V_norm*Cma*alpha + (Cmad+Cmq)*q*C/2)...
        -Iyydot*q/Iyy;
    dr = (Ixx-Iyy)/Izz*p*q + qdynL_V/Izz*(V_norm*Cnb*beta + (Cnr*r+Cnp*p)*C/2)...
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
