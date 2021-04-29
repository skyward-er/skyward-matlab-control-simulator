function acc = accelerometerAscent(t, Y, settings, c, uw, vw, ww, uncert)
%{ 

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
email: adriano.filippo.inno@skywarder.eu
Release date: 30/11/2020

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

%% QUATERION ATTITUDE
Q = [q0 q1 q2 q3];
normQ = norm(Q);

if abs(normQ - 1) > 0.1
    Q = Q/normQ;
end

%% ADDING WIND (supposed to be added in NED axes);
if settings.wind.input
    [uw, vw, ww] = wind_input_generator(settings, z, uncert);    
end

dcm = quatToDcm(Q);
wind = dcm*[uw; vw; ww];

% Relative velocities (plus wind);
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

% Body to Inertial velocities
V_norm = norm([ur vr wr]);

%% ATMOSPHERE DATA
if -z < 0     % z is directed as the gravity vector
    z = 0;
end

h = -z + settings.z0;

atmC = [9.80665, 1.4, 287.0531, 0.0065, 11000, 20000, ...
            1.225, 101325, 288.15]; % atmosisa constants:

h(h > atmC(6)) = atmC(6);
h(h < 0) = 0;
hGrThHTS = (h > atmC(5));

h_tmp = h;
h_tmp(hGrThHTS) = atmC(5);

T = atmC(9) - atmC(4)*h_tmp;

expon = ones(size(h));
expon(hGrThHTS) = exp(atmC(1)./(atmC(3)*T(hGrThHTS)).*(atmC(5) - h(hGrThHTS)));

a = sqrt(T*atmC(2)*atmC(3));

theta = T/atmC(9);

rho = atmC(7)*theta.^((atmC(1)/(atmC(4)*atmC(3)))-1.0).*expon;  

M = V_norm/a;

%% CONSTANTS
S = settings.S;              % [m^2] cross surface
CoeffsE = settings.CoeffsE;  % Empty Rocket Coefficients
CoeffsF = settings.CoeffsF;  % Full Rocket Coefficients
g = 9.80655;                 % [N/kg] module of gravitational field at zero
tb = settings.tb;            % [s]     Burning Time
OMEGA = settings.OMEGA;      % [rad] Elevation Angle in the launch pad

%% TIME-DEPENDENTS VARIABLES

if t < tb
    m = settings.ms + interp1(settings.motor.exp_time, settings.motor.exp_m, t);
    T = interp1(settings.motor.exp_time, settings.motor.exp_thrust, t);
else 
    m = settings.ms;
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

%% DATCOM COEFFICIENTS
A_datcom = settings.Alphas*pi/180;
B_datcom = settings.Betas*pi/180;
H_datcom = settings.Altitudes;
M_datcom = settings.Machs;
C_datcom = settings.Controls;

cellT = {A_datcom, M_datcom, B_datcom, H_datcom};
inst = [alpha, M, beta, -z];

index = zeros(4, 1);
for i = 1:4
    [~, index(i)] = min(abs(cellT{i} - inst(i)));
end


coeffsNames = {'CA', 'CYB', 'CY', 'CNA', 'CN'};
coeffsValues = zeros(5, 1);

for i = 1:5
    
    CmatE = CoeffsE.(coeffsNames{i});
    CmatF = CoeffsF.(coeffsNames{i});
    
    if c == 0
        VE = CmatE(index(1), index(2), index(3), index(4), 1);
    else
        c_cmp = C_datcom(c > C_datcom); 
        n0 = length(c_cmp);
        n1 = n0 + 1;
        c0 = c_cmp(end);
        c1 = C_datcom(n1);
        C0 =  CmatE(index(1), index(2), index(3), index(4), n0);
        C1 =  CmatE(index(1), index(2), index(3), index(4), n1);
        VE = C1 + ((C1 - C0)/(c1 - c0))*(c - c1);
    end

    if t <= tb
        VF = CmatF(index(1), index(2), index(3), index(4));
        
        coeffsValues(i) =  t/tb*(VE-VF)+VF;
    else 
        coeffsValues(i) = VE;
    end

end

% Retrieve Coefficients 
CA = coeffsValues(1); CYB = coeffsValues(2); CY0 = coeffsValues(3);
CNA = coeffsValues(4); CN0 = coeffsValues(5);

alpha0 = A_datcom(index(1)); beta0 = B_datcom(index(3));

CN = (CN0 + CNA*(alpha-alpha0));
CY = (CY0 + CYB*(beta-beta0));

if -z < settings.lrampa*sin(OMEGA)      % No torque on the Launch
    
    Fg = m*g*sin(OMEGA);                % [N] force due to the gravity
    X = 0.5*rho*V_norm^2*S*CA;
    F = -Fg +T -X;
    du = F/m;
    
    dv = 0;
    dw = 0;
    
    if T < Fg                           % No velocity untill T = Fg
        du = 0;
    end
    
else
    
    %% FORCES
    % first computed in the body-frame reference system
    qdyn = 0.5*rho*V_norm^2;        %[Pa] dynamics pressure

    X = qdyn*S*CA;              %[N] x-body component of the aerodynamics force
    Y = qdyn*S*CY;            %[N] y-body component of the aerodynamics force
    Z = qdyn*S*CN;           %[N] z-body component of the aerodynamics force
%    Fg = quatrotate(Q,[0 0 m*g])';  %[N] force due to the gravity in body frame
    
    F = [-X+T,+Y,-Z]';          %[N] total forces vector
    
    %% STATE DERIVATIVES
    % velocity
    du = F(1)/m;
    dv = F(2)/m;
    dw = F(3)/m;
    
    
end

acc = [du, dv, dw];

