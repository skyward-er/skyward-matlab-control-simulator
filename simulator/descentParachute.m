function dY = descentParachute(~, Y, settings, uw, vw, ww, para, uncert)
%{ 

descentParachute - ode function of the 3DOF Rigid Rocket-Paraachute Model

INPUTS:      
            - t, integration time;
            - Y, state vector, [ x y z | u v w ]:

                                * (x y z), NED{north, east, down} horizontal frame; 
                                * (u v w), horizontal frame velocities.

            - settings, rocket data structure;
            - uw, wind component along x;
            - vw, wind component along y;
            - ww, wind component along z;
            - uncert, wind uncertanties;
            - Hour, hour of the day of the needed simulation;
            - Day, day of the month of the needed simulation;

OUTPUTS:    
            - dY, state derivatives;

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

%% ADDING WIND (supposed to be added in NED axes);

if settings.wind.input 
    [uw, vw, ww] = wind_input_generator(settings, z, uncert);    
end

wind = [uw vw ww];

% Relative velocities (plus wind);
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

V_norm = norm([ur vr wr]);

%% CONSTANTS
S = settings.para(para).S;                                               % [m^2]   Surface
CD = settings.para(para).CD;                                             % [/] Parachute Drag Coefficient
CL = settings.para(para).CL;                                             % [/] Parachute Lift Coefficient
if para == 1
    pmass = 0 ;                                                          % [kg] detached mass
else
    pmass = sum(settings.para(1:para-1).mass) + settings.mnc;
end

g = 9.80655;                                                             % [N/kg] magnitude of the gravitational field at zero
m = settings.ms - pmass;                                                 % [kg] descend mass

%% ATMOSPHERE DATA

if -z < 0
    z = 0;
end

[~, ~, ~, rho] = atmosisa(-z+settings.z0);

%% REFERENCE FRAME
% The parachutes are approximated as rectangular surfaces with the normal
% vector perpendicular to the relative velocity

t_vect = [ur vr wr];                     % Tangenzial vector
h_vect = [-vr ur 0];                     % horizontal vector    

if all(abs(h_vect) < 1e-8)
    h_vect = [-vw uw 0];
end

t_vers = t_vect/norm(t_vect);            % Tangenzial versor
h_vers = -h_vect/norm(h_vect);           % horizontal versor

n_vect = cross(t_vers, h_vers);          % Normal vector
n_vers = n_vect/norm(n_vect);            % Normal versor

if (n_vers(3) > 0)                       % If the normal vector is downward directed
    n_vect = cross(h_vers, t_vers);
    n_vers = n_vect/norm(n_vect);
end

%% FORCES
D = 0.5*rho*V_norm^2*S*CD*t_vers';       % [N] Drag vector
L = 0.5*rho*V_norm^2*S*CL*n_vers';       % [N] Lift vector
Fg = m*g*[0 0 1]';                       % [N] Gravitational Force vector
F = -D + L + Fg;                         % [N] total forces vector

%% STATE DERIVATIVES
% velocity
du = F(1)/m;
dv = F(2)/m;
dw = F(3)/m;

%% FINAL DERIVATIVE STATE ASSEMBLING
dY(1:3) = [u v w]';
dY(4) = du;
dY(5) = dv;
dY(6) = dw;

dY = dY';

