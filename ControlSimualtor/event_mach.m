function [value, isterminal, direction] = event_mach(t, Y, settings, uw, vw, ww, varargin)
%{

EVENT_MACH - Event function to stop simulation at certain mach checking when a value tends to zero;
               the value taken is to account is the mach tends to certain
               value setted in the config.

INPUTS:     
            - t, integration time;
            - Y, state vector, [ x y z | u v w | p q r | q0 q1 q2 q3 | m | Ixx Iyy Izz ]:

                                * (x y z), NED{north, east, down} horizontal frame; 
                                * (u v w), body frame velocities;
                                * (p q r), body frame angular rates;
                                * m , total mass;
                                * (Ixx Iyy Izz), Inertias;
                                * (q0 q1 q2 q3), attitude unit quaternion.

            - settings, rocket data structure.

OUTPUTS:        
            - isterminal, logical input to stop the integration;
            - direction, to select the sign that the function must have when tends to zero, 1 = positive;
            - value, selected value to check if the integration has to be stopped (vertical velocity).

Author: Matteo Pozzoli
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: matteo.pozzoli@skywarder.eu

%}

Q = Y(10:13)';
z = Y(3);
u = Y(4);
v = Y(5);
w = Y(6);


if settings.wind.model
    
    [uw,vw,ww] = wind_matlab_generator(settings,z,t);
        
elseif settings.wind.input
    
    [uw,vw,ww] = wind_input_generator(settings,z,uncert);
    
end

wind = quatrotate(Q, [uw vw ww]);

% Relative velocities (plus wind);
ur = u - wind(1);
vr = v - wind(2);
wr = w - wind(3);

% Body to Inertial velocities
V_norm = norm([ur vr wr]);

[~, a, ~, ~] = atmosisa(-z+settings.z0);
M = V_norm/a;


% Stop checking if I'm in Propulsion Phase
if t > settings.tb
    value =  M - settings.Mach_control;
else
    value = 1;
end

isterminal = 1;
direction = 0;


end
