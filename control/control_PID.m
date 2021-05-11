function [alpha_degree, Vz_setpoint, z_setpoint, pid, U_linear, Cd, delta_S, csett] = control_PID(z, Vz, V_mod, csett)

% Author: Leonardo Bertelli
% Co-Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@kywarder.eu
% email: leonardo.bertelli@skywarder.eu, alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

%{
CONTROL_ALGORITHM  Finds trejectory (z-Vz) to follow and uses a PI
controler to follow the trejectory and then transfere it with a to a force


  INPUTS:
  z               acutal hight of the rocket
  Vz              actual vertical velocity of the rocket
  V_mod           actual roket velocity in the direction of the main axis
  sample_time     sample time of the control system

  OUTPUTS:
  alpha_degree    output angle for servo
  Vz_setpoint     setpoint vertical velocity from trejectory
  z_setpoint      setpoint hight
  pid             PI control output
  U_linear        linearized PI controler output
  Cd              resulting drag coefficiant 
  delta_S         resulting force
%}

%% CHOOSE THE SETPOINTS

[z_setpoint, Vz_setpoint, csett] = set_Trajectory(z, Vz, csett);

%% PI ALGORITHM 

% If e>0 the rocket is too fast. I slow it down with Fx>0 --> open aerobrakes
% If e<0 the rocket is too slow. I speed it up with Fx<0 --> close aerobrakes

% Parameters
ro = getRho(z);

% Control variable limits
Umin = 0;     
Umax = 0.5*ro*csett.S0*Vz*V_mod; % Cd limit check

% Input for PI controler
error = (Vz - Vz_setpoint); % > 0 (in teoria)

% P part of controler
P = csett.Kp_1*error;

% I part of controler
if csett.saturation == false
    csett.I = csett.I + csett.Ki_1*error;
end

% Combining PI controler
U = P + csett.I;

% Anti-windup
[U, csett.saturation] = Saturation(U, Umin, Umax);

%% TRANSFORMATION FROM U to delta_S 

% Get the Cd for each possible aerobrake surface
Cd_available = 1:length(csett.delta_S_available);
for ind = 1:length(csett.delta_S_available)
    Cd_available(ind) = getDrag(V_mod,z,csett.delta_S_available(ind), csett.coeff_Cd);
end
Cd_available = Cd_available';

% For all possible delta_S compute Fdrag
% Then choose the delta_S which gives an Fdrag which has the minimum error if compared with F_drag_pid
[~, index_minimum] = min( abs(U - 0.5*ro*csett.S0*Cd_available*Vz*V_mod) );
delta_S = csett.delta_S_available(index_minimum); 

% Just for plotting
pid = U;
U_linear = 0.5*ro*csett.S0*Cd_available(index_minimum)*Vz*V_mod;
Cd = Cd_available(index_minimum);

%% TRANSFORMATION FROM delta_S to SERVOMOTOR ANGLE DEGREES 

alpha_rad      = (-csett.b + sqrt(csett.b^2 + 4*csett.a*delta_S)) / (2*csett.a);

% Alpha saturation 
[alpha_rad, ~] = Saturation(alpha_rad, 0, 0.89);

alpha_degree   = (alpha_rad*180)/pi;

%% LIMIT THE RATE OF THE CONTROL VARIABLE 

alpha_degree = rate_Limiter(alpha_degree, csett.alpha_degree_prec, csett.rate_limiter, csett.sample_time);

alpha_degree = smooth_Control(alpha_degree, csett.alpha_degree_prec, csett.filter_coeff);
csett.alpha_degree_prec = alpha_degree;

end