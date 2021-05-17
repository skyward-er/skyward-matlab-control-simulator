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
  z               actual hight of the rocket
  Vz              actual vertical velocity of the rocket
  V_mod           modulus of the rocket velocity
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
S0 = csett.S0;

cdbar = getDrag(V_mod, z, csett.Sbar, csett.coeff_Cd);
Ubar = 0.5*ro*cdbar*S0*Vz*V_mod;

% Control variable limits
camin = getDrag(V_mod, z, csett.delta_S_available(1), csett.coeff_Cd);
Umin  = 0.5*ro*camin*S0*Vz*V_mod;     
dUmin = Umin - Ubar;

camax = getDrag(V_mod, z, csett.delta_S_available(end), csett.coeff_Cd);
Umax  = 0.5*ro*camax*S0*Vz*V_mod; % Cd limit check
dUmax = Umax - Ubar;



% Input for PI controler
error = (Vz - Vz_setpoint); % > 0 (in teoria)

% P part of controler
P = csett.Kp_1*error;

% I part of controler
if csett.saturation == false
    csett.I = csett.I + csett.Ki_1*error;
end

% Combining PI controler
dU = P + csett.I;

% Anti-windup
[dU, csett.saturation] = Saturation(dU, dUmin, dUmax);

%% TRANSFORMATION FROM U to delta_S 

% Get the Cd for each possible aerobrake surface
ca_available = zeros(1,length(csett.delta_S_available));
dU_available  = zeros(1,length(csett.delta_S_available));
for ind = 1:length(csett.delta_S_available)
    ca_available(ind) = getDrag(V_mod, z, csett.delta_S_available(ind), csett.coeff_Cd)';
    dU_available(ind)  = 0.5*ro*ca_available(ind)*S0*Vz*V_mod - Ubar;
end

% For all possible delta_S compute Fdrag
% Then choose the delta_S which gives an Fdrag which has the minimum error if compared with F_drag_pid
[~, index_minimum] = min( abs(dU - dU_available ));
delta_S = csett.delta_S_available(index_minimum); 

% Just for plotting
pid = dU + Ubar;
U_linear = 0.5*ro*ca_available(index_minimum)*S0*Vz*V_mod;
Cd = ca_available(index_minimum);

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