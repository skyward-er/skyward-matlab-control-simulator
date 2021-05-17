function [alpha_degree, Vz_setpoint, z_setpoint, U, formula, Cd, delta_S, csett] = control_Lin(z, Vz, V_mod, csett)

% Author: Leonardo Bertelli
% Co-Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@kywarder.eu
% email: leonardo.bertelli@skywarder.eu, alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

%{
CONTROL_ALGORITHM_LINEARIZED  Finds trejectory (z-Vz) to follow and uses a
PI controler to follow the trejectory but uses linearization to calculate delta_S

  INPUTS:
  z               acutal hight of the rocket
  Vz              actual vertical velocity of the rocket
  V_mod           actual roket velocity in the direction of the main axis
  sample_time     sample time of the control system

  OUTPUTS:
  alpha_degree    output angle for servo
  Vz_setpoint     setpoint vertical velocity from trejectory
  z_setpoint      setpoint hight
  U               PI control output
  formula         linearized PI controler output
  Cd              resulting drag coefficiant 
  delta_S         resulting force
%}

%%%%%%%%%%%%%%%%%%%% TRAJECTORY SELECTION and REFERENCES COMPUTATION %%%%%%%%%%%%%%%%%%%%

%% Choose the nearest trajectory ( only at the first iteration )

[z_setpoint, Vz_setpoint, csett] = set_Trajectory(z, Vz, csett);

%%%%%%%%%%%%%%%%%%%% PID ALGORITHM %%%%%%%%%%%%%%%%%%%%

% If e>0 the rocket is too fast. I slow it down with Fx>0 --> open aerobrakes
% If e<0 the rocket is too slow. I speed it up with Fx<0 --> close aerobrakes

ro = getRho(z);

% Control variable limits
delta_S_max = 0.01;
Cd_min = getDrag(V_mod,z,0, csett.coeff_Cd);
Cd_max = getDrag(V_mod,z,delta_S_max, csett.coeff_Cd);
Umax = -csett.m*csett.g - 0.5*ro*Cd_min*csett.S0*Vz*V_mod   
Umin = -csett.m*csett.g - 0.5*ro*Cd_max*csett.S0*Vz*V_mod

% PID
error = (Vz_setpoint - Vz); % Changed the signum

P = csett.Kp_2*error;
if csett.saturation == false
    csett.I = csett.I + csett.Ki_2*error*csett.sample_time;
end

U = P + csett.I; 

% Anti-windup
[U, csett.saturation] = Saturation(U, Umin, Umax);


%%%%%%%%%%%%%%%%%%%% TRANSFORMATION FROM U to delta_S %%%%%%%%%%%%%%%%%%%%

% Get the Cd for each possible aerobrake surface
Cd_available = 1:length(csett.delta_S_available);
for ind = 1:length(csett.delta_S_available)
    Cd_available(ind) = getDrag(V_mod,z,csett.delta_S_available(ind), csett.coeff_Cd);
end
Cd_available = Cd_available';

% For all possible delta_S compute U
% Then choose the delta_S which gives an U which has the minimum error if compared with U_pid
U_linearization = -csett.m*csett.g-0.5*ro*csett.S0*Cd_available*Vz*V_mod;
[~, index_minimum] = min( abs(abs(U) - abs(U_linearization)) ); 
delta_S = csett.delta_S_available(index_minimum);  

% Just for plotting
pid = U;
formula = -csett.m*csett.g-0.5*ro*csett.S0*Cd_available(index_minimum)*Vz*V_mod;
Cd = Cd_available(index_minimum);



%%%%%%%%%%%%%%%%%%%% TRANSFORMATION FROM delta_S to SERVOMOTOR ANGLE DEGREES %%%%%%%%%%%%%%%%%%%%

alpha_rad      = (-csett.b + sqrt(csett.b^2 + 4*csett.a*delta_S)) / (2*csett.a);

% Alpha saturation ( possibili problemi per azione integrale ? )
[alpha_rad, ~] = Saturation(alpha_rad, 0, 0.89);

alpha_degree   = (alpha_rad*180)/pi;

%%%%%%%%%%%%%%%%%%%% LIMIT THE RATE OF THE CONTROL VARIABLE %%%%%%%%%%%%%%%%%%%%

alpha_degree = rate_Limiter(alpha_degree, csett.alpha_degree_prec, csett.rate_limiter, csett.sample_time);

alpha_degree = smooth_Control(alpha_degree, csett.alpha_degree_prec, csett.filter_coeff);
csett.alpha_degree_prec = alpha_degree;
end