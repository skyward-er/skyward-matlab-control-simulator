function [alpha_degree, Vz_setpoint, z_setpoint, pid, U_linear, Cd, delta_S, csett] = controlAlgorithm(z, Vz, V_mod, sample_time, csett)
%CONTROL_ALGORITHM  Finds trejectory (z-Vz) to follow and uses a PI
%controler to follow the trejectory and then transfere it with a to a force
%
%
%   INPUTS:
%   z               acutal hight of the rocket
%   Vz              actual vertical velocity of the rocket
%   V_mod           actual roket velocity in the direction of the main axis
%   sample_time     sample time of the control system
%
%   OUTPUTS:
%   alpha_degree    output angle for servo
%   Vz_setpoint     setpoint vertical velocity from trejectory
%   z_setpoint      setpoint hight
%   pid             PI control output
%   U_linear        linearized PI controler output
%   Cd              resulting drag coefficiant 
%   delta_S         resulting force

%% Choose the nearest trajectory ( only at the first iteration )
if csett.iteration_flag == 1
    
    best_min   = inf;
    best_index = inf;

    for ind = 1:length(csett.data_trajectories)
       
        % Select a z trajectory and a Vz trajectory (to speed up select only the first values, not ALL)
        z_ref  = csett.data_trajectories(ind).Z_ref(1:50); 
        Vz_ref = csett.data_trajectories(ind).V_ref(1:50); 
        err    = (z_ref-z).^2 + (Vz_ref-Vz).^2; 

        % Find the nearest point to the current trajectory
        [min_value, index_min_value] = min( err ); 

        if (min_value < best_min)
             best_min = min_value;
             best_index = index_min_value;
             csett.chosen_trajectory = ind;  
        end

    end

    csett.index_min_value = best_index; % Save the actual index to speed up the research
    csett.iteration_flag  = 0; % Don't enter anymore the if condition
    
    % I select the reference altitude and the reference vertical velocity
    z_setpoint  =  csett.data_trajectories(csett.chosen_trajectory).Z_ref(csett.index_min_value);
    Vz_setpoint =  csett.data_trajectories(csett.chosen_trajectory).V_ref(csett.index_min_value);

    
%% For the following iterations keep tracking the chosen trajectory
else

    % Select the z trajectory and the Vz trajectory 
    % To speed up the research, I reduce the vector at each iteration (add if-else for problem in index limits)
    z_ref  = csett.data_trajectories(csett.chosen_trajectory).Z_ref(csett.index_min_value:end);  
    Vz_ref = csett.data_trajectories(csett.chosen_trajectory).V_ref(csett.index_min_value:end);  

    % 1) Find the value of the altitude in z_reference nearer to z_misured 
    [~, csett.index_min_value] = min( abs(z_ref - z) );

    % 2) Find the reference using Vz(z)
    z_setpoint  =  z_ref(csett.index_min_value);
    Vz_setpoint = Vz_ref(csett.index_min_value);
end  



%%%%%%%%%%%%%%%%%%%% PI ALGORITHM %%%%%%%%%%%%%%%%%%%%

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
if ( U < Umin)  
    U = Umin; % fully close
    csett.saturation = true;                                         
elseif ( U > Umax) 
    U = Umax; % fully open
    csett.saturation = true;                          
else
    csett.saturation = false;
end



%%%%%%%%%%%%%%%%%%%% TRANSFORMATION FROM U to delta_S %%%%%%%%%%%%%%%%%%%%


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



%%%%%%%%%%%%%%%%%%%% TRANSFORMATION FROM delta_S to SERVOMOTOR ANGLE DEGREES %%%%%%%%%%%%%%%%%%%%

% delta_S [m^2] = (-9.43386 * alpha^2 + 19.86779 * alpha) * 10^(-3). Alpha belongs to [0 ; 0.89 rad]

alpha_rad = (-csett.b + sqrt(csett.b^2 + 4*csett.a*delta_S)) / (2*csett.a);

% Alpha saturation ( possibili problemi per azione integrale ? )
if (alpha_rad < 0)
    alpha_rad = 0;
elseif (alpha_rad > 0.89)
    alpha_rad = 0.89;
end

alpha_degree = (alpha_rad*180)/pi;



%%%%%%%%%%%%%%%%%%%% LIMIT THE RATE OF THE CONTROL VARIABLE %%%%%%%%%%%%%%%%%%%%

rate = (alpha_degree - csett.alpha_degree_prec) / sample_time;

if (rate > csett.rate_limiter_max)
    alpha_degree = sample_time*csett.rate_limiter_max + csett.alpha_degree_prec;
elseif (rate < csett.rate_limiter_min)
    alpha_degree = sample_time*csett.rate_limiter_min + csett.alpha_degree_prec;
end

% Smooth the control variable with a filter

alpha_degree = csett.filter_coeff*alpha_degree + (1-csett.filter_coeff)*csett.alpha_degree_prec;

alpha_degree = round(alpha_degree);
csett.alpha_degree_prec = alpha_degree;

end