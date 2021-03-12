function [alpha_degree, Vz_setpoint, z_setpoint, Cd, delta_S] = controlAlgorithmServoDegree(z,Vz,V_mod,sample_time)
%CONTROL_ALGORITHM_SERVO_DEGREE  Finds trejectory (z-Vz) to follow and uses a PI controler to follow the trejectory
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
%   Cd              resulting drag coefficiant 
%   delta_S         resulting force

% Define global variables
global data_trajectories coeff_Cd 
global Kp_3 Ki_3 I alpha_degree_prec index_min_value iteration_flag chosen_trajectory saturation


%%%%%%%%%%%%%%%%%%%% TRAJECTORY SELECTION and REFERENCES COMPUTATION %%%%%%%%%%%%%%%%%%%%

%% Choose the nearest trajectory ( only at the first iteration )
if iteration_flag == 1
    
    best_min   = inf;
    best_index = inf;

    for ind = 1:length(data_trajectories)
       
        % Select a z trajectory and a Vz trajectory (to speed up select only the first values, not ALL)
        z_ref  = data_trajectories(ind).Z_ref(1:50); 
        Vz_ref = data_trajectories(ind).V_ref(1:50); 
        distances_from_current_state = (z_ref-z).^2 + (Vz_ref-Vz).^2; 

        % Find the nearest point to the current trajectory
        [min_value, index_min_value] = min( distances_from_current_state ); 

        if (min_value < best_min)
             best_min = min_value;
             best_index = index_min_value;
             chosen_trajectory = ind;  
        end

    end

    index_min_value = best_index; % Save the actual index to speed up the research
    iteration_flag  = 0; % Don't enter anymore the if condition
    
    % I select the reference altitude and the reference vertical velocity
    z_setpoint  =  data_trajectories(chosen_trajectory).Z_ref(index_min_value);
    Vz_setpoint =  data_trajectories(chosen_trajectory).V_ref(index_min_value);

    
%% For the following iterations keep tracking the chosen trajectory
else

    % Select the z trajectory and the Vz trajectory 
    % To speed up the research, I reduce the vector at each iteration (add if-else for problem in index limits)
    z_ref  = data_trajectories(chosen_trajectory).Z_ref(index_min_value:end);  
    Vz_ref = data_trajectories(chosen_trajectory).V_ref(index_min_value:end);  

    % Find the value of the altitude in z_reference nearer to z_misured 
    [~, index_min_value] = min( abs(z_ref - z) );

    z_setpoint  =  z_ref(index_min_value);
    Vz_setpoint = Vz_ref(index_min_value);

end  



%%%%%%%%%%%%%%%%%%%% PID ALGORITHM %%%%%%%%%%%%%%%%%%%%

% Control variable limits
Umin = 0;  % degrees   
Umax = 48; % degrees

% PID
error = (Vz - Vz_setpoint); % > 0 (in teoria)

P = Kp_3*error;
if saturation == false
    I = I + Ki_3*error;
end

U = P + I;
    
if ( U < Umin)  
    U = Umin; % fully close
    saturation = true;                                         
elseif ( U > Umax) 
    U = Umax; % fully open
    saturation = true;                          
else
    saturation = false;
end


alpha_degree = U;



%%%%%%%%%%%%%%%%%%%% LIMIT THE RATE OF THE CONTROL VARIABLE %%%%%%%%%%%%%%%%%%%%

rate_limiter_max =  60/0.2; % datasheet: 60deg/0.13s --> increased for robustness
rate_limiter_min = -60/0.2;

rate = (alpha_degree - alpha_degree_prec) / sample_time;

if (rate > rate_limiter_max)
    alpha_degree = sample_time*rate_limiter_max + alpha_degree_prec;
elseif (rate < rate_limiter_min)
    alpha_degree = sample_time*rate_limiter_min + alpha_degree_prec;
end

%%%%%%%%%%%%%%%%%%%% Smooth the control variable with a filter %%%%%%%%%%%%%%%%%%%%
filter_coeff = 0.9;
alpha_degree = filter_coeff*alpha_degree + (1-filter_coeff)*alpha_degree_prec;

alpha_degree = round(alpha_degree);
alpha_degree_prec = alpha_degree;

end