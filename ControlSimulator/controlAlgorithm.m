function [alpha_degree] = controlAlgorithm(z,Vz,V_mod,sample_time)

% Define global variables
global coeff_Cd data_1 data_2 data_3 data_4 data_5 data_6 data_7 data_8 data_9 data_10 data_11
global Kp Ki I alpha_degree_prec iteration_flag chosen_trajectory saturation


%% TRAJECTORY SELECTION and REFERENCES COMPUTATION

data_trajectories = [data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11];

if iteration_flag == 1 % Choose the nearest trajectory ( only at the first iteration )
    
best_min = inf;
best_index = inf;

for ind = 1:length(data_trajectories)

% Select a z trajectory and a Vz trajectory
z_ref = data_trajectories(ind).Z_ref;
Vz_ref = data_trajectories(ind).V_ref;

% Find the value of z_reference nearer to z_misured
[min_value, index_min_value] = min( abs(z_ref - z) ); 

if (min_value < best_min)
    best_min = min_value;
    best_index = index_min_value;
    chosen_trajectory = ind;
end

end

% I select the reference altitude and the reference vertical velocity
z_setpoint = z_ref(best_index);
Vz_setpoint = Vz_ref(best_index);

iteration_flag = 0; % Don't enter anymore the if condition

else  % For the following iterations keep tracking the chosen trajectory
    
% Select the z trajectory and the Vz trajectory
z_ref =  data_trajectories(chosen_trajectory).Z_ref;
Vz_ref = data_trajectories(chosen_trajectory).V_ref;

% Find the value of z_reference nearer to z_misured
[~, index_min_value] = min( abs(z_ref - z) ); 

% I select the reference altitude and the reference vertical velocity
z_setpoint = z_ref(index_min_value);
Vz_setpoint = Vz_ref(index_min_value);

end


%% PID ALGORITHM

Umin = 0;      % F_drag_min = 0
Umax = 1000;   % F_drag_max = 0.5*1.225*(0.0201+0.01)*1*250^2
dt = 0.1; % ASK THE FINAL STEP TIME !!!!!!!!!!

error = (Vz - Vz_setpoint);

P = Kp*error;

if saturation == false
I = I + Ki*error*dt;
end

U = P + I;
    
if ( U < Umin)  
U=Umin; 
saturation = true;                                         
elseif ( U > Umax) 
U=Umax; 
saturation = true;                          
else
saturation = false;
end

%% TRANSFORMATION FROM U to delta_S

% If I forecast an overshoot, e>0, u>0, Fx>0 --> closed aerobrakes
% If I forecast an undershoot, e<0, u<0, Fx<0 --> open aerobrakes

% Parameters
ro = getRho(z);
m = 22;
g = 9.81;
diameter = 0.15; 
S0 = (pi*diameter^2)/4;  

% Range of values for the control variable
delta_S_available = 0.0:0.001:0.01; 

% Get the Cd for each possible aerobrake surface
Cd_available = 1:length(delta_S_available);
for ind = 1:length(delta_S_available)
Cd_available(ind) = getDrag(V_mod,z,delta_S_available(ind), coeff_Cd);
end

delta_S_available = delta_S_available';
Cd_available = Cd_available';

% For all possible delta_S compute Fdrag.
% Then choose the delta_S which gives an Fdrag which has the minimum error if compared with F_drag_pid
[~, index_min_value] = min( abs(U - 0.5*ro*S0*Cd_available*Vz*V_mod) ); 

% Cd_available = Cd_available(index_min_value)
delta_S = delta_S_available(index_min_value);  % delta_S belongs to [0; 0.01]

%% TRANSFORMATION FROM delta_S to SERVOMOTOR ANGLE DEGREES

% delta_S [m^2] = (-9.43386 * alpha^2 + 19.86779 * alpha) * 10^(-3). Alpha belongs to [0 ; 0.89 rad]
a = -9.43386/1000;
b = 19.86779/1000;

alpha_rad = (-b + sqrt(b^2 + 4*a*delta_S)) / (2*a);
% alpha_rad_rad = (-b - sqrt(b^2 + 4*a*delta_S)) / (2*a);

% Alpha saturation
if (alpha_rad < 0)
    alpha_rad = 0;
elseif (alpha_rad > 0.89)
    alpha_rad = 0.89;
end

alpha_degree = (alpha_rad*180)/pi

%% LIMIT THE RATE OF THE CONTROL VARIABLE

rate_limiter_max = 60/0.13; % 60deg/0.13s
rate_limiter_min = -60/0.13;

rate = (alpha_degree - alpha_degree_prec) / sample_time;

if (rate > rate_limiter_max)
    alpha_degree = sample_time*rate_limiter_max + alpha_degree_prec;
elseif (rate < rate_limiter_min)
    alpha_degree = sample_time*rate_limiter_min + alpha_degree_prec;
end

alpha_degree_prec = alpha_degree;

% Testing:
% alpha_degree = 25;

end



