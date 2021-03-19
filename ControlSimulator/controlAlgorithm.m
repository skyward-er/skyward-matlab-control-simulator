function [alpha_degree, delta_S, Vz_setpoint, z_setpoint, Vx_setpoint, Cd] = controlAlgorithm(z,Vz,x,Vx,V_mod,sample_time)

% Define global variables
global data_trajectories coeff_Cd 
global delta_S_prec alpha_degree_prec index_min_value iteration_flag chosen_trajectory 

%% TRAJECTORY SELECTION and REFERENCES COMPUTATION

if iteration_flag == 1 % Choose the nearest trajectory ( only at the first iteration )
    
    best_min   = inf;
    best_index = inf;

    for ind = 1:length(data_trajectories)
       
        % Select a z trajectory and a Vz trajectory (to speed up select only the first values, not ALL)
        z_ref  = data_trajectories(ind).Z_ref(1:50); 
        Vz_ref = data_trajectories(ind).VZ_ref(1:50); 
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
    Vz_setpoint =  data_trajectories(chosen_trajectory).VZ_ref(index_min_value);
    Vx_setpoint =  data_trajectories(chosen_trajectory).VX_ref(index_min_value);

else  % For the following iterations keep tracking the chosen trajectory

    % Select the z trajectory and the Vz trajectory 
    % To speed up the research, I reduce the vector at each iteration (add if-else for problem in index limits)
    z_ref  = data_trajectories(chosen_trajectory).Z_ref(index_min_value:end);  
    Vz_ref = data_trajectories(chosen_trajectory).VZ_ref(index_min_value:end);   
    Vx_ref = data_trajectories(chosen_trajectory).VX_ref(index_min_value:end);  

    % 1) Find the value of the altitude in z_reference nearer to z_misured 
    [~, index_min_value] = min( abs(z_ref - z) );

    z_setpoint  = z_ref(index_min_value);
    Vz_setpoint = Vz_ref(index_min_value);
    Vx_setpoint = Vx_ref(index_min_value);

end  


%% PARAMETERS

T        = sample_time;
m        = 22;
g        = 9.81;
rho      = getRho(z);
diameter = 0.15; 
S0       = (pi*diameter^2)/4;   

S        = S0 + delta_S_prec;
Cd_fake  = getDrag(V_mod,z,delta_S_prec, coeff_Cd); % Però questo non è Cd, perchè tiene conto anche di delta_S
Cd       = (S0*Cd_fake)/S;
    
%% LQR ALGORITHM

% States: z,Vz,Vx

%% Original matrix
Q = [0.7,    0,           0,   0;
       0,    1,           0,   0;
       0,    0,       0.001,   0;
       0,    0,           0,   1];
   
R = 65000; 

%% Test
% S_max_squared = 0.01^2;
% 
% Q = [0.7,    0,            0,           0;
%        0,    5,            0,           0;
%        0,    0,       0.00005,            0;
%        0,    0,            0,      0.0001];
%    
% R = 0.6/S_max_squared; 

% Linearize the system around the current state
% % % A = [1,                                                                                        T, 0,                                                                                        0;
% % %      0, 1 - T*((Cd*S*rho*(Vx^2 + Vz^2)^(1/2))/(2*m) + (Cd*S*Vz^2*rho)/(2*m*(Vx^2 + Vz^2)^(1/2))), 0,                                            -(Cd*S*T*Vx*Vz*rho)/(2*m*(Vx^2 + Vz^2)^(1/2));
% % %      0,                                                                                        0, 1,                                                                                        T;
% % %      0,                                            -(Cd*S*T*Vx*Vz*rho)/(2*m*(Vx^2 + Vz^2)^(1/2)), 0, 1 - (Cd*S*T*Vx^2*rho)/(2*m*(Vx^2 + Vz^2)^(1/2)) - (Cd*S*T*rho*(Vx^2 + Vz^2)^(1/2))/(2*m)];
% % %  
% % % B = [                                        0;
% % %       -(Cd*T*Vz*rho*(Vx^2 + Vz^2)^(1/2))/(2*m);
% % %                                              0;
% % %        -(Cd*T*Vx*rho*(Vx^2 + Vz^2)^(1/2))/(2*m)];
   
A = [1,                                                                                        T,                                                                                        0, 0;
     0, 1 - T*((Cd*S*rho*(Vx^2 + Vz^2)^(1/2))/(2*m) + (Cd*S*Vz^2*rho)/(2*m*(Vx^2 + Vz^2)^(1/2))),                                            -(Cd*S*T*Vx*Vz*rho)/(2*m*(Vx^2 + Vz^2)^(1/2)), 0;
     0,                                            -(Cd*S*T*Vx*Vz*rho)/(2*m*(Vx^2 + Vz^2)^(1/2)), 1 - T*((Cd*S*rho*(Vx^2 + Vz^2)^(1/2))/(2*m) + (Cd*S*Vx^2*rho)/(2*m*(Vx^2 + Vz^2)^(1/2))), 0;
    -1,                                                                                        0,                                                                                        0, 0];
 
B = [                                       0;
     -(Cd*T*Vz*rho*(Vx^2 + Vz^2)^(1/2))/(2*m);
     -(Cd*T*Vx*rho*(Vx^2 + Vz^2)^(1/2))/(2*m);
                                            0];
                                             
x_measured  = [z, Vz, Vx, 0]';
x_reference = [z_setpoint, Vz_setpoint, Vx_setpoint, 0]';
x_error     =  x_measured - x_reference

% Solve Riccati equation
P       = Q;   % Initial guess for P    
maxiter = 100;
eps     = 0.01;

for i=1:maxiter
    Pn = A' * P * A - A' * P * B * inv(R + B' * P * B) * B' * P * A + Q;
    if (max(max((abs(Pn - P))))) < eps % Continue to compute P until the actual and previous solution are almost equal
        break
    end
    P = Pn;
end

K = inv(B' * P * B + R) * B' * P * A;
U = -K*x_error

% Debug
J_z  = Q(1,1)*x_error(1)^2
J_Vz = Q(2,2)*x_error(2)^2
J_Vx = Q(4,4)*x_error(4)^2

J_Q = x_error'*Q*x_error
J_R = U'*R*U

% Control variable limits
Umin = 0;     
Umax = 0.01;

if ( U < Umin)  
    U = Umin; % fully close                                      
elseif ( U > Umax) 
    U = Umax; % fully open                       
end

filter_coeff = 0.9;
delta_S = filter_coeff*U + (1-filter_coeff)*delta_S_prec;

%% TRANSFORMATION FROM delta_S to SERVOMOTOR ANGLE DEGREES

% delta_S [m^2] = (-9.43386 * alpha^2 + 19.86779 * alpha) * 10^(-3). Alpha belongs to [0 ; 0.89 rad]
a = -9.43386/1000;
b = 19.86779/1000;

alpha_rad    = (-b + sqrt(b^2 + 4*a*delta_S)) / (2*a);
alpha_degree = (alpha_rad*180)/pi;

%% LIMIT THE RATE OF THE CONTROL VARIABLE

rate_limiter_max =  60/0.2; % datasheet: 60deg/0.13s --> increased for robustness
rate_limiter_min = -60/0.2;

rate = (alpha_degree - alpha_degree_prec) / sample_time;

if (rate > rate_limiter_max)
    alpha_degree = sample_time*rate_limiter_max + alpha_degree_prec;
elseif (rate < rate_limiter_min)
    alpha_degree = sample_time*rate_limiter_min + alpha_degree_prec;
end

alpha_degree      = round(alpha_degree);
alpha_degree_prec = alpha_degree;

alpha_rad    = (alpha_degree*pi)/180;
delta_S_prec = a * alpha_rad^2 + b * alpha_rad;

end