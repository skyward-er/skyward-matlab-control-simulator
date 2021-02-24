%% PARAMETER DEFINITION OF RTHE Full_state_Kalman_filter.m file

% %------Rocket parameters-------
% 
% M           =   [0;1e-3;0];         %[N·m]Free torque dynamics
% 
% Ix          =   0.05;               %[kg·m^2]
% 
% Iy          =   0.05;               %[kg·m^2]
% 
% Iz          =   0.01;               %[kg·m^2]
% 
% I           =   diag([Ix,Iy,Iz]);   %consideration of the body axis to be aligned 
%                                 %with the principal intertia axis
% m           =   22;                 %Dry mass of the rocket
% 
% rho         =   1.25;               %[kg/m^3] Density of air at ground 
%                                 %(overestimation of drag)
% A_d         =   pi*(0.75/2)^2;      %Cross area of the rocket assuming a diameter of 0.75 m
% 
% CD          =   0.25;               %Approximated drag coefficient for circular section                          
% 
% %--------Simulation parameters-------
% 
% t_sim       =   40;                 %time of simulation

%--------Enviromental parameters---------

m_inertial  =   [1 0 0];           %Initial vector direction of the magnetic field --> the inertial frame is the 

%--------Initial conditions--------

% w0          =   [0;0.02;0];         %[rad/s] Initial condition of the angular velocity
% 
% h0          =   I*w0;               %Initial value of angular momentum
% 
% q0          =   [1;2;3;1];          %Initial value of the attitude (quaternion)
% 
% q0          =   q0/norm(q0);        %Normalisation of the initial quaternion
% 
% r0          =   [0;0;-500];          %[m] Initial position vector
% 
% v0          =   [2;0;-300];          %[m/s] Initial velocity

%---------Measurement parameters----------------

% f_sample1   =   100;                %[Hz] Frequency of the measurement of the IMU

% sigma_w     =   0.05;               % Variance of the angular velocity
% 
% sigma_a     =   0.3;               % Variance of the linear acceleration
 
% f_sample_mag=20;             %[Hz] Frequency of measurement of magnetometer 
% 
% f_sample_h=20;               %[Hz] Frequency of measurement of barometer
% 
% f_sample_pitot=20;           %[Hz] Frequency of measurement of pitot tube 
% 
% f_sample_GPS=5;             %[Hz] Frequency of measurement of GPS

                                    
sigma_mag   =   0.05;               % This is the variance of the direction 
                                    % of the magnetic field, not of the magnitude 
                                    % itself, so the variance used in the 
                                    % initSensors file is not valid. Since 
                                    % the transition from modulus to direction 
                                    % is non-linear, and the sigmas are only 
                                    % weighting parameters, i think this can 
                                    % be used as a parameter for the kalman                                    
                                    
sigma_baro     =   4;               % Variance of the altitude measurement 
                                    % with barometer
                                    % Since the avilable sigma is on
                                    % pressure measurements and the
                                    % tranformation is also non-linear,
                                    % this is also set as parameter

sigma_GPS   =   2;                  % Variance of position due to GPS. This 
                                    % should be the same number as the sqrt of 
                                    % GPS_NEOM9N.noiseVariance from
                                    % initSensors.m
                                    
f_sample1   =   100;                %[Hz] Frequency of the measurement of the IMU

Dt_gyro     = 1/f_sample1;

sigma_w     =   100*(1000*pi/180)^2;   % Variance of the angular velocity;
                                       % This must be the number in line 39
                                       % of initSensors GYRO_LSM9DS1.noiseVariance
                                       % multiplied by (1000*pi/180)^2
                                    
sigma_mag   =   1e-3;               %Variance of the magnetometer                                    
                                    
sigma_beta  =   1e-4;               %Variance of the bias
      
sigma_g   = sqrt(sigma_w^2/Dt_gyro+sigma_beta^2*Dt_gyro/12);




%---------Estimation parameters-------------------

dt_k        =   0.01;               %Time step for the estimation

QLinear     =   100*...
                [1     0     0      0      0      0;   %Noise covariance matrix
                 0     1     0      0      0      0;    %for the linear dynamics
                 0     0     1      0      0      0;
                 0     0     0      0.1    0      0;
                 0     0     0      0      0.1    0;
                 0     0     0      0      0      0.1];
             

% EXPRESSION OF Q OBTAINED FROM ''OPTIMAL ESTIMATION OF DYNAMIC SYSTEMS'', 
% CRASSIDIS AND JUNKINS (pg. 459) 

Qq = [(sigma_w^2*dt_k+(1/3)*sigma_beta^2*dt_k^3)*eye(3)     0.5*sigma_beta^2*dt_k^2*eye(3);
        0.5*sigma_beta^2*dt_k^2*eye(3)                sigma_beta^2*dt_k*eye(3)];

                  
                  
                  
                  
                  
                  
