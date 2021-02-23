function [x_c,P_c]=run_kalman(x_prev,P_prev,t_v,a_v,w_v,t_baro,baro,sigma_baro,...
                              t_mag,mag,sigma_mag,...
                              t_GPS,GPS,sigma_GPS,n_sats,fix,Q)
%23/02/2021 ANY QUESTIONS CAN BE DIRECTED TO ALEJANDRO MONTERO FROM SKYWARD

%-----------DESCRIPTION OF FUNCTION:------------------
% This function takes the angular velocity, the acceleration and the
% sensors outputs from the previous integration and estimates the current
% vector of state of the rocket (x_c) and its covariance matrix (P_c)
%     -INPUTS: 
%         -x_prev:  1x10 VECTOR OF PREVIOUS VALUES --> 3 FIRST STATES
%                   ARE X , Y AND H, THE FOLLOWING THREE ARE VX, VY AND VZ
%                   AND THE LAST 4 ARE THE QUATERNION COMPONENTS
%
%         -P_prev:  10x10 MATRIX OF PREVIOUS COVARIANCE OF STATE
%
%         -t_v:     TIME VECTOR OF THE PREVIOUSLY INTEGRATED INSTANTS. SPANS 
%                   0.1 SECONDS AND SHOULD HAVE A TIME STEP OF 0.01 SINCE 
%                   KALMAN FILTER RUNS AT 100 HZ (10 TIME INTANTS IN TOTAL)
%                   s
%
%         -a_v:     ACCELERATION VECTOR MEASURED DURING THE INTEGRATED PERIOD.
%                   SINCE IT ALSO RUNS AT 100 HZ AND THE INTEGRATION PERIOD 
%                   SPANS 0.1 SECONDS, THIS VECTOR HAS SIZE 10x3. m/s^2
%
%         -w_v:     ANGULAR VELOCITY VECTOR MEASURED DURING THE INTEGRATED 
%                   PERIOD. SINCE IT ALSO RUNS AT 100 HZ AND THE INTEGRATION 
%                   PERIOD SPANS 0.1 SECONDS, THIS VECTOR HAS SIZE 10x3.
%                   rad/s
% 
%         -t_baro:  VECTOR OF TIME INSTANTS AT WHICH THE BAROMETER TOOK 
%                   SAMPLES INSIDE THE 0.1 INTEGRATION PERIOD. SINCE IT RUNS
%                   AT 20 HZ, THIS IS A 2x1 VECTOR. s
%
%         -baro:    CORRESPONDING ALTITUDE MEASUREMENTS FROM THE BAROMETER.
%                   2x1. m
%
%         -sigma_baro:  STANDARD DEVIATION OF THE BAROMETER (SQRT OF
%                       VARIANCE) 1x1
%
%         -t_mag:   VECTOR OF TIME INSTANTS AT WHICH THE MGNETOMETER TOOK 
%                   SAMPLES INSIDE THE 0.1 INTEGRATION PERIOD. SINCE IT RUNS
%                   AT 20 HZ, THIS IS A 2x1 VECTOR. s
%
%         -mag:     CORRESPONDING MAGNETIC FIELD MEASUREMENTS FROM THE 
%                   MAGNETOMETER. 2x3 [-] --> SHOULD BE A NON-DIMENSIONAL
%                   VECTOR
%
%         -sigma_mag:   STANDARD DEVIATION OF THE MAGNETOMETER (SQRT OF
%                       VARIANCE) 1x1
%
%         -t_GPS:   VECTOR OF TIME INSTANTS AT WHICH THE GPS TOOK 
%                   SAMPLES INSIDE THE 0.1 INTEGRATION PERIOD. SINCE IT RUNS
%                   AT 10 HZ, THIS IS A 1x1 VECTOR. s
%
%         -GPS:     CORRESPONDING POSITION MEASUREMENTS FROM THE PITOT.
%                   1x3 m
%
%         -sigma_GPS:   STANDARD DEVIATION OF THE GPS (SQRT OF
%                       VARIANCE) 1x1
%
%         -n_sats: NUMBER OF SATELLITES VAILABLE FOR THE GPS 1x1 [-]
%       
%         -fix:     FIX OF THE GPS; BINARY VARIABLE TO DETERMINE WHETHER TO
%                   TRUST THE GPS OR NOT [-]
%
%         -Q:       COVARIANCE MATRIX OF PROCESS NOISE. 10x10   
%
%       -OUTPUTS:
%         -x_c: CORRECTED VECTOR OF STATES OF THE ROCKET. CONTAINS ALL THE
%               ESTIMATIONS FOR EACH TIME INSTANT IN t_v --> 10x10
%
%         -P_c: CORRECTED COVARIANCE MATRIX FOR EACH OF THE ESTIMATION TIME
%               INSTANTS. 10x10x10
% -----------------------------------------------------------------------
dt_k        = t_v(1)-t_v(2); %Time step of the kalman
x_c         = zeros(length(t_v),10); %Pre-allocation of corrected estimation
x_c(1,:)    = x_prev;                 %Allocation of the initial value
P_c         = zeros(10,10,length(t_v)); %Pre-allocation of the covariance matrix
P_c(:,:,1)  = P_prev;
index_GPS=1;
index_bar=1;
index_mag=1;
for i=2:length(t_v)
    %Prediction part
    [x_c(i,:),P_c(:,:,i)] = kalmanFilterPrediction(x_c(i-1,:),dt_k,...
                            P_c(:,:,i-1),a_v(i-1,:),w_v(i-1,:),Q);
    
    %Corrections
     if t_v(i)>=t_GPS(index_GPS)  %Comparison to see the there's a new measurement
       [x_c(i,:),P_c(:,:,i),~]     = correctionGPS(x_c(i,:),P_c(:,:,i),GPS(index_GPS,1),...
                            GPS(index_GPS,2),GPS(index_GPS,3),sigma_GPS,n_sats,fix);
        index_GPS   =  index_GPS + 1;
     end
    
    if t_v(i)>=t_baro(index_bar) %Comparison to see the there's a new measurement
       [x_c(i,:),P_c(:,:,i),~]     = correctionBarometer(x_c(i,:),P_c(:,:,i),baro(index_bar),sigma_baro);
        index_bar   =  index_bar + 1;     
    end
%      
%     if t_v(i)>=t_mag(index_mag) %Comparison to see the there's a new measurement
%        [x_c(i,:),P_c(:,:,i),~]     = correctionMagnetometer(x_c(i,:),P_c(:,:,i),mag(index_mag,:),sigma_mag);
%        index_mag    =  index_mag + 1;  
%     end
%     norm(x_c(i,7:10))
end
end







