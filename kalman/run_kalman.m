function [x_c, P_c, kalman] = run_kalman(x_prev, P_prev, sp, kalman, mag_NED)

% Author: Alejandro Montero
% Co-Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@kywarder.eu
% email: alejandro.montero@skywarder.eu, alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

%{
-----------DESCRIPTION OF FUNCTION:------------------
This function implement a "Multiplicative extended kalman filter" that takes 
the angular velocity, the acceleration and the sensors outputs from the
previous integration and estimates the current vector of state of the rocket (x_c)
and its covariance matrix (P_c).
For more information check the navigation system report 
    -INPUTS: 
        - x_prev: 1x13 VECTOR OF PREVIOUS VALUES -->
                  STATES: [ N, E, D, Vn, Ve, Vd, q1, q2, q3, q4, b1, b2, b3]

        - P_prev: 13x13 MATRIX OF PREVIOUS COVARIANCE OF STATE

        - sp      STRUCT THAT CONTAIN ALL THE MEASUREMENTS THAT THE KALMAN
                  NEED.

            - t_v:   TIME VECTOR OF THE PREVIOUSLY INTEGRATED INSTANTS. SPANS 
                     0.1 SECONDS AND SHOULD HAVE A TIME STEP OF 0.01 SINCE 
                     KALMAN FILTER RUNS AT 100 HZ (10 TIME INTANTS IN TOTAL)
    
            - a_v:    ACCELERATION VECTOR MEASURED DURING THE INTEGRATED PERIOD.
                      SINCE IT ALSO RUNS AT 100 HZ AND THE INTEGRATION PERIOD 
                      SPANS 0.1 SECONDS: [10x3] [m/s^2]

            - w_v:    ANGULAR VELOCITY VECTOR MEASURED DURING THE INTEGRATED 
                      PERIOD. SINCE IT ALSO RUNS AT 100 HZ AND THE INTEGRATION 
                      PERIOD SPANS 0.1 SECONDS: [10x3] [rad/s]


            - t_baro: VECTOR OF TIME INSTANTS AT WHICH THE BAROMETER TOOK 
                      SAMPLES INSIDE THE 0.1 INTEGRATION PERIOD. SINCE IT RUNS
                      AT 20 HZ:[2x1] [s] 

            - baro:   CORRESPONDING ALTITUDE MEASUREMENTS FROM THE
                      BAROMETER [2x1] [m]


            - t_mag:  VECTOR OF TIME INSTANTS AT WHICH THE MGNETOMETER TOOK 
                      SAMPLES INSIDE THE 0.1 INTEGRATION PERIOD. SINCE IT RUNS
                      AT 20 HZ: [2x1] [s] 

            - mag:    CORRESPONDING MAGNETIC FIELD MEASUREMENTS FROM THE 
                      MAGNETOMETER: [2x3] [-] --> NORMALIZED INSIDE THE MAG
                      CORRECTION

            - t_GPS:  VECTOR OF TIME INSTANTS AT WHICH THE GPS TOOK 
                      SAMPLES INSIDE THE 0.1 INTEGRATION PERIOD. SINCE IT RUNS
                      AT 10 HZ: [1x1] [s]

            - GPS:    CORRESPONDING POSITION MEASUREMENTS FROM THE GPS:
                      [1x3] [m]

            - vGPS:   CORRESPONDING VELOCITY MEASUREMENTS FROM THE GPS.
                      [1x3] [m/s]

        - kalman:     STRUCT THAT CONTAIN ALL THE KALMAN PARAMETERS:

            - sigma_GPS:  STANDARD DEVIATION OF THE GPS (SQRT OF
                          VARIANCE) [1x1]

            - sigma_baro: STANDARD DEVIATION OF THE BAROMETER (SQRT OF
                          VARIANCE) [1x1]

            - sigma_mag:  STANDARD DEVIATION OF THE MAGNETOMETER (SQRT OF
                          VARIANCE) [1x1]

            - Qlin:       COVARIANCE MATRIX OF (LINEAR) PROCESS NOISE [6x6]

            - Qq:         COVARIANCE MATRIX OF (QUATERNION) PROCESS NOISE [4x4]  

        - mag_NED: INERTIAL DIRECTION OF THE MAGNETIC FIELD: [3x1]




        - nsat :   NUMBER OF SATELLITES VAILABLE FOR THE GPS 1x1 [-]
      
        - fix:     FIX OF THE GPS; BINARY VARIABLE TO DETERMINE WHETHER TO
                   TRUST THE GPS OR NOT [-]



      -OUTPUTS:
        - x_c:   CORRECTED VECTOR OF STATES OF THE ROCKET. CONTAINS ALL THE
                 ESTIMATIONS FOR EACH TIME INSTANT IN t_v --> [10x13]

        - P_c:   CORRECTED COVARIANCE MATRIX FOR EACH OF THE ESTIMATION TIME
                 INSTANTS. [12x12x10]
-----------------------------------------------------------------------
%}
tv          =   sp.t_acc;
[fix, nsat] =   gpsFix(sp.accel);

dt_k        =   tv(2)-tv(1);                 % Time step of the kalman
x_lin       =   zeros(length(tv),6);         % Pre-allocation of corrected estimation
xq          =   zeros(length(tv),7);         % Pre-allocation of quaternions and biases
x_c         =   zeros(length(tv),13);

P_c         =   zeros(12,12,length(tv));
P_lin       =   zeros(6,6,length(tv));       %Pre-allocation of the covariance matrix
P_q         =   zeros(6,6,length(tv));

x_lin(1,:)  =   x_prev(1:6);                 % Allocation of the initial value
xq(1,:)     =   x_prev(7:13);
x_c(1,:)    =   [x_lin(1,:),xq(1,:)];

P_lin(:,:,1)=   P_prev(1:6,1:6);
P_q(:,:,1)  =   P_prev(7:12,7:12);
P_c(:,:,1)  =   P_prev;

index_GPS=1;
index_bar=1;
index_mag=1;

for i=2:length(tv)
    %Prediction part
    [x_lin(i,:),P_lin(:,:,i)] = predictorLinear(x_lin(i-1,:),dt_k,...
                                                P_lin(:,:,i-1),sp.accel(i-1,:),xq(i-1,1:4),kalman.QLinear);
    
    [xq(i,:),P_q(:,:,i)]       = predictorQuat(xq(i-1,:),P_q(:,:,i-1),...
                                               sp.gyro(i-1,:),dt_k,kalman.Qq);                   
    %Corrections
     if tv(i) >= sp.t_gps(index_GPS)              %Comparison to see the there's a new measurement
       [x_lin(i,:),P_lin(:,:,i),~]     = correctionGPS2(x_lin(i,:),P_lin(:,:,i),sp.gps(index_GPS,1:2),...
                                                        sp.gpsv(index_GPS,1:2),kalman.sigma_GPS,nsat,fix);
        index_GPS   =  index_GPS + 1;
     end
    
    if tv(i) >= sp.t_baro(index_bar)              %Comparison to see the there's a new measurement
       [x_lin(i,:),P_lin(:,:,i),~]     = correctionBarometer(x_lin(i,:),P_lin(:,:,i),sp.h_baro(index_bar),kalman.sigma_baro);
        index_bar   =  index_bar + 1;     
   end
     
    if tv(i) >= sp.t_mag(index_mag)               %Comparison to see the there's a new measurement
       [xq(i,:),P_q(:,:,i),~,~]    = correctorQuat(xq(i,:),P_q(:,:,i),sp.mag(index_mag,:),kalman.sigma_mag,mag_NED);
       index_mag    =  index_mag + 1;  
    end
    x_c(i,:) = [x_lin(i,:),xq(i,:)];
    P_c(1:6,1:6,i)   = P_lin(:,:,i);
    P_c(7:12,7:12,i) = P_q(:,:,i);
    
    if kalman.flag_apo  == false
            if -x_c(i,6) < kalman.v_thr && -x_c(i,3) > 100
                kalman.counter = kalman.counter + 1;
            else
                kalman.counter = 0;
            end
            if kalman.counter >= kalman.count_thr
               kalman.t_kalman = tv(i);
               kalman.flag_apo = true;
            end
    end
end
end
