function [xp, xv, Pout, ada]   =  run_ADA(xin, Pin, sensorData, ada)

% Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@skywarder.eu
% email: alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

%{ 
-----------DESCRIPTION OF FUNCTION:------------------
This function simulate the Apogee Detection Algorithm (ADA), the
algorithm is a 3rd order linear kalman filter with this states:
- x(1): Vertical position
- x(2): Vertical velocity
- x(3): Vertical acceleration

The filter model is a constant acceleration point mass, with trial
computation we will end up with this discrete time system:
- p(t+1) = pt + v*dt +0.5*a*dt^2
- v(t+1) = vt + a*dt
- a(t+1) = const
with observation on the first state y(t) = p(t)

INPUTS: 
    - xin:        [1x3] VECTOR OF PREVIOUS VALUES --> p,v,a

    - Pin:        [3x3] MATRIX OF PREVIOUS COVARIANCE OF STATE

    - t_baro:     [1x2] VECTOR OF TIME INSTANTS AT WHICH THE BAROMETER TOOK 
                        SAMPLES INSIDE THE 0.1 INTEGRATION PERIOD 

    - p_baro:     [1x2] PRESSURE MEASUREMENTS FROM THE BAROMETER

    - ada:              STRUCT WITH THE ADA PARAMETERS:

        - Q:          [3x3] COVARIANCE MATRIX OF PROCESS NOISE

        - R:          [1x1] COVARIANCE MATRIX OF OBSERVATION NOISE

        - flag_apo:   BOOLEAN VARIABLE THAT BECOME TRUE AFTER THE APOGEE IS
                      DETECTED count_threshold TIMES

        - t_ada:      TIME OF THE PREDICTED APOGEE 

        - counter:    COUNTER OF THE NUMBER OF CONSECUTIVES APOGEE
                      DETECTED

        - v_thr:      VELOCITY THRESHOLD FOR THE APOGEE DETECTION

        - count_thr:  COUNTER THRESHOLD FOR THE APOGEE DETECTION

      OUTPUTS:

        - xp:        [1x3] FILTERED STATES OF THE ADA 

        - xv:        [1x2] VERTICAL POSITION AND VELOCITY

        - Pout:      [3x3] PROPAGATED COVARIANCE MATRIX FOR EACH OF THE ESTIMATION TIME
                           INSTANTS. 

        - ada

-----------------------------------------------------------------------
%}

    dt = sensorData.barometer.time(2) - sensorData.barometer.time(1);
    
    x  = xin';
    
    At = [ 1       dt    0.5*dt^2;
           0       1         dt;
           0       0         1;];

    Ct = [ 1     0     0 ];
    
    xp = zeros(length(sensorData.barometer.time),3);
    xv = zeros(length(sensorData.barometer.time),2);
    
    for ii = 1:length(sensorData.barometer.time)
    
        % Prediction step:
        x      =   At * x;
    
        % Prediction variance propagation:
        P      =   ada.Q + At * Pin * At';
    
        % Correction step:
        S      =   Ct * P * Ct' + ada.R;
        K      =   P * Ct' /S;
        x      =   x + K*(sensorData.barometer.measures(ii) - Ct*x);
        Pout   =  (eye(3) - K*Ct) * P;

        xp(ii,:)  =   x';
        
        xv(ii,1)  =   getaltitude(xp(ii,1),ada.temp_ref, ada.p_ref);
        xv(ii,2)  =   getvelocity(xp(ii,1),xp(ii,2),ada.temp_ref, ada.p_ref);
        
        if ada.flag_apo  == false
            if xv(ii,2) < ada.v_thr
                ada.counter = ada.counter + 1;
            else
                ada.counter = 0;
            end
            if ada.counter >= ada.count_thr
            	ada.t_ada = sensorData.barometer.time(ii);
                ada.flag_apo = true;
            end
        end
    end
    
end

function h = getaltitude(p, temp_ref, p_ref)
a  = 0.0065;
n  = 9.807/(287.05*a);

h  = temp_ref / a * (1 - (p / p_ref)^(1/n));
end

function v = getvelocity(p, dpdt, temp_ref, p_ref)
a  = 0.0065;
n  = 9.807/(287.05*a);

v  = -(temp_ref * dpdt * (p / p_ref)^ (1/n)) / (a * n * p);
end

        