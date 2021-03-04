function [x_ada, P_ada, flag_ADA, t_ADA, count_ADA]   =  run_ADA(x_ada, P_ada, h_baro, t_baro, Q_ada, R_ada, N, count_ADA, flag_ADA, t_ADA)

% Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@skywarder.eu
% email: alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

%-----------DESCRIPTION OF FUNCTION:------------------
% This function simulate the Apogee Detection Algorithm (ADA), the
% algorithm is a 3rd order linear kalman filter with this states:
% - x_ada(1): Vertical position
% - x_ada(2): Vertical velocity
% - x_ada(3): Vertical acceleration
% The filter model is a constant acceleration point mass, with trial
% computation we will end up with this discrete time system:
% - p(t+1) = pt + v*dt +0.5*a*dt^2
% - v(t+1) = vt + a*dt
% - a(t+1) = const
% with observation on the first state y(t) = p(t)
%      INPUTS: 
%         - x_ada:      1x3 VECTOR OF PREVIOUS VALUES --> p,v,a
%
%         - P_ada:      3x3 MATRIX OF PREVIOUS COVARIANCE OF STATE
%
%         - t_baro:     VECTOR OF TIME INSTANTS AT WHICH THE BAROMETER TOOK 
%                       SAMPLES INSIDE THE 0.1 INTEGRATION PERIOD. SINCE IT RUNS
%                       AT 20 HZ, THIS IS A 2x1 VECTOR. s
%
%         - h_baro:     CORRESPONDING ALTITUDE MEASUREMENTS FROM THE BAROMETER.
%                       2x1. m
%
%         - Q_ada:      COVARIANCE MATRIX OF PROCESS NOISE. 3x3   
%
%         - R_ada:      COVARIANCE MATRIX OF OBSERVATION NOISE. 1x1 
%
%         - N:          INTEGRATION WINDOWS FOR THE PREDICTION
%
%         - count_ADA:  COUNTER OF THE NUMBER OF CONSECUTIVES APOGEE
%                      	DETECTED
%
%         - flagADA:    BOOLEAN VARIABLE THAT BECOME TRUE AFTER THE APOGEE IS
%                       DETECTED count_threshold TIMES
%
%         - t_ADA:      TIME OF THE PREDICTED APOGEE

%       OUTPUTS:

%         - x_ada:      FILTERED STATES OF THE ADA 
% 
%         - P_ada:      PROPAGATED COVARIANCE MATRIX FOR EACH OF THE ESTIMATION TIME
%                       INSTANTS. 3x3x2
%
%         - flagADA:    BOOLEAN VARIABLE THAT BECOME TRUE AFTER THE APOGEE IS
%                       DETECTED count_threshold TIMES
%
%         - t_ADA:      TIME OF THE PREDICTED APOGEE 
%
%         - count_ADA:  COUNTER OF THE NUMBER OF CONSECUTIVES APOGEE
%                      	DETECTED
% -----------------------------------------------------------------------
    count_threshold = 5;

    dt = t_baro(2) - t_baro(1);
    x  = x_ada';
    At = [ 1    dt 0.5*dt^2;
           0     1     dt;
           0     0     1;];

    Ct = [ 1     0     0 ];
    
    for ii = 1:length(t_baro)
    
        % Prediction step:
        x      =   At*x;
    
        % Prediction variance propagation:
        P_ada  =   Q_ada + At*P_ada*At';
    
        % Correction step:
        S      =   Ct*P_ada*Ct' + R_ada;
        K      =   P_ada*Ct'/S;
        x      =   x + K*(h_baro(ii) - Ct*x);
        P_ada  =  (eye(3) - K*Ct)*P_ada;

        x_ada(ii,:)  =   x';
    % Prediction N state ahead and check if the apogee is reached
        if flag_ADA == false
            xapo = x;
%            xapo = (At^N)*xapo;
            if xapo(2) < 0 
                count_ADA = count_ADA + 1;
            else
                count_ADA = 0;
            end
            if count_ADA >= count_threshold
            	t_ADA = t_baro(ii);
                flag_ADA = true;
            end
        end
    end
end

    
        