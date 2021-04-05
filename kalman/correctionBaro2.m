function [x,P,p_pred,y_res] = correctionBaro2(x_pred, P_pred, p, p_pred, sigma_p)

% Author: Alejandro Montero
% Co-Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@skywarder.eu
% email: alejandro.montero@skywarder.eu, alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

%-----------DESCRIPTION OF FUNCTION:------------------

%STATE SPACE ESTIMATOR (CORRECTION STEP FOR BAROMETER) FOR LINEAR MOVEMENT OF
%ROCKET AND ATTITUDE DYNAMICS
%THE DYNAMIC SYSTEM DESCRIPTION IS:
%       x' = f(x,u) + w         F=df/dx --> F IS THE GRADIENT OF f
%                                           EVALUATED AT EACH ESTIMATION 
%                               w is process noise --> Q IS ITS COVARIANCE
%       z  = h(x,u) + v         H=dh/dx --> H IS THE GRADIENT OF h
%                                           EVALUATED AT EACH ESTIMATION
%                               v is measurement noise --> R IS ITS
%                               COVARIANCE
%       -INPUTS:
%           -x_pred:    1x6 VECTOR OF PREDICTED VALUES --> 3 FIRST STATES
%                       ARE X , Y AND H, THE FOLLOWING THREE ARE VX, VY AND VZ
%                       
%           -P_pred:    6x6 MATRIX OF PREDICTED COVARIANCE OF STATE
%           -h_sam:     MEASUREMENT OF ALTITUDE FROM BAROMETER AT TIME T --> 1x1
%           -sigma_h:   VARIANCE OF THE BAROMETER
%
%       -OUTPUTS:
%           -x_es:      STATE ESTIMATION CORRECTED AT T. VECTOR WITH 6 COLUMNS
%           -P:         MATRIX OF VARIANCE OF THE STATE AT T, CORRECTED--> IS A
%                       6 x 6 matrix
%           -y_res:     VECTOR OF DIFFERENCES BETWEEN THE CORRECTED ESTIMATION 
%                       OF THE OUTPUT AND THE MEASSURE; ONLY FOR CHECKING
%                       --> 1x1
%---------------------------------------------------------------------------
temp_ref       =   288.1500;
p_ref          =   101325;
dt             =   0.05;

threshold      =   10e-11;

H              =   [0 0 1 0 0 0;
                    0 0 0 0 0 1;];
                
R              =   diag([sigma_p^2, sigma_p^2]);

z              =   H*x_pred';

S              =   H*P_pred*H'+R;                %Matrix necessary for the correction factor

h              =   getaltitude(p, temp_ref, p_ref);

dpdt           =   grad(p,p_pred,dt);

vz             =   getvelocity(p, dpdt, temp_ref, p_ref);

y              =   [-h;-vz];
   if cond(S) > threshold 
       
       e       =   y - z;
       K       =   P_pred*H'/S;                   %Kalman correction factor

       x       =   x_pred + (K*e)';               %Corrector step of the state
       

       P       =   (eye(6) - K*H)*P_pred;        %Corrector step of the state covariance
    else
       x       =   x_pred;
       P       =   P_pred;
   end
   
   
z_corr         =   H*x';                          %Corrected output expectation

y_res          =   y - z_corr;
   
p_pred         =   p;
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

function dudt = grad(u, u_prev, dt)
dudt = (u - u_prev)/dt;
end