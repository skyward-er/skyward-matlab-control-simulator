function [x_pred,P_pred]=predictorQuat(x,P,w,dt,Q)

% Author: Alejandro Montero
% Co-Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@skywarder.eu
% email: alejandro.montero@skywarder.eu, alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

%-----------DESCRIPTION OF FUNCTION:------------------

%STATE SPACE ESTIMATOR (PREDICTION STEP) FOR ATTITUDE DYNAMICS
%THE DYNAMIC SYSTEM DESCRIPTION IS:
%       x' = f(x,u) + w         F=df/dx --> F IS THE GRADIENT OF f
%                                           EVALUATED AT EACH ESTIMATION 
%                               w is process noise --> Q IS ITS COVARIANCE
%       -INPUTS:
%           -x:         1x7 VECTOR OF PREVIOUS VALUES --> 4 FIRST STATES
%                       ARE QUATERNION AND THE FOLLOWING THREE ARE BIASES
%           -P:         6x6 MATRIX OF PREVIOUS COVARIANCE OF STATE
%                       ONLY 6 BECAUSE OF THE SIMPLIFICATION IN THE ERROR
%                       QUATERNION
%           -dt:        TIME STEP
%           -w:         VECTOR OF ANGULAR VELOCITY MEASUREMENT AT T --> 1X3
%           -Q:         COVARIANCE MATRIX OF PROCESS NOISE
%
%       -OUTPUTS:
%           -x_pred:      STATE ESTIMATION AT T. VECTOR WITH 7 COLUMNS
%           -P_pred:      MATRIX OF VARIANCE OF THE STATE AT T--> IS A
%                         6 x 6 matrix
%---------------------------------------------------------------------------
q_prev      = x(1:4);                                 %Definition of the previous quaternion
beta_prev   = x(5:7);                                 %Definition of the previous bias

omega       = w - beta_prev;                            %Computation of real w (no bias)
% omega = w;
%------------------------------------------------------------------------
%COMPUTATION OF PROPAGATION MATRIX FOR QUATERNION.
omega_mat   = [ 0      -omega(3)   omega(2);
               omega(3)  0       -omega(1);
               -omega(2)  omega(1)   0;];
Omega       = [ -omega_mat  omega';
                -omega      0];
%-------------------------------------------------------------------------
%PROPAGATION OF QUATERNION AND BIAS (no dynamics) +  ASSIGMENT TO X          
q_pred      = (eye(4) + 0.5*Omega*dt)*q_prev';

q_pred      = q_pred/norm(q_pred);

x_pred(1:4) = q_pred';

x_pred(5:7) = beta_prev;
%-------------------------------------------------------------------------
%PROPAGATION OF COVARIANCE ERROR MATRIX P. DISCRETE APPROXIMATION.
G           = [ -eye(3)       zeros(3,3);
                 zeros(3,3)   eye(3)];
F           = [ -eye(3)       -eye(3)*dt;
                 zeros(3,3)   eye(3)];
P_pred      = F*P*F'+G*Q*G';
end