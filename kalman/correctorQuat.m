function [x_c,P_c,e,z]=correctorQuat(x_pred,P_pred,mag_sam,sigma_mag,mag_NED)

% Author: Alejandro Montero
% Co-Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@skywarder.eu
% email: alejandro.montero@skywarder.eu, alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

%-----------DESCRIPTION OF FUNCTION:------------------

%STATE SPACE ESTIMATOR (CORRECTION STEP) FOR ATTITUDE DYNAMICS
%       -INPUTS:
%           -x_pred:    1x7 VECTOR OF PREDICTED VALUES --> 4 FIRST STATES
%                       ARE QUATERNION AND THE FOLLOWING THREE ARE BISES
%           -P_pred:    6x6 MATRIX OF PREDICTED COVARIANCE OF STATE
%                       ONLY 6 BECAUSE OF THE SIMPLIFICATION IN THE ERROR
%                       QUATERNION
%           -dt:        TIME STEP
%           -w:         VECTOR OF ANGULAR VELOCITY MEASUREMENT AT T --> 1X3
%           -Q:         COVARIANCE MATRIX OF PROCESS NOISE
%
%       -OUTPUTS:
%           -x_c:       STATE CORRECTION AT T. VECTOR WITH 7 COLUMNS
%           -P_c:       MATRIX OF VARIANCE CORRECTED AT T--> IS A
%                       6 x 6 matrix
%---------------------------------------------------------------------------
% Computation of the covariance matrix of the noise
R       =  sigma_mag^2*eye(3);
mag_sam =  mag_sam/norm(mag_sam);

%--------------------------------------------------------------------------
% Computation of the output equation and innovation term of the filter
A       = [x_pred(1)^2 - x_pred(2)^2 - x_pred(3)^2 + x_pred(4)^2,               2*(x_pred(1)*x_pred(2) + x_pred(3)*x_pred(4)),                 2*(x_pred(1)*x_pred(3) - x_pred(2)*x_pred(4));
                 2*(x_pred(1)*x_pred(2) - x_pred(3)*x_pred(4)),      -x_pred(1)^2 + x_pred(2)^2 - x_pred(3)^2 + x_pred(4)^2,                2*(x_pred(2)*x_pred(3) + x_pred(1)*x_pred(4)) ;
                 2*(x_pred(1)*x_pred(3) + x_pred(2)*x_pred(4)),               2*(x_pred(2)*x_pred(3) - x_pred(1)*x_pred(4)),       -x_pred(1)^2 - x_pred(2)^2 + x_pred(3)^2 + x_pred(4)^2];

z       = A*mag_NED;                   %Magnetic vector in the body axis (estimated)

z_mat   = [ 0      -z(3)   z(2);
           z(3)     0     -z(1);
          -z(2)     z(1)   0;];        %Matrix needed to obtain the derivative H
%-------------------------------------------------------------------------
% Computation of the derivative matrix of the output equation (H) and
% kalman gain (K)
H       = [z_mat zeros(3,3)];

S       = H*P_pred*H'+R;

K       = P_pred*H'/S;
%------------------------------------------------------------------------
% Correction
% Innovation computation
% mag_sam
% z
e       = mag_sam' - z;       %Difference between measured (mag_sam) and 
                              %estimated (z) magnetic vectors
delta_x = (K*e)';             % Error correction computation (delta_x)

r       = [0.5*delta_x(1:3), sqrt(1-0.25*delta_x(1:3)*delta_x(1:3)')];

u       = quatProd(r',x_pred(1:4)')';      %The correction of the quaternion
                                           %is done through a multiplication 
                                           %(from which the multiplicative 
                                           %filter gets its name)
                                           
% Assigment of the new quaternion to the state vector
x_c(1:4)= u/norm(u);                       %Re-normalisation of the quaternion 
                                           %to avoid issues
x_c(5:7)= x_pred(5:7)+delta_x(4:6);        %Correction of the bias only with a sum

P_c     = (eye(6) - K*H)*P_pred*(eye(6) - K*H)'+ K*R*K';           %Correction of the covariance matrix
end

function quat = quatProd( quat1, quat2 )
%	Calculates the Hemiltonian product between two quaternions
%
%   quat = quatProd( quat1, quat2 )
%
%   This function compute the hamiltonian product between two quaternions
%   written as 4-by-1 vectors with the scalar component as the fourth
%   element of the vector.
%
%   References:
%	[1] Markley, F. Landis. "Attitude error representations for Kalman filtering." 
%       Journal of guidance control and dynamics 26.2 (2003): 311-317.

qv1 = quat1(1:3);
qs1 = quat1(4);

qv2 = quat2(1:3);
qs2 = quat2(4);

quat = [qs1 * qv2 + qs2 * qv1 - cross( qv1, qv2 ) ;
              qs1 * qs2 - dot( qv1, qv2 )        ];

quat = quat / norm(quat);
          
end