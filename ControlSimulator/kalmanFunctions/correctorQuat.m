function [x_c,P_c,e,z]=correctorQuat(x_pred,P_pred,mag_sam,sigma_mag)
%15/12/2020  ANY QUESTIONS CAN BE DIRECTED TO ALEJANDRO MONTERO FROM SKYWARD

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

%--------------------------------------------------------------------------
% Computation of the output equation and innovation term of the filter
A       = [x_pred(1)^2 - x_pred(2)^2 - x_pred(3)^2 + x_pred(4)^2,               2*(x_pred(1)*x_pred(2) + x_pred(3)*x_pred(4)),                 2*(x_pred(1)*x_pred(3) - x_pred(2)*x_pred(4));
                 2*(x_pred(1)*x_pred(2) - x_pred(3)*x_pred(4)),      -x_pred(1)^2 + x_pred(2)^2 - x_pred(3)^2 + x_pred(4)^2,                2*(x_pred(2)*x_pred(3) + x_pred(1)*x_pred(4)) ;
                 2*(x_pred(1)*x_pred(3) + x_pred(2)*x_pred(4)),               2*(x_pred(2)*x_pred(3) - x_pred(1)*x_pred(4)),       -x_pred(1)^2 - x_pred(2)^2 + x_pred(3)^2 + x_pred(4)^2];

z       = A*[1;0;0];       %Magnetic vector in the body axis (estimated)

z_mat   = [ 0      -z(3)   z(2);
           z(3)     0     -z(1);
          -z(2)     z(1)      0;];        %Matrix needed to obtain the derivative H
%-------------------------------------------------------------------------
% Computation of the derivative matrix of the output equation (H) and
% kalman gain (K)
H       = [z_mat zeros(3,3)];

S       = H*P_pred*H'+R;

K       = P_pred*H'/S;
%------------------------------------------------------------------------
% Correction
% Innovation computation
e       = mag_sam' - z;     %Difference between measured (mag_sam) and 
                            %estimated (z) magnetic vectors
% Error correction computation (delta_x)
delta_x = (K*e)';

% Since the definition of the quaternions is [q_vec;q4] but the function
% quatmultiply uses [q1;q_vec], it is necessary to assemble q2 changing the
% order of the components
q2      = [x_pred(4),x_pred(1:3)];

% Assembly of the error quaternion from Euler angles alpha --> delta_x(1:3)
r       = [sqrt(1-0.25*delta_x(1:3)*delta_x(1:3)'), 0.5*delta_x(1:3)];

u       = quatmultiply(r,q2);              %The correction of the quaternion
                                           %is done through a multiplication 
                                           %(from which the multiplicative 
                                           %filter gets its name)
                                           
% Re-assembly of the previously used notation
aux     = u(1);
u(1:3)  = u(2:4);
u(4)    = aux;
% Assigment of the new quaternion to the state vector
x_c(1:4)= u/norm(u);                       %Re-normalisation of the quaternion 
                                           %to avoid issues
x_c(5:7)= x_pred(5:7)+delta_x(4:6);        %Correction of the bias only with a sum

P_c     = (eye(6) - K*H)*P_pred*(eye(6) - K*H)'+ K*R*K';           %Correction of the covariance matrix
end