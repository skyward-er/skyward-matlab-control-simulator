function [x,P] = kalmanFilterPrediction(x_prev,dt,P_prev,ab,w,Q)
%13/11/2020  ANY QUESTIONS CAN BE DIRECTED TO ALEJANDRO MONTERO FROM SKYWARD

%-----------DESCRIPTION OF FUNCTION:------------------

%STATE SPACE ESTIMATOR (PREDICTION STEP) FOR LINEAR MOVEMENT OF
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
%           -x_prev:    1x10 VECTOR OF PREVIOUS VALUES --> 3 FIRST STATES
%                       ARE X , Y AND H, THE FOLLOWING THREE ARE VX, VY AND VZ
%                       AND THE LAST 4 ARE THE QUATERNION COMPONENTS
%           -P_prev:    10x10 MATRIX OF PREVIOUS COVARIANCE OF STATE
%           -dt:        TIME STEP
%           -t_sam:     VECTOR OF INSTANTS OF MEASUREMENT OF INPUTS
%           -ab:        VECTOR OF LINEAR ACCELERATION MEASUREMENT AT T --> 1X3
%           -w:         VECTOR OF ANGULAR VELOCITY MEASUREMENT AT T --> 1X3
%           -Q:         COVARIANCE MATRIX OF PROCESS NOISE
%
%       -OUTPUTS:
%           -x_es:      STATE ESTIMATION AT T. VECTOR WITH 10 COLUMNS
%           -P:         MATRIX OF VARIANCE OF THE STATE AT T--> IS A
%                       10 x 10 matrix
%---------------------------------------------------------------------------

    
               
A              =   [x_prev(7)^2 - x_prev(8)^2 - x_prev(9)^2 + x_prev(10)^2,                   2*(x_prev(7)*x_prev(8) + x_prev(9)*x_prev(10)),               2*(x_prev(7)*x_prev(9) - x_prev(8)*x_prev(10));
                            2*(x_prev(7)*x_prev(8) - x_prev(9)*x_prev(10)),          -x_prev(7)^2 + x_prev(8)^2 - x_prev(9)^2 + x_prev(10)^2,               2*(x_prev(8)*x_prev(9) + x_prev(7)*x_prev(10));
                            2*(x_prev(7)*x_prev(9) + x_prev(8)*x_prev(10)),                   2*(x_prev(8)*x_prev(9) - x_prev(7)*x_prev(10)),     -x_prev(7)^2 - x_prev(8)^2 + x_prev(9)^2 + x_prev(10)^2;];




a             =   (A'*ab')';                    %Rotation of the acceleration 
                                                %from body axis to inertial frame 
                                                %to use the inertial equations of motion                                                


w_mat         =   [0          w(3)    -w(2)    w(1); %Defition of the
                  -w(3)       0        w(1)    w(2); %matrix needed
                   w(2)      -w(1)     0       w(3); %for the quaternion
                  -w(1)      -w(2)    -w(3)    0];       %propagation 
               
x_dot(1:6)    =   [x_prev(4)  x_prev(5)  x_prev(6)   a+[0,0,9.81]];   %The derivatives 
                                                                    %of the position are the velocities of the previous 
                                                                    %instant and the derivatives of the
                                                                    %velocities are the acceleration
                                                                    %measured in the previous instant

x_dot(7:10)   =   0.5*(w_mat*x_prev(7:10)')';  %The derivative of the last 4 states
                                                %which are the
                                                %quaternnions, are obtained
                                                %multiplying the matrix of
                                                %angular velocities by the
                                                %previous quaternion

x             =   x_prev + dt*x_dot; %The following estimated state is 
                                        %the previous one plus the derivative 
                                        %multiplied by the time step ---> 
                                        %FORWARD EULER SIMPLE PROPAGATION

x(7:10)       =   x(7:10)/(norm(x(7:10)));               %Quaternions have 
                                                       %to be normalised at each time step
%-------------------- Jacobianv -----------------------------------------------
F              =   sparse(10,10);                %Pre-allocation of gradient 
                                                %of the derivative of the state
                                                
F(1,4)         =   1;                          
F(2,5)         =   1;
F(3,6)         =   1;                           %definition of the gradient 
                                                %for the linear dynamics --> CONSTANT
F(4,7:10)     =   2*[[x_prev(7)     x_prev(8)     x_prev(9)]*ab' ...
                    [-x_prev(8)     x_prev(7)     x_prev(10)]*ab'...
                    [-x_prev(9)    -x_prev(10)    x_prev(7)]*ab'...
                    [x_prev(10)    -x_prev(9)     x_prev(8)]*ab']; %Update
                                                                %of the gradient for the linear
                                                                %dynamics due to the rotation of
                                                                %acceleration (row 4)

F(5,7:10)     =  2*[[x_prev(8)     -x_prev(7)    -x_prev(10)]*ab'...
                    [x_prev(7)       x_prev(8)     x_prev(9)]*ab'...
                    [x_prev(10)     -x_prev(9)     x_prev(8)]*ab'...
                    [x_prev(9)       x_prev(10)   -x_prev(7)]*ab']; %Update
                                                                %of the gradient for the linear
                                                                %dynamics due to the rotation of
                                                                %acceleration (row 5)                                    
F(6,7:10)     =   2*[[x_prev(9)     x_prev(10)   -x_prev(7)]*ab'...
                    [-x_prev(10)    x_prev(9)    -x_prev(8)]*ab'...
                    [ x_prev(7)     x_prev(8)     x_prev(9)]*ab'...
                    [-x_prev(8)     x_prev(7)     x_prev(10)]*ab']; %Update
                                                                %of the gradient for the linear
                                                                %dynamics due to the rotation of
                                                                %acceleration (row 6)

F(7:10,7:10)  =   0.5*w_mat;              %Update of the gradient 
                                          %for the quaternion part of the state                                               

F             =   eye(10)+dt*F;                                                
                                               
% Covariance propagation       

P             =   F*P_prev*F' + Q;          %Prediction of the covariance 
                                            %matrix of the state 
end
                                              