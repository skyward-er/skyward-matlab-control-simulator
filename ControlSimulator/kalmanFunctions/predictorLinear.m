function [x,P] = predictorLinear(x_prev,dt,P_prev,ab,q,Q)
%24/02/2021  ANY QUESTIONS CAN BE DIRECTED TO ALEJANDRO MONTERO FROM SKYWARD

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
%           -x_prev:    1x6 VECTOR OF PREVIOUS VALUES --> 3 FIRST STATES
%                       ARE X , Y AND H, THE FOLLOWING THREE ARE VX, VY AND VZ
%                       
%           -P_prev:    6x6 MATRIX OF PREVIOUS COVARIANCE OF STATE
%           -dt:        TIME STEP
%           -ab:        VECTOR OF LINEAR ACCELERATION MEASUREMENT AT T --> 1X3
%           -q:         QUATERNION AT PREVIOUS TIME STEP; 1x4; [q;q0] format
%           -Q:         COVARIANCE MATRIX OF PROCESS NOISE
%
%       -OUTPUTS:
%           -x_es:      STATE ESTIMATION AT T. VECTOR WITH 6 COLUMNS
%           -P:         MATRIX OF VARIANCE OF THE STATE AT T--> IS A
%                       6 x 6 matrix
%---------------------------------------------------------------------------

    
               
A              =   [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2,                   2*(q(1)*q(2) + q(3)*q(4)),               2*(q(1)*q(3) - q(2)*q(4));
                            2*(q(1)*q(2) - q(3)*q(4)),          -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,               2*(q(2)*q(3) + q(1)*q(4));
                            2*(q(1)*q(3) + q(2)*q(4)),                   2*(q(2)*q(3) - q(1)*q(4)),     -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2;];




a             =   (A'*ab')';                    %Rotation of the acceleration 
                                                %from body axis to inertial frame 
                                                %to use the inertial equations of motion                                                

               
x_dot(1:6)    =   [x_prev(4)  x_prev(5)  x_prev(6)   a+[0,0,9.81]];   %The derivatives 
                                                                    %of the position are the velocities of the previous 
                                                                    %instant and the derivatives of the
                                                                    %velocities are the acceleration
                                                                    %measured in the previous instant


x             =   x_prev + dt*x_dot; %The following estimated state is 
                                        %the previous one plus the derivative 
                                        %multiplied by the time step ---> 
                                        %FORWARD EULER SIMPLE PROPAGATION

%-------------------- Jacobianv -----------------------------------------------
F              =   sparse(10,10);                %Pre-allocation of gradient 
                                                %of the derivative of the state
                                                
F(1,4)         =   1;                          
F(2,5)         =   1;
F(3,6)         =   1;                           %definition of the gradient 
                                                %for the linear dynamics --> CONSTANT
                                              

F             =   eye(6)+dt*F;                                                
                                               
% Covariance propagation       

P             =   F*P_prev*F' + Q;          %Prediction of the covariance 
                                            %matrix of the state 
end