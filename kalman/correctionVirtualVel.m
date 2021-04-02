function [x,P,y_res] = correctionVirtualVel(x_pred,P_pred,v_sam,sigma_vv)

% Author: Alejandro Montero
% Co-Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@skywarder.eu
% email: alejandro.montero@skywarder.eu, alessandro.delduca@skywarder.eu
% Release date: 02/04/2021

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
%           -v_sam:     MEASUREMENT OF VERTICAL VELOCITY FROM BAROMETER AT 
%                       TIME T (ESTIMATED WITH BACKWARD DIFFERENCES) --> 1x1
%           -sigma_vv:  VARIANCE OF THE VIRTUAL VELOCITY
%
%       -OUTPUTS:
%           -x_es:      STATE ESTIMATION CORRECTED AT T. VECTOR WITH 6 COLUMNS
%           -P:         MATRIX OF VARIANCE OF THE STATE AT T, CORRECTED--> IS A
%                       6 x 6 matrix
%           -y_res:     VECTOR OF DIFFERENCES BETWEEN THE CORRECTED ESTIMATION 
%                       OF THE OUTPUT AND THE MEASSURE; ONLY FOR CHECKING
%                       --> 1x1
%---------------------------------------------------------------------------
threshold      =   10e-3;
H              =   sparse(1,6);                %Pre-allocation of gradient 
                                                %of the output function
                                                
R              =   sigma_vv^2;

z              =   x_pred(6);

H(6)           =   1;                            %Update of the matrix H 

S              =   H*P_pred*H'+R;                %Matrix necessary for the correction factor

   if cond(S) > threshold 
       
       e       =   v_sam - z;
       K       =   P_pred*H'/S;                   %Kalman correction factor

       x       =   x_pred + (K*e)';               %Corrector step of the state
       

       P       =   (eye(6) - K*H)*P_pred;        %Corrector step of the state covariance
    else
       x       =   x_pred;
       P       =   P_pred;
   end
   
   
z_corr         =   x(6);                          %Corrected output expectation

y_res          =   v_sam - z_corr;
   
    end