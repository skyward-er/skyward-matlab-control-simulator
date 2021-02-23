function [x,P,y_res] = correctionMagnetometer(x_pred,P_pred,mag_sam,sigma_mag)
%13/11/2020  ANY QUESTIONS CAN BE DIRECTED TO ALEJANDRO MONTERO FROM SKYWARD

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
%           -x_pred:      1x10 VECTOR OF PREDICTED VALUES --> 3 FIRST STATES
%                         ARE X , Y AND H, THE FOLLOWING THREE ARE VX, VY AND VZ
%                         AND THE LAST 4 ARE THE QUATERNION COMPONENTS
%           -P_pred:      10x10 MATRIX OF PREDICTED COVARIANCE OF STATE
%           -mag_sam:     MEASUREMENT OF THE MAGNETIC FIELD DIRECTION FROM 
%                         MAGNETOMETER AT TIME T --> 3x1
%           -sigma_mag:   VARIANCE OF THE MAGNETOMETER
%
%       -OUTPUTS:
%           -x_es:        STATE ESTIMATION CORRECTED AT T. VECTOR WITH 10 COLUMNS
%           -P:           MATRIX OF VARIANCE OF THE STATE AT T, CORRECTED--> IS A
%                         10 x 10 matrix
%           -y_res:       VECTOR OF DIFFERENCES BETWEEN THE CORRECTED ESTIMATION 
%                         OF THE OUTPUT AND THE MEASSURE; ONLY FOR CHECKING
%                         --> 3x1
%---------------------------------------------------------------------------
threshold      =   10e-3;

                                                
R              =  sigma_mag^2*eye(3);

               
A   =   [x_pred(7)^2 - x_pred(8)^2 - x_pred(9)^2 + x_pred(10)^2,               2*(x_pred(7)*x_pred(8) + x_pred(9)*x_pred(10)),                 2*(x_pred(7)*x_pred(9) - x_pred(8)*x_pred(10));
                 2*(x_pred(7)*x_pred(8) - x_pred(9)*x_pred(10)),      -x_pred(7)^2 + x_pred(8)^2 - x_pred(9)^2 + x_pred(10)^2,                2*(x_pred(8)*x_pred(9) + x_pred(7)*x_pred(10)) ;
                 2*(x_pred(7)*x_pred(9) + x_pred(8)*x_pred(10)),               2*(x_pred(8)*x_pred(9) - x_pred(7)*x_pred(10)),       -x_pred(7)^2 - x_pred(8)^2 + x_pred(9)^2 + x_pred(10)^2];

z=A*[1;0;0];

H   =   2*[zeros(1,6) x_pred(7) -x_pred(8)  -x_pred(9)   x_pred(10) ;
           zeros(1,6) x_pred(8)  x_pred(7)  -x_pred(10) -x_pred(9) ;
           zeros(1,6) x_pred(9)  x_pred(10)  x_pred(7)   x_pred(8) ;];
 
e=mag_sam'-z;

S   =   H*P_pred*H'+R;
    if cond(S)>threshold
        K       = P_pred*H'/S;
        x       = x_pred + (K*e)';
        x(7:10) = x(7:10)/norm(x(7:10));
        P       = (eye(10) - K*H)*P_pred;
    else

        x=x_pred;
        P=P_pred;
    end

A       =   [x(7)^2 - x(8)^2 - x(9)^2 + x(10)^2,               2*(x(7)*x(8) + x(9)*x(10)),                 2*(x(7)*x(9) - x(8)*x(10));
                     2*(x(7)*x(8) - x(9)*x(10)),      -x(7)^2 + x(8)^2 - x(9)^2 + x(10)^2,                 2*(x(8)*x(9) + x(7)*x(10));
                     2*(x(7)*x(9) + x(8)*x(10)),               2*(x(8)*x(9) - x(7)*x(10)),       -x(7)^2 - x(8)^2 + x(9)^2 + x(10)^2];

z_corr  =    A*[1;0;0];

y_res   =    mag_sam'-z_corr;

    end