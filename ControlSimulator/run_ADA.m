function [x_ada,P_ada,flagADA,t_ADA]   =  run_ADA(x_ada, P_ada, h_baro, t_baro, Q_ada, R_ada)
    dt = t_baro(2) - t_baro(1);
    x  = x_ada';
    At = [ 1    dt 0.5*dt^2;
           0     1     dt;
           0     0     1;];

    Ct = [ 1     0     0 ];
for ii = 1:length(t_baro)
    % Prediction step: free fall model
    x      =   At*x;
    % Prediction variance propagation 
    P_ada  =   Q_ada + At*P_ada*At';
    % Correction step: altitude measure from the barometer
    S      =   Ct*P_ada*Ct' + R_ada;
    K      =   P_ada*Ct'/S;
    x      =   x + K*(h_baro(ii) - Ct*x);
    P_ada  =  (eye(3) - K*Ct)*P_ada;
    
    x_ada(ii,:)  =   x';
    % Prediction N state ahead and check if the apogee is reached
end
flagADA = false;
t_ADA   = false;
end
