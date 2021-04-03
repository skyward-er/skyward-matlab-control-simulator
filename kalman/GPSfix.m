function flagGPS_fix=GPSfix(accel)
%Author: Alejandro Montero
%email: alejandro.montero@skywarder.eu

%This function computes the fix of the GPS based on the acceleration norm.
%If the acceleration modulus is bigger than a certain threshold, the
%probablity of losing the GPS fix is bigger. Otherwise, the fix should be
%available most of the time.
%INPUT:
%   -accel: Acceleration vector at the computation time. 3x1.
%
%OUPUT:
%   -flagGPS_fix: binary variable that decides whether the GPS can be
%   trusted or not
g=9.81;
if norm(accel)/g>2
    if rand(1)>0.8
        flagGPS_fix=1;
    else
        flagGPS_fix=0;
    end
else
    if rand(1)>0.01
        flagGPS_fix=1;
    else
        flagGPS_fix=0;
    end
end

end