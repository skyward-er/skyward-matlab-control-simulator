function [fix, n_sat] = gpsFix(accel)
% Routine to simulate the data acquisition from the sensors

% Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept
% email: alessandro.delduca@skywarder.eu
% Revision date: 18/03/2021
if norm(accel)/9.81 < 2
    if rand(1)>0.05
        fix   = true;
        n_sat = 4;
    else
        fix   = 0;
        n_sat = 0;
    end
else
    if rand(1)>0.95
        fix   = true;
        n_sat = 4;
    else
        fix   = 0;
        n_sat = 0;
    end
end
end

