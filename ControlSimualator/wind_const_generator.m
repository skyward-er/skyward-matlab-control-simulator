function [uw, vw, ww, Az] = wind_const_generator(AzMin, AzMax, ElMin, ElMax, MagMin, MagMax)
%{

wind_const_generator - function that generates wind components in NED axes

INPUTS:
            - AzMin, Minimum angle of Azimuth from North;
            - AzMax, Maximum angle of Azimuth from North;
            - ElMin, Minimum angle of Elevation;
            - ElMax, Maximum angle of Elevatiom.

OUTPUTS:
            - uw, wind component along x;
            - vw, wind component along y;
            - ww, wind component along z;
            - Az, angle of Azimuth from North.

Author: Ruben Di Battista
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu
April 2014; Last revision: 25.IV.2014

%}

% Generating random values for orientation and magnitude
Az = AzMin + (AzMax-AzMin)*rand;
El = ElMin + (ElMax-ElMin)*rand;
Mag = MagMin + (MagMax-MagMin)*rand;

% Random Wind Vector
R = Mag*angle2dcm(Az,El,0,'ZYX');
R(abs(R) < 1e-4) = 0;

uw = R(1,1);
vw = R(1,2);
ww = R(1,3);

if abs(uw) < 1e-3
    uw = 0;
end

if abs(vw) < 1e-3
    vw = 0;
end

if abs(ww) < 1e-3
    ww = 0;
end




end

