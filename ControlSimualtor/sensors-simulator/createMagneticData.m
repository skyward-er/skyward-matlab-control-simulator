function [outX,outY,outZ] = createMagneticData(magneticDensity, q0, q1, q2, q3)
%CREATE_MAGNETIC_DATA Uses quaternions to rotate the earth magnetic feel
%and therefor create magnetic data from a magnetometer
%
%   Inputs:
%   magneticDensity: density of the earth magnetic field (between 250 to
%   650 mgauss) https://en.wikipedia.org/wiki/Earth%27s_magnetic_field
%   q0, q1, q2, q3: quaternion parameter to rotate
%   
%   Outputs:
%   outX,outY,outZ: Output coordinates after rotation

% earth magnetic
m = [magneticDensity, 0, 0];

% rotate
xNew = quatrotate([q0, q1, q2, q3], m);

outX=xNew(1);
outY=xNew(2);
outZ=xNew(3);
end

