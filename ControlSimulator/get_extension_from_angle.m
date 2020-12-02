function [x] = get_extension_from_angle(alpha_degree)

% Transform the angle in radiants
alpha_rad = (alpha_degree*pi)/180;

% Obtain delta_S given the servomotor angle
a = -9.43386/1000;
b = 19.86779/1000;
delta_S = a*(alpha_rad^2) + b*alpha_rad;
% delta_S [m^2] = (-9.43386 * alpha^2 + 19.86779 * alpha) * 10^(-3), with alpha in [rad] from 0 to 0.89 rad.

% Obtain aerobrakes extension given the delta_S
a = -1.04034;
b = 0.30548;
x = (-b + sqrt(b^2 + 4*a*delta_S)) / (2*a);
% A [m^2] = -1.04034 * x^2 + 0.30548 * x, with x in [m] from 0 to 0.03866 m;

end