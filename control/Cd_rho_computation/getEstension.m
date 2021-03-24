function x = getEstension(delta_S)

% A [m^2] = -1.04034 * x^2 + 0.30548 * x, with x in [m] from 0 to 0.03866 m;
% A è il delta_S o l'area totale?
a = -1.04034;
b = 0.30548;

% Obtain aerobrakes extension given the delta_S
x = (-b + sqrt(b^2 + 4*a*delta_S)) / (2*a);

end


