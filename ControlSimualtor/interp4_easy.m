function [V] = interp4_easy(X1, X2, X3, X4, F, x1, x2, x3, x4)

%{ 

interp4_easy - This function interpolate with nearest-neighbor method a R4->R function

INPUTS:    
            - F = F(x1...x4), [NxMxLxI] Matrix that discretize the function;
            - X1,...,X4: column vectors of the variables the function depends on;
            - (x1....x4): single values for each variable the interpolation is wanted.

OUTPUTS:    
            - V interpolated coefficients.

Author: Ruben Di Battista
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: ruben.dibattista@skywarder.eu
April 2014; Last revision: 25.IV.2014

%}

M = {X1 X2 X3 X4};
m = [x1 x2 x3 x4];

index = zeros(4,1);

for i = 1:4
    [~, index(i)] = min(abs(M{i} - m(i)));
end

V = F(index(1),index(2),index(3),index(4));
% V = interpn(X1,X2,X3,X4,F,x1,x2,x3,x4);



end

