function [V] = InterpAero(settings, AoAvec, Mvec, Bvec, hvec, cvec, CmatF, CmatE, AoA, M, B, h, c, t)

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

tb = settings.tb;
T = {AoAvec Mvec Bvec hvec};
inst = [AoA M B h];

index = zeros(4,1);

for i = 1:4
    [~, index(i)] = min(abs(T{i} - inst(i)));
end



if c == 0
    VE = CmatE(index(1), index(2), index(3), index(4), 1);
else
    c_cmp = cvec(c > cvec); 
    n0 = length(c_cmp);
    n1 = n0 + 1;
    c0 = c_cmp(end);
    c1 = cvec(n1);
    C0 =  CmatE(index(1), index(2), index(3), index(4), n0);
    C1 =  CmatE(index(1), index(2), index(3), index(4), n1);
    VE = C1 + ((C1 - C0)/(c1 - c0))*(c - c1);
end

if t <= tb
    VF = CmatF(index(1), index(2), index(3), index(4));
    V =  t/tb*(VE-VF)+VF;
else 
    V = VE;
end



