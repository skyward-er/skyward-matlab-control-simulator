function alfa = getAngle(S)
a = -9.43386*2e-3;
b =  19.86779e-3 ;
alfa = -b/a + sqrt(b^2+2*a*S)/a;
end
