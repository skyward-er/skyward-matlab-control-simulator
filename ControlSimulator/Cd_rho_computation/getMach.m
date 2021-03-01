function mach = getMach(v,z)

co   = 340.3;
alfa = -3.871e-3;

c = co + alfa*z;
mach = v./c;

end
