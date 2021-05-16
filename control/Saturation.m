function [u, saturation] = Saturation(u, umin, umax)

if (u < umin)
    u = umin;
    saturation = true;   
elseif (u > umax)
    u = umax;
    saturation = true;    
else 
    saturation = false;    
end
end
