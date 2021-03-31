function y = exp_mean(u,alfa)
y = 0;
for i = 1:length(u)
    y = alfa*u(i) + (1 - alfa)*y;
end
end

 
