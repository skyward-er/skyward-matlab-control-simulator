function y = simBasedObjectiveFcn(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Close figures

save('interationData.mat','x')
run( '../start_simulation.m');
load('interationData.mat');
y=(x(1)^2+x(2)^2);
end

