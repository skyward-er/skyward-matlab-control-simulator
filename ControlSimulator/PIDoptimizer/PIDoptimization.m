clear;
close all;

objectiveFcn = @simBasedObjectiveFcn;

x0=[77 5];


%% 
opts = optimoptions(@patternsearch,'PlotFcn',{@psplotbestf,@psplotfuncount});

[X1,Fval,Exitflag,Output] = patternsearch(objectiveFcn,x0,[],[],[],[],[],[],[],opts);
fprintf('The number of iterations is: %d\n', Output.iterations);
fprintf('The number of function evaluations is: %d\n', Output.funccount);
fprintf('The best function value found is: %g\n', Fval);


trapz 


