clear;
close all;

objectiveFcn = @simBasedObjectiveFcn;

x0=[300 100];
%% 
opts = optimoptions(@patternsearch,'PlotFcn',{@psplotbestf,@psplotfuncount});

[X1,Fval,Exitflag,Output] = patternsearch(objectiveFcn,x0,[],[],[],[],[],[],[],opts);
fprintf('The number of iterations is: %d\n', Output.iterations);
fprintf('The number of function evaluations is: %d\n', Output.funccount);
fprintf('The best function value found is: %g\n', Fval);
% %% 
% [X1,Fval,Exitflag,Output] = fminsearch(objectiveFcn,x0);
% %opts = optimoptions(@fminsearch,'PlotFcn',{@psplotbestf,@psplotfuncount});
% fprintf('The number of iterations is: %d\n', Output.iterations);
% fprintf('The best function value found is: %g\n', Fval);
% fprintf('Output quality (0 lack of convergence): %g\n',Exitflag);



