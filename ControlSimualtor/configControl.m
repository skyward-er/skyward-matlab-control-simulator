%% CONTROL PARAMETER
data=load('ControllerLookup.mat');
settings.control.all_Vz=data.all_Vz;
settings.control.all_alpha=data.all_alpha;
settings.control.hightInterval=data.hightInterval;

settings.control.alphaDefault=data.all_alpha(1);

settings.control.start1state=230; % velocity when control should start
settings.control.start2state=100; % velocity when control should end and alpha keeps last control value

settings.control.dt=0.1;