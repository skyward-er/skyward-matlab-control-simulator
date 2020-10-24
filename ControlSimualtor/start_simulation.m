%{

START_SIMULATION - this is the main script; it runs the simulation that has been chosen in config.m

Author: Francesco Colombi
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: francesco.colombi@skywarder.eu
Release date: 16/04/2016

%}

close all
clear 
clc

path = genpath(pwd);
addpath(path);

%% LOAD DATA
run('config.m');

%% START THE CHOSEN SIMULATION
% T = vector of time used by ODE, [s] also for Tf Ta
% Y = State = ( x y z | u v w | p q r | q0 q1 q2 q3 | thetax thetay thetaz | ) also for Ya,Yf corresponding to T

tic

[Yf,Tf] = std_run_control(settings);

toc

% %% DATA-PRINTING
% 
% Na = length(Ya(:,1));
% 
% % POSITIONS
% xa = Ya(:,1);
% ya = Ya(:,2);
% za = -Ya(:,3);
% Xa = [xa, ya, za];
% 
% x = Yf(:,1);
% y = Yf(:,2);
% z = -Yf(:,3);
% 
% % VELOCITIES
% ua = Ya(:,4);
% va = Ya(:,5);
% wa = -Ya(:,6);
% Va = [ua, va, wa];
% 
% % MAXIMUM POSITIONS, VELOCITIES AND ACCELERATION
% load('ascent_plot.mat');
% A = data_ascent.accelerations.body_acc;
% 
% % pre-allocation
% abs_X = zeros(Na, 1); abs_V = abs_X; abs_A = abs_X;
% 
% % determine the norm of every row element
% for k = 1:Na
%     abs_X(k) = norm(Xa(k, :));
%     abs_V(k) = norm(Va(k, :));
%     abs_A(k) = norm(A(k, :));
% end
% 
% [max_dist, imax_dist] = max(abs_X);
% [max_v, imax_v] = max(abs_V);
% [max_a, imax_a] = max(abs_A);
% iexit = find(abs_X <= settings.lrampa);  % checking where the missile is undocked from the hook of the launch pad
% iexit = iexit(end);
% 
% % TEMPERATURE AND MACH NUMBER
% 
% % pre-allocation
% M = zeros(Na, 1); Tamb = M;
% 
% for i = 1:Na
%     [Tamb(i), a, ~, ~] = atmosisa(za(i));
%     M(i) = abs_V(i)/a;
% end
% 
% % determine the maximum Mach number
% [max_M, imax_M] = max(M);
% 
% % determine latitude and longitude of the landing point
% [lat_LP, lon_LP, ~] = ned2geodetic(xa(end), ya(end), 0, settings.lat0, settings.lon0, 0, wgs84Ellipsoid);
% 
% % DATA RECORD (display)
% 
% disp(' ')
% disp('DATA RECORD:')
% fprintf('apogee reached: %g [m] \n', za(end));
% 
% fprintf('time: %g [sec] \n\n', Ta(end))
% 
% fprintf('max speed reached: %g [m/s] \n', max_v)
% fprintf('altitude: %g [m] \n', za(imax_v))
% fprintf('Mach: %g [-] \n', M(imax_v))
% fprintf('time: %g [sec] \n\n', Ta(imax_v))
% 
% fprintf('max Mach reached: %g [-] \n', max_M)
% fprintf('altitude: %g [m] \n', za(imax_M))
% fprintf('velocity: %g [m/s] \n', abs_V(imax_M))
% fprintf('time: %g [sec] \n\n', Ta(imax_M))
% 
% fprintf('max acceleration reached: %g [m/s2] = %g [g] \n', max_a, max_a/9.81)
% fprintf('velocity: %g [m/s] \n', abs_V(imax_a))
% fprintf('time: %g [sec] \n\n', Ta(imax_a))
% 
% fprintf('run on launch pad: %g [m] \n', abs_X(iexit))
% fprintf('speed at launch pad exit: %g [m/s] \n', abs_V(iexit))
% fprintf('time: %g [sec] \n\n', Ta(iexit))
% 
% fprintf('latitude of landing point: %10.8f [deg] \n',lat_LP);
% fprintf('longitude of landing point: %10.8f [deg] \n\n',lon_LP);
% 
% 
% % %% PLOTS
% % if settings.plots
% %     run('plots.m')
% % end
% % 
% % if settings.stoch.N == 1
% %     V_apo = norm(data_ascent.state.Y(end, 4:6) - ...
% %         data_ascent.wind.body_wind(1:3, end)');
% %     fprintf('apogee velocity relative to wind: %g [m/s] \n', V_apo);
% %     delete('ascent_plot.mat')
% % end
% % 
% % clearvars -except T data_ascent data_para data_bal flag LP