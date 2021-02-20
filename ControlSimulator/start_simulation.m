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

if settings.electronics
    [Yf, Tf, cpuTimes, flagMatr] = std_run_control(settings);
else
    [Yf, Tf, cpuTimes, flagMatr, data_flight] = std_run_control(settings);
end

%% DATA-PRINTING

Na = length(Yf(:,1));

% POSITIONS
xa = Yf(:,1);
ya = Yf(:,2);
za = -Yf(:,3);
Xa = [xa, ya, za];

[max_z, i_apo] = max(za);
T_apo = Tf(i_apo);

% VELOCITIES
ua = Yf(:,4);
va = Yf(:,5);
wa = -Yf(:,6);
Va = [ua, va, wa];

% MAXIMUM POSITIONS, VELOCITIES AND ACCELERATION
abs_X = vecnorm(Xa');
abs_V = vecnorm(Va');

[max_dist, imax_dist] = max(abs_X);
[max_v, imax_v] = max(abs_V);

% DATA RECORD (display)
fprintf('OUTCOMES:\n\n\n')

fprintf('total computational Time: %.3f [s]: \n', sum(cpuTimes))
fprintf('mean step computational Time: %.3f [s]: \n', mean(cpuTimes))
fprintf('max step omputational Time: %.3f [s]: \n\n', max(cpuTimes))

fprintf('apogee: %.1f [m] \n', max_z);
fprintf('@time: %g [sec] \n\n', T_apo)

fprintf('max speed reached: %g [m/s] \n', max_v)
fprintf('@altitude: %g [m] \n', za(imax_v))
fprintf('@time: %g [sec] \n\n', Tf(imax_v))

if not(settings.electronics)
    M = data_flight.interp.M;
    [max_M, imax_M] = max(M);
    A = data_flight.accelerations.body_acc;
    abs_A = vecnorm(A');
    [max_a, imax_a] = max(abs_A);
    
    fprintf('max Mach reached: %g [-] \n', max_M)
    fprintf('@altitude: %g [m] \n', za(imax_M))
    fprintf('@velocity: %g [m/s] \n', abs_V(imax_M))
    fprintf('@time: %g [sec] \n\n', Tf(imax_M))
    
    fprintf('max acceleration reached: %g [m/s2] = %g [g] \n', max_a, max_a/9.81)
    fprintf('@velocity: %g [m/s] \n', abs_V(imax_a))
    fprintf('@time: %g [sec] \n\n', Tf(imax_a))

end

%% PLOT 

if settings.plots && not(settings.electronics)
    
    % AERO FORCES
    figure('Name', 'Forces - ascent Phase', 'NumberTitle', 'off');
    plot(Tf(flagMatr(:, 2)), data_flight.forces.AeroDyn_Forces(:, 1)),grid on;
    xlabel('Time [s]'); ylabel('X-body force [N]')
    
    % CA
    figure('Name', 'Axial Drag Coefficient - ascent Phase', 'NumberTitle', 'off');
    plot(Tf(flagMatr(:, 2)), data_flight.coeff.CA), title('Axial Drag Coefficient vs time'), grid on;
    xlabel('Time [s]'); ylabel('CA [/]')
    
    % ACCELERATION
    figure('Name', 'Velocity-Abs, Acceleration-Abs - ascent Phase', 'NumberTitle', 'off');
    subplot(1, 2, 1)
    plot(Tf(flagMatr(:, 2)), abs_V(flagMatr(:, 2))), grid on;
    xlabel('time [s]'), ylabel('|V| [m/s]');
    
    subplot(1, 2, 2)
    plot(Tf(flagMatr(:, 2)), abs_A/9.81), grid on;
    xlabel('time [s]'), ylabel('|A| [g]');
       
end

clearvars -except Yf data_flight settings
