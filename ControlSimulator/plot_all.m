function plot_all(c, csett)

if c.plot_control
%% PLOT THE RESULTS  
c.ctr_end = 0.1*(length(c.alpha_degree_tot) - 1) + c.ctr_start;
time = c.ctr_start:0.1:c.ctr_end;
% time_setpoint = 0:0.05:0.05*(length(csett.data_trajectories(csett.chosen_trajectory).V_ref(csett.starting_index:end))-1); % 0.05 discretization step fro trajectory generation

%% Control variable: servo angle
figure('Name','Servo angle after burning phase');
plot(time, c.alpha_degree_tot);grid on;xlabel('time [s]');ylabel('|alfa| [°]');
title('Servo control action');

if c.flagPID ~= 3
    % Control variable: pid vs linearization
    figure('Name','Linearization of the control variable');
    plot(time, c.U_lin_tot, 'DisplayName','Linearized','LineWidth',0.8), grid on;
    hold on
    plot(time, c.pid_tot, 'DisplayName','PID','LineWidth',0.8), grid on;
    xlabel('time [s]'), ylabel('U [N]');
    hold off
    legend('Location','northeast')

    % delta_S
    figure('Name','Delta_S','NumberTitle','off');
    plot(time, c.dS_tot), grid on;
    xlabel('time [s]'), ylabel('A [m^2]');

    % Cd
    figure('Name','Cd','NumberTitle','off');
    plot(time, c.Cd_tot), grid on;
    xlabel('time [s]'), ylabel('Cd []');
end

% Altitude real vs setpoint
figure('Name','Altitude real vs setpoint after burning phase');
plot(time, c.z_tot,'DisplayName','real','LineWidth',0.8), grid on;
hold on
plot(time, c.z_setpoint_tot,'DisplayName','setpoint','LineWidth',0.8), grid on;
% plot(time_setpoint, csett.data_trajectories(csett.chosen_trajectory).Z_ref(csett.starting_index:end),'DisplayName','setpoint choosen trajectory','LineWidth',0.8), grid on;
xlabel('time [s]'), ylabel('z [m]');
hold off
legend('Location','southeast')
    
% Vertical velocity real vs setpoint
figure('Name','Vertical velocity real vs setpoint after burning phase');
plot(time, c.vz_tot,'DisplayName','real','LineWidth',0.8), grid on;
hold on
plot(time, c.vz_setpoint_tot, 'DisplayName','setpoint', 'LineWidth',0.8), grid on;xlabel('time [s]'), ylabel('Vz [m/s]');
% plot(time_setpoint, csett.data_trajectories(csett.chosen_trajectory).V_ref(csett.starting_index:end),'DisplayName','setpoint choosen trajectory','LineWidth',0.8), grid on;
hold off
legend

% V(z) real vs setpoint
figure('Name','V(z) real vs setpoint after burning phase');
plot(c.z_tot, c.vz_tot,'DisplayName','real','LineWidth',0.8), grid on;
hold on
plot(c.z_setpoint_tot, c.vz_setpoint_tot, 'DisplayName','setpoint', 'LineWidth',0.8), grid on;
xlabel('z [m]'), ylabel('Vz [m/s]');
hold off
legend
end

if c.plot_sensors
%% FIGURE: Barometer reads
fbaro = 20;
tp = 0:1/fbaro:1/fbaro*(length(c.pn_tot)-1);
figure('Name','Barometer reads')
subplot(2,1,1);plot(tp,c.pn_tot',c.Tf_tot,c.p_tot');grid on;xlabel('time [s]');ylabel('|P| [mBar]');
legend('Pressure','Ground-truth','location','southeast');
title('Barometer pressure reads');
subplot(2,1,2);plot(tp,-c.hb_tot',c.Tf_tot,-c.Yf_tot(:,3));grid on;xlabel('time [s]');ylabel('|h| [m]');
legend('Altitude','Ground-truth','location','northeast');
title('Barometer altitude reads');

%% FIGURE: Accelerometer reads
faccel = 100;
ta = 0:1/faccel:1/faccel*(length(c.accel_tot)-1);
for i = 1:length(ta)
    q(1) = c.x_est_tot(i,7);
    q(2) = c.x_est_tot(i,8);
    q(3) = c.x_est_tot(i,9);
    q(4) = c.x_est_tot(i,10);
    A              =   [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2,                   2*(q(1)*q(2) + q(3)*q(4)),               2*(q(1)*q(3) - q(2)*q(4));
                            2*(q(1)*q(2) - q(3)*q(4)),          -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,               2*(q(2)*q(3) + q(1)*q(4));
                            2*(q(1)*q(3) + q(2)*q(4)),                   2*(q(2)*q(3) - q(1)*q(4)),     -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2;];

    c.accel_tot(i,:) = A'*c.accel_tot(i,:)';
    eul(i,:) = quat2eul([q(4) q(1) q(2) q(3)],'ZYX')*180/pi; 
end

    eulsim = quat2eul(c.Yf_tot(1:c.i_apo,10:13),'ZYX')*180/pi; 


figure
subplot(3,1,1);plot(ta',eul(:,1),c.Tf_tot(1:c.i_apo),eulsim(:,1));xlabel('time[s]');ylabel('|roll| [°]'); title('Roll');
subplot(3,1,2);plot(ta',eul(:,2),c.Tf_tot(1:c.i_apo),eulsim(:,2));xlabel('time[s]');ylabel('|pitch| [°]'); title('Pitch');
subplot(3,1,3);plot(ta',eul(:,3),c.Tf_tot(1:c.i_apo),eulsim(:,3));xlabel('time[s]');ylabel('|yaw| [°]'); title('Jaw');


figure('Name','Accelerometer reads')
subplot(3,1,1);plot(ta,c.accel_tot(:,1)');grid on;xlabel('time[s]');ylabel('|ax| [g]'); title('Accelerometer reads along x');
subplot(3,1,2);plot(ta,c.accel_tot(:,2)');grid on;xlabel('time[s]');ylabel('|ay| [g]'); title('Accelerometer reads along y');
subplot(3,1,3);plot(ta,c.accel_tot(:,3)');grid on;xlabel('time[s]');ylabel('|az| [g]'); title('Accelerometer reads along z');

%% FIGURE: Gyroscope reads 
figure('Name','Gyroscope reads')
subplot(3,1,1);plot(ta,c.gyro_tot(:,1)*180/pi');grid on;xlabel('time[s]');ylabel('|wx| [°/s]'); title('Gyroscope reads along x');
subplot(3,1,2);plot(ta,c.gyro_tot(:,2)*180/pi');grid on;xlabel('time[s]');ylabel('|wy| [°/s]'); title('Gyroscope reads along y');
subplot(3,1,3);plot(ta,c.gyro_tot(:,3)*180/pi');grid on;xlabel('time[s]');ylabel('|wz| [°/s]'); title('Gyroscope reads along z'); 

%% FIGURE:Magnetometer reads
figure('Name','Magnetometer reads')
subplot(3,1,1);plot(ta,c.mag_tot(:,1)');grid on;xlabel('time [s]');ylabel('|mx| [Gauss]'); title('Magnetometer readsalong x'); 
subplot(3,1,2);plot(ta,c.mag_tot(:,2)');grid on;xlabel('time[s]');ylabel('|my| [Gauss]'); title('Magnetometer reads along y');
subplot(3,1,3);plot(ta,c.mag_tot(:,3)');grid on;xlabel('time[s]');ylabel('|mz| [Gauss]'); title('Magnetometer reads along z'); 

%% FIGURE: Gps reads 
fgps = 10; 
tgps = 0:1/fgps:1/fgps*(length(c.gps_tot)-1); 
figure('Name','Gps position reads')
subplot(3,1,1);plot(tgps, c.gps_tot(:,1)');grid on;xlabel('time [s]');ylabel('|Pn| [m]'); title('GPS position  North'); 
subplot(3,1,2);plot(tgps, c.gps_tot(:,2)');grid on;xlabel('time [s]');ylabel('|Pe| [m]'); title('GPS position  East');
subplot(3,1,3);plot(tgps,-c.gps_tot(:,3)');grid on;xlabel('time[s]');ylabel('|Pu| [m]'); title('GPS position Upward'); 
figure('Name','Gps position reads')
subplot(3,1,1);plot(tgps,c.gpsv_tot(:,1)');grid on;xlabel('time[s]');ylabel('|Velocity N| [m/s]');
subplot(3,1,2);plot(tgps,c.gpsv_tot(:,2)');grid on;xlabel('time[s]');ylabel('|Velocity E| [m/s]');
subplot(3,1,3);plot(tgps,c.gpsv_tot(:,3)');grid on;xlabel('time[s]');ylabel('|Velocity D| [m/s]');
title('GPS velocity reads');
end


if c.plot_kalman
%% FIGURE: Estimated position vs ground-truth
figure('Name','Estimated position vs ground-truth')
subplot(3,1,1);plot(c.t_est_tot(1:c.i_apo_est),c.x_est_tot(1:c.i_apo_est,1),c.Tf_tot(1:c.i_apo), c.Yf_tot(1:c.i_apo,1));grid on;xlabel('time[s]');ylabel('|Pn| [m]');legend('North','Ground-truth','location','best');
title('Estimated Northposition vs ground-truth');
subplot(3,1,2);plot(c.t_est_tot(1:c.i_apo_est),c.x_est_tot(1:c.i_apo_est,2),c.Tf_tot(1:c.i_apo), c.Yf_tot(1:c.i_apo,2));grid on;xlabel('time[s]');ylabel('|Pe| [m]');legend('East','Ground-truth','location','best'); 
title('Estimated East position vs ground-truth'); 
subplot(3,1,3);plot(c.t_est_tot(1:c.i_apo_est),-c.x_est_tot(1:c.i_apo_est,3),c.Tf_tot(1:c.i_apo), -c.Yf_tot(1:c.i_apo,3));grid on;xlabel('time [s]');ylabel('|Pu| [m]');legend('Upward','Ground-truth','location','best'); 
title('Estimated Upward position vs ground-truth'); 

%% FIGURE: Estimated NED velocities vs ground-truth 
figure('Name','Estimated velocities vs ground-truth')
subplot(3,1,1);plot(c.t_est_tot(1:c.i_apo_est),c.x_est_tot(1:c.i_apo_est,4),c.Tf_tot(1:c.i_apo), c.v_ned_tot(1:c.i_apo,1));grid on;xlabel('time [s]');ylabel('|Vn| [m/s]');
legend('North','Ground-truth','location','best'); title('Estimated North velocity vs ground-truth'); 
subplot(3,1,2);plot(c.t_est_tot(1:c.i_apo_est),c.x_est_tot(1:c.i_apo_est,5),c.Tf_tot(1:c.i_apo), c.v_ned_tot(1:c.i_apo,2));grid on;xlabel('time [s]');ylabel('|Ve| [m/s]');
legend('East','Ground-truth','location','best'); title('Estimated Eastvelocity vs ground-truth');
subplot(3,1,3);plot(c.t_est_tot(1:c.i_apo_est+1),-c.x_est_tot(1:c.i_apo_est+1,6),c.Tf_tot(1:c.i_apo),-c.v_ned_tot(1:c.i_apo,3));grid on;xlabel('time [s]');ylabel('|Vu| [m/s]');
legend('Upward','Ground-truth','location','best'); title('EstimatedUpward velocity vs ground-truth'); 
%% FIGURE: Estimated body velocities vs ground-truth 
figure('Name','Estimated body velocities vs ground-truth')
subplot(3,1,1);plot(c.t_est_tot(1:c.i_apo_est),c.vels_tot(1:c.i_apo_est,1),c.Tf_tot(1:c.i_apo), c.Yf_tot(1:c.i_apo,4));grid on;xlabel('time [s]');ylabel('|Vn| [m/s]');
legend('x','Ground-truth','location','best'); title('Estimated x velocity vs ground-truth'); 
subplot(3,1,2);plot(c.t_est_tot(1:c.i_apo_est),c.vels_tot(1:c.i_apo_est,2),c.Tf_tot(1:c.i_apo), c.Yf_tot(1:c.i_apo,5));grid on;xlabel('time [s]');ylabel('|Ve| [m/s]');
legend('y','Ground-truth','location','best'); title('Estimated y velocity vs ground-truth');
subplot(3,1,3);plot(c.t_est_tot(1:c.i_apo_est+1),c.vels_tot(1:c.i_apo_est+1,3),c.Tf_tot(1:c.i_apo),c.Yf_tot(1:c.i_apo,6));grid on;xlabel('time [s]');ylabel('|Vu| [m/s]');
legend('z','Ground-truth','location','best'); title('Estimated z velocity vs ground-truth'); 

%% FIGURE: Estimated quaternions vs ground-truth 
figure('Name','Estimated quaternions vs ground-truth')
subplot(4,1,1);plot(c.t_est_tot(1:c.i_apo_est),c.x_est_tot(1:c.i_apo_est,10),c.Tf_tot(1:c.i_apo),c.Yf_tot(1:c.i_apo,10));grid on;ylabel('|q0| [-]'); 
legend('Estimatedq0','Ground-truth','location','northeast'); title('Estimated q0 vsground-truth');
subplot(4,1,2);plot(c.t_est_tot(1:c.i_apo_est),c.x_est_tot(1:c.i_apo_est,7),c.Tf_tot(1:c.i_apo),c.Yf_tot(1:c.i_apo,11));grid on;ylabel('|q1| [-]'); 
legend('Estimatedq1','Ground-truth','location','northeast'); title('Estimated q1 vsground-truth');
subplot(4,1,3);plot(c.t_est_tot(1:c.i_apo_est),c.x_est_tot(1:c.i_apo_est,8),c.Tf_tot(1:c.i_apo),c.Yf_tot(1:c.i_apo,12));grid on;ylabel('|q2| [-]'); 
legend('Estimatedq2','Ground-truth','location','northeast'); title('Estimated q2 vsground-truth');
subplot(4,1,4);plot(c.t_est_tot(1:c.i_apo_est),c.x_est_tot(1:c.i_apo_est,9),c.Tf_tot(1:c.i_apo),c.Yf_tot(1:c.i_apo,13));grid on;ylabel('|q3| [-]'); 
legend('Estimatedq3','Ground-truth','location','northeast'); title('Estimated q3 vsground-truth');



end


if c.plot_ada 
%% ADA States
figure('Name','ADA States')
subplot(3,1,1);plot(c.t_ada_tot,c.xp_ada_tot(:,1));grid on;xlabel('time [s]');ylabel('|P| [mBar]');
title('ADA pressure estimation');
subplot(3,1,2);plot(c.t_ada_tot,c.xp_ada_tot(:,2));grid on;xlabel('time [s]');ylabel('|P_dot| [mBar/s]');
title('ADA velocity estimation');
subplot(3,1,3);plot(c.t_ada_tot,c.xp_ada_tot(:,3));grid on;xlabel('time [s]');ylabel('|P_dot^2| [mBar/s^2]');
title('ADA acceleration estimation');
%% ADA vertical position and velocity
figure('Name','ADA vertical position and velocity')
subplot(2,1,1);plot(c.t_ada_tot,c.xv_ada_tot(:,1));grid on;xlabel('time [s]');ylabel('|P| [mBar]');
title('ADA altitude estimation');
subplot(2,1,2);plot(c.t_ada_tot,c.xv_ada_tot(:,2));grid on;xlabel('time [s]');ylabel('|P_dot| [mBar/s]');
title('ADA vertical velocity estimation');
end
end