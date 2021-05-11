function plot_all(tot, csett)

if tot.plot_control
%% PLOT THE RESULTS  
tot.ctr_end = 0.1*(length(tot.alpha_degree_tot) - 1) + tot.ctr_start;
time = tot.ctr_start:0.1:tot.ctr_end;
% time_setpoint = 0:0.05:0.05*(length(csett.data_trajectories(csett.chosen_trajectory).V_ref(csett.starting_index:end))-1); % 0.05 discretization step fro trajectory generation

%% Control variable: servo angle
figure('Name','Servo angle after burning phase');
plot(time, tot.alpha_degree_tot);grid on;xlabel('time [s]');ylabel('|alfa| [°]');
title('Servo control action');

if tot.flagPID ~= 3
    % Control variable: pid vs linearization
    figure('Name','Linearization of the control variable');
    plot(time, tot.U_lin_tot, 'DisplayName','Linearized','LineWidth',0.8), grid on;
    hold on
    plot(time, tot.pid_tot, 'DisplayName','PID','LineWidth',0.8), grid on;
    xlabel('time [s]'), ylabel('U [N]');
    hold off
    legend('Location','northeast')

    % delta_S
    figure('Name','Delta_S','NumberTitle','off');
    plot(time, tot.dS_tot), grid on;
    xlabel('time [s]'), ylabel('A [m^2]');

    % Cd
    figure('Name','Cd','NumberTitle','off');
    plot(time, tot.Cd_tot), grid on;
    xlabel('time [s]'), ylabel('Cd []');
end

% Altitude real vs setpoint
figure('Name','Altitude real vs setpoint after burning phase');
plot(time, tot.z_tot,'DisplayName','real','LineWidth',0.8), grid on;
hold on
plot(time, tot.z_setpoint_tot,'DisplayName','setpoint','LineWidth',0.8), grid on;
% plot(time_setpoint, csett.data_trajectories(csett.chosen_trajectory).Z_ref(csett.starting_index:end),'DisplayName','setpoint choosen trajectory','LineWidth',0.8), grid on;
xlabel('time [s]'), ylabel('z [m]');
hold off
legend('Location','southeast')
    
% Vertical velocity real vs setpoint
figure('Name','Vertical velocity real vs setpoint after burning phase');
plot(time, tot.vz_tot,'DisplayName','real','LineWidth',0.8), grid on;
hold on
plot(time, tot.vz_setpoint_tot, 'DisplayName','setpoint', 'LineWidth',0.8), grid on;xlabel('time [s]'), ylabel('Vz [m/s]');
% plot(time_setpoint, csett.data_trajectories(csett.chosen_trajectory).V_ref(csett.starting_index:end),'DisplayName','setpoint choosen trajectory','LineWidth',0.8), grid on;
hold off
legend

% V(z) real vs setpoint
figure('Name','V(z) real vs setpoint after burning phase');
plot(tot.z_tot, tot.vz_tot,'DisplayName','real','LineWidth',0.8), grid on;
hold on
plot(tot.z_setpoint_tot, tot.vz_setpoint_tot, 'DisplayName','setpoint', 'LineWidth',0.8), grid on;
xlabel('z [m]'), ylabel('Vz [m/s]');
hold off
legend
end

if tot.plot_sensors
%% FIGURE: Barometer reads
fbaro = 20;
tp = 0:1/fbaro:1/fbaro*(length(tot.pn_tot)-1);
figure('Name','Barometer reads')
subplot(2,1,1);plot(tp,tot.pn_tot',tot.Tf_tot,tot.p_tot');grid on;xlabel('time [s]');ylabel('|P| [mBar]');
legend('Pressure','Ground-truth','location','southeast');
title('Barometer pressure reads');
subplot(2,1,2);plot(tp,-tot.hb_tot',tot.Tf_tot,-tot.Yf_tot(:,3));grid on;xlabel('time [s]');ylabel('|h| [m]');
legend('Altitude','Ground-truth','location','northeast');
title('Barometer altitude reads');

%% FIGURE: Accelerometer reads
faccel = 100;
ta = 0:1/faccel:1/faccel*(length(tot.accel_tot)-1);
for i = 1:length(ta)
    q(1) = tot.x_est_tot(i,7);
    q(2) = tot.x_est_tot(i,8);
    q(3) = tot.x_est_tot(i,9);
    q(4) = tot.x_est_tot(i,10);
    A              =   [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2,                   2*(q(1)*q(2) + q(3)*q(4)),               2*(q(1)*q(3) - q(2)*q(4));
                            2*(q(1)*q(2) - q(3)*q(4)),          -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,               2*(q(2)*q(3) + q(1)*q(4));
                            2*(q(1)*q(3) + q(2)*q(4)),                   2*(q(2)*q(3) - q(1)*q(4)),     -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2;];

    tot.accel_tot(i,:) = A'*tot.accel_tot(i,:)';
    eul(i,:) = quat2eul([q(4) q(1) q(2) q(3)],'ZYX')*180/pi; 
end

    eulsim = quat2eul(tot.Yf_tot(1:tot.i_apo,10:13),'ZYX')*180/pi; 


figure
subplot(3,1,1);plot(ta',eul(:,1),tot.Tf_tot(1:tot.i_apo),eulsim(:,1));xlabel('time[s]');ylabel('|roll| [°]'); title('Roll');
subplot(3,1,2);plot(ta',eul(:,2),tot.Tf_tot(1:tot.i_apo),eulsim(:,2));xlabel('time[s]');ylabel('|pitch| [°]'); title('Pitch');
subplot(3,1,3);plot(ta',eul(:,3),tot.Tf_tot(1:tot.i_apo),eulsim(:,3));xlabel('time[s]');ylabel('|yaw| [°]'); title('Jaw');


figure('Name','Accelerometer reads')
subplot(3,1,1);plot(ta,tot.accel_tot(:,1)');grid on;xlabel('time[s]');ylabel('|ax| [g]'); title('Accelerometer reads along x');
subplot(3,1,2);plot(ta,tot.accel_tot(:,2)');grid on;xlabel('time[s]');ylabel('|ay| [g]'); title('Accelerometer reads along y');
subplot(3,1,3);plot(ta,tot.accel_tot(:,3)');grid on;xlabel('time[s]');ylabel('|az| [g]'); title('Accelerometer reads along z');

%% FIGURE: Gyroscope reads 
figure('Name','Gyroscope reads')
subplot(3,1,1);plot(ta,tot.gyro_tot(:,1)*180/pi');grid on;xlabel('time[s]');ylabel('|wx| [°/s]'); title('Gyroscope reads along x');
subplot(3,1,2);plot(ta,tot.gyro_tot(:,2)*180/pi');grid on;xlabel('time[s]');ylabel('|wy| [°/s]'); title('Gyroscope reads along y');
subplot(3,1,3);plot(ta,tot.gyro_tot(:,3)*180/pi');grid on;xlabel('time[s]');ylabel('|wz| [°/s]'); title('Gyroscope reads along z'); 

%% FIGURE:Magnetometer reads
figure('Name','Magnetometer reads')
subplot(3,1,1);plot(ta,tot.mag_tot(:,1)');grid on;xlabel('time [s]');ylabel('|mx| [Gauss]'); title('Magnetometer readsalong x'); 
subplot(3,1,2);plot(ta,tot.mag_tot(:,2)');grid on;xlabel('time[s]');ylabel('|my| [Gauss]'); title('Magnetometer reads along y');
subplot(3,1,3);plot(ta,tot.mag_tot(:,3)');grid on;xlabel('time[s]');ylabel('|mz| [Gauss]'); title('Magnetometer reads along z'); 

%% FIGURE: Gps reads 
fgps = 10; 
tgps = 0:1/fgps:1/fgps*(length(tot.gps_tot)-1); 
figure('Name','Gps position reads')
subplot(3,1,1);plot(tgps, tot.gps_tot(:,1)');grid on;xlabel('time [s]');ylabel('|Pn| [m]'); title('GPS position  North'); 
subplot(3,1,2);plot(tgps, tot.gps_tot(:,2)');grid on;xlabel('time [s]');ylabel('|Pe| [m]'); title('GPS position  East');
subplot(3,1,3);plot(tgps,-tot.gps_tot(:,3)');grid on;xlabel('time[s]');ylabel('|Pu| [m]'); title('GPS position Upward'); 
figure('Name','Gps position reads')
subplot(3,1,1);plot(tgps,tot.gpsv_tot(:,1)');grid on;xlabel('time[s]');ylabel('|Velocity N| [m/s]');
subplot(3,1,2);plot(tgps,tot.gpsv_tot(:,2)');grid on;xlabel('time[s]');ylabel('|Velocity E| [m/s]');
subplot(3,1,3);plot(tgps,tot.gpsv_tot(:,3)');grid on;xlabel('time[s]');ylabel('|Velocity D| [m/s]');
title('GPS velocity reads');
end


if tot.plot_kalman
%% FIGURE: Estimated position vs ground-truth
figure('Name','Estimated position vs ground-truth')
subplot(3,1,1);plot(tot.t_est_tot(1:tot.i_apo_est),tot.x_est_tot(1:tot.i_apo_est,1),tot.Tf_tot(1:tot.i_apo), tot.Yf_tot(1:tot.i_apo,1));grid on;xlabel('time[s]');ylabel('|Pn| [m]');legend('North','Ground-truth','location','best');
title('Estimated Northposition vs ground-truth');
subplot(3,1,2);plot(tot.t_est_tot(1:tot.i_apo_est),tot.x_est_tot(1:tot.i_apo_est,2),tot.Tf_tot(1:tot.i_apo), tot.Yf_tot(1:tot.i_apo,2));grid on;xlabel('time[s]');ylabel('|Pe| [m]');legend('East','Ground-truth','location','best'); 
title('Estimated East position vs ground-truth'); 
subplot(3,1,3);plot(tot.t_est_tot(1:tot.i_apo_est),-tot.x_est_tot(1:tot.i_apo_est,3),tot.Tf_tot(1:tot.i_apo), -tot.Yf_tot(1:tot.i_apo,3));grid on;xlabel('time [s]');ylabel('|Pu| [m]');legend('Upward','Ground-truth','location','best'); 
title('Estimated Upward position vs ground-truth'); 

%% FIGURE: Estimated NED velocities vs ground-truth 
figure('Name','Estimated velocities vs ground-truth')
subplot(3,1,1);plot(tot.t_est_tot(1:tot.i_apo_est),tot.x_est_tot(1:tot.i_apo_est,4),tot.Tf_tot(1:tot.i_apo), tot.v_ned_tot(1:tot.i_apo,1));grid on;xlabel('time [s]');ylabel('|Vn| [m/s]');
legend('North','Ground-truth','location','best'); title('Estimated North velocity vs ground-truth'); 
subplot(3,1,2);plot(tot.t_est_tot(1:tot.i_apo_est),tot.x_est_tot(1:tot.i_apo_est,5),tot.Tf_tot(1:tot.i_apo), tot.v_ned_tot(1:tot.i_apo,2));grid on;xlabel('time [s]');ylabel('|Ve| [m/s]');
legend('East','Ground-truth','location','best'); title('Estimated Eastvelocity vs ground-truth');
subplot(3,1,3);plot(tot.t_est_tot(1:tot.i_apo_est+1),-tot.x_est_tot(1:tot.i_apo_est+1,6),tot.Tf_tot(1:tot.i_apo),-tot.v_ned_tot(1:tot.i_apo,3));grid on;xlabel('time [s]');ylabel('|Vu| [m/s]');
legend('Upward','Ground-truth','location','best'); title('EstimatedUpward velocity vs ground-truth'); 
%% FIGURE: Estimated body velocities vs ground-truth 
figure('Name','Estimated body velocities vs ground-truth')
subplot(3,1,1);plot(tot.t_est_tot(1:tot.i_apo_est),tot.vels_tot(1:tot.i_apo_est,1),tot.Tf_tot(1:tot.i_apo), tot.Yf_tot(1:tot.i_apo,4));grid on;xlabel('time [s]');ylabel('|Vn| [m/s]');
legend('x','Ground-truth','location','best'); title('Estimated x velocity vs ground-truth'); 
subplot(3,1,2);plot(tot.t_est_tot(1:tot.i_apo_est),tot.vels_tot(1:tot.i_apo_est,2),tot.Tf_tot(1:tot.i_apo), tot.Yf_tot(1:tot.i_apo,5));grid on;xlabel('time [s]');ylabel('|Ve| [m/s]');
legend('y','Ground-truth','location','best'); title('Estimated y velocity vs ground-truth');
subplot(3,1,3);plot(tot.t_est_tot(1:tot.i_apo_est+1),tot.vels_tot(1:tot.i_apo_est+1,3),tot.Tf_tot(1:tot.i_apo),tot.Yf_tot(1:tot.i_apo,6));grid on;xlabel('time [s]');ylabel('|Vu| [m/s]');
legend('z','Ground-truth','location','best'); title('Estimated z velocity vs ground-truth'); 

%% FIGURE: Estimated quaternions vs ground-truth 
figure('Name','Estimated quaternions vs ground-truth')
subplot(4,1,1);plot(tot.t_est_tot(1:tot.i_apo_est),tot.x_est_tot(1:tot.i_apo_est,10),tot.Tf_tot(1:tot.i_apo),tot.Yf_tot(1:tot.i_apo,10));grid on;ylabel('|q0| [-]'); 
legend('Estimatedq0','Ground-truth','location','northeast'); title('Estimated q0 vsground-truth');
subplot(4,1,2);plot(tot.t_est_tot(1:tot.i_apo_est),tot.x_est_tot(1:tot.i_apo_est,7),tot.Tf_tot(1:tot.i_apo),tot.Yf_tot(1:tot.i_apo,11));grid on;ylabel('|q1| [-]'); 
legend('Estimatedq1','Ground-truth','location','northeast'); title('Estimated q1 vsground-truth');
subplot(4,1,3);plot(tot.t_est_tot(1:tot.i_apo_est),tot.x_est_tot(1:tot.i_apo_est,8),tot.Tf_tot(1:tot.i_apo),tot.Yf_tot(1:tot.i_apo,12));grid on;ylabel('|q2| [-]'); 
legend('Estimatedq2','Ground-truth','location','northeast'); title('Estimated q2 vsground-truth');
subplot(4,1,4);plot(tot.t_est_tot(1:tot.i_apo_est),tot.x_est_tot(1:tot.i_apo_est,9),tot.Tf_tot(1:tot.i_apo),tot.Yf_tot(1:tot.i_apo,13));grid on;ylabel('|q3| [-]'); 
legend('Estimatedq3','Ground-truth','location','northeast'); title('Estimated q3 vsground-truth');



end


if tot.plot_ada 
%% ADA States
figure('Name','ADA States')
subplot(3,1,1);plot(tot.t_ada_tot,tot.xp_ada_tot(:,1));grid on;xlabel('time [s]');ylabel('|P| [mBar]');
title('ADA pressure estimation');
subplot(3,1,2);plot(tot.t_ada_tot,tot.xp_ada_tot(:,2));grid on;xlabel('time [s]');ylabel('|P_dot| [mBar/s]');
title('ADA velocity estimation');
subplot(3,1,3);plot(tot.t_ada_tot,tot.xp_ada_tot(:,3));grid on;xlabel('time [s]');ylabel('|P_dot^2| [mBar/s^2]');
title('ADA acceleration estimation');
%% ADA vertical position and velocity
figure('Name','ADA vertical position and velocity')
subplot(2,1,1);plot(tot.t_ada_tot,tot.xv_ada_tot(:,1));grid on;xlabel('time [s]');ylabel('|P| [mBar]');
title('ADA altitude estimation');
subplot(2,1,2);plot(tot.t_ada_tot,tot.xv_ada_tot(:,2));grid on;xlabel('time [s]');ylabel('|P_dot| [mBar/s]');
title('ADA vertical velocity estimation');
end
end