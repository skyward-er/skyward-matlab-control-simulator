%% KALMAN FILTER WITH DIRECT KINEMATIC MODEL
%---------Generation of the data-----------
close all;clearvars;clc

addpath('Rocket_simulation\Rocket_simplified_dynamics')

run('param.m')                                      %Generation of necessary parameters

simout  =   sim('rocket_state_sim.slx');            %Running of simulation

%---------Allocation of variables:-------------------

t           =   simout.tout;                        %Time instants of the simulation

v           =   simout.v;                           %Matrix of velocity of the rocket
                                                    %(each row is a vector for each time instant)
acc         =   simout.acc-[0,0,9.81];                         %Matrix of acceleration of the rocket
                                                    %(each row is a vector for each time instant)
w           =   simout.w;                           %Matrix of angular velocity of the rocket
                                                    %(each row is a vector for each time instant)
q           =   simout.q;                           %Matrix of quaternions of the rocket
                                                    %(each row is a vector for each time instant)
h           =   simout.h;                           %Vector of altitude of the rocket
                                                    %(each row corresponds to a time instant)
x           =   simout.x;                           %Vector of x position of the rocket
                                                    %(each row corresponds to a time instant)
y           =   simout.y;                           %Vector of y position of the rocket
                                                    %(each row corresponds to a time instant)

A           =   zeros(3,3,length(t));               %Pre-allocation of rotation matrix
acc_body    =   zeros(length(t),3);                 %Pre-allocation of acceleration in body axes
v_mod       =   zeros(length(t),1);                 %Pre-allocation of velocity modulus
mag_body    =   zeros(length(t),3);                 %Pre-allocation of magnetic field in body axes

for i=1:length(t)
    q_mat           =  [0        -q(i,3)     q(i,2);         %Matrix needed for the
                        q(i,3)    0         -q(i,1);         %definition of rotation
                       -q(i,2)  q(i,1)      0     ;];        %matrix
                   
    A(:,:,i)        =   (q(i,4)^2-q(i,1:3)*q(i,1:3)')*eye(3) + 2*q(i,1:3)'*q(i,1:3)...
                        - 2*q(i,4)*q_mat;                   %Rotation matrix to
                                                            %body axis from inertial
    acc_body(i,:)   =   (A(:,:,i)*acc(i,:)')';
    
    v_mod(i)        =   (norm(v(i,:)))^2;                   %Obtention of the modulus of the velocity  

    mag_body(i,:)   =   (A(:,:,i)*m_inertial')';            %magnetic field in the body axis
end



%% ---------Transformation into measurements-------------

%These up are the ''real'' values with which I will compare the integration
%Now they have to be sampled

t_sam       =   0:1/f_sample1:t_sim;                 %Sampling instants for the IMU

w_sam       =   interp1(t,w,t_sam);                  %Sampled angular velocities

rng('default')                                       %Definition of the seed 
                                                     %of randn function this
                                                     %is for repeatability of simulation

noise(:,1)  =   sigma_w*randn(length(w_sam(:,1)),1); %Generation of the noise 
                                                     %for component 1 of
                                                     %the angular velocity
noise(:,2)  =   sigma_w*randn(length(w_sam(:,1)),1); %Generation of the noise 
                                                     %for component 1 of
                                                     %the angular velocity
noise(:,3)  =   sigma_w*randn(length(w_sam(:,1)),1); %Generation of the noise 
                                                     %for component 1 of
                                                     %the angular velocity
                                                     
%The three components of the noise are assumed to be uncorrelated

w_sam       =   noise + w_sam;                         %Addition of the noise

% figure
% title('Comparison of real angular velocity and noisy measurements')
% plot(t,w,'DisplayName','Real w')
% hold on
% plot(t_sam,w_sam,'DisplayName','Noisy measurement of w')
% xlabel('t [s]')
% ylabel('w [rad/s]')
% grid on
% legend('Location','best')


acc_sam     =    interp1(t,acc_body,t_sam);                 %Sampled linear acceleration

noise2(:,1) =   sigma_a*randn(length(acc_sam(:,1)),1); %Generation of the noise 
                                                       %for component 1 of
                                                       %the linear
                                                       %acceleration
noise2(:,2) =   sigma_a*randn(length(acc_sam(:,1)),1); %Generation of the noise 
                                                       %for component 2 of
                                                       %the linear
                                                       %acceleration
noise2(:,3) =   sigma_a*randn(length(acc_sam(:,1)),1); %Generation of the noise 
                                                       %for component 3 of
                                                       %the linear
                                                       %acceleration
                                                     
%The three components of the noise are assumed to be uncorrelated

acc_sam     =   acc_sam + noise2;                         %Addition of noise


% figure
% title('Comparison of real linear acceleration and noisy measurements')
% plot(t,acc_body,'DisplayName','Real a')
% hold on
% plot(t_sam,acc_sam,'DisplayName','Noisy measurement of a')
% xlabel('t [s]')
% ylabel('a [m/s^2]')
% grid on
% legend('Location','best')


%With this, the measurements that are going to be used as the inputs are
%generated. Now it's necesasry to generate the measurements with which
%correct the prediction

t_sam_mag      =   0:1/f_sample_mag:t_sim;               %Sampling instants 
                                                         %for the
                                                         %magnetometer
                                                         
                                                         
% q_sam       =   interp1(t,q,t_sam2);                     %Sampled quaternions
% 
% 
% noise3(:,1) =   sigma_q*randn(length(q_sam(:,1)),1);   %Generation of the noise 
%                                                        %for component 1 of
%                                                        %the quaternion
% noise3(:,2) =   sigma_q*randn(length(q_sam(:,1)),1);   %Generation of the noise 
%                                                        %for component 2 of
%                                                        %the quaternion
% noise3(:,3) =   sigma_q*randn(length(q_sam(:,1)),1);   %Generation of the noise 
%                                                        %for component 3 of
%                                                        %the quaternion
% noise3(:,4) =   sigma_q*randn(length(q_sam(:,1)),1);   %Generation of the noise 
%                                                        %for component 4 of
%                                                        %the quaternion
%                                                        
% q_sam       =   q_sam + noise3;                        %Addition of noise


% figure
% title('Comparison of real quaternion and noisy measurements')
% plot(t,q,'DisplayName','Real q')
% hold on
% plot(t_sam2,q_sam,'DisplayName','Noisy measurement of q')
% xlabel('t [s]')
% ylabel('q [-]')
% grid on
% legend('Location','best')


mag_sam       =   interp1(t,mag_body,t_sam_mag);            %Sampled magentic field in body frame


noise3(:,1) =   sigma_mag*randn(length(mag_sam(:,1)),1);   %Generation of the noise 
                                                           %for component 1 of
                                                           %the quaternion
noise3(:,2) =   sigma_mag*randn(length(mag_sam(:,1)),1);   %Generation of the noise 
                                                           %for component 2 of
                                                           %the quaternion
noise3(:,3) =   sigma_mag*randn(length(mag_sam(:,1)),1);   %Generation of the noise 
                                                           %for component 3 of
                                                           %the quaternion
                                                       
mag_sam       =   mag_sam + noise3;                        %Addition of noise

for i=1:length(t_sam_mag)
    mag_sam(i,:) = mag_sam(i,:)/norm(mag_sam(i,:));        %Normalisation of the 
                                                           %magnetic field direction after adding noise
end

% figure
% title('Comparison of real magentic field in body axis and noisy measurements')
% plot(t,mag_body,'DisplayName','Real magnetic vector')
% hold on
% plot(t_sam_mag,mag_sam,'DisplayName','Noisy measurement of magnetic field')
% xlabel('t [s]')
% ylabel('m [H]')
% grid on
% legend('Location','best')

t_sam_h      =   0:1/f_sample_h:t_sim;               %Sampling instants 
                                                         %for the
                                                         %barometer

h_sam       =   interp1(t,h,t_sam_h)';                  %Sampling of altitude  

noise4      =   sigma_h*randn(length(h_sam(:,1)),1);     %Generation of the noise 
                                                       %for altitude
h_sam       =   h_sam + noise4;                          %Addition of noise

% figure
% title('Comparison of real altitude and noisy measurements')
% plot(t,h,'DisplayName','Real h')
% hold on
% plot(t_sam_h,h_sam,'DisplayName','Noisy measurement of h')
% xlabel('t [s]')
% ylabel('h [m]')
% grid on
% legend('Location','best')

t_sam_GPS      =   0:1/f_sample_GPS:t_sim;               %Sampling instants 
                                                         %for the GPS

x_sam       =   interp1(t,x,t_sam_GPS)';                 %Position x sampling

noise5      =   sigma_GPS*randn(length(x_sam(:,1)),1);   %Generation of the noise
                                                         %for the position x
x_sam       =   x_sam + noise5;                          %Addition of noise

% figure
% title('Comparison of real position x and noisy measurements')
% plot(t,x,'DisplayName','Real x')
% hold on
% plot(t_sam_GPS,x_sam,'DisplayName','Noisy measurement of x')
% xlabel('t [s]')
% ylabel('x [m]')
% grid on
% legend('Location','best')

y_sam       =   interp1(t,y,t_sam_GPS)';                %Position y sampling

noise6      =   sigma_GPS*randn(length(y_sam(:,1)),1);  %Generation of the noise
                                                        %for the position y
y_sam       =   y_sam + noise6;                         %Addition of the noise

% figure
% title('Comparison of real position y and noisy measurements')
% plot(t,y,'DisplayName','Real y')
% hold on
% plot(t_sam_GPS,y_sam,'DisplayName','Noisy measurement of y')
% xlabel('t [s]')
% ylabel('x [m]')
% grid on
% legend('Location','best')

h_sam2       =   interp1(t,h,t_sam_GPS)';                   %Position z sampling from GPS

noise8      =   sigma_GPS*randn(length(h_sam2(:,1)),1);     %Generation of the noise
                                                            %for the
                                                            %position z for
                                                            %the GPS
h_sam2       =   h_sam2 + noise8;                           %Addition of noise


t_sam_pitot      =   0:1/f_sample_pitot:t_sim;               %Sampling instants 
                                                            %for the
                                                            %Pitot tube

v_mod_sam    =   interp1(t,v_mod,t_sam_pitot)';                   %Velocity modulus sampling

noise7       =   sigma_pitot*randn(length(v_mod_sam(:,1)),1);  %Generation of the noise
                                                        %for the velocity
                                                        %modulus
v_mod_sam       =   v_mod_sam + noise7;                         %Addition of the noise

% figure
% title('Comparison of real velocity modulus and noisy measurements')
% plot(t,v_mod,'DisplayName','Real v^2')
% hold on
% plot(t_sam_pitot,v_mod_sam,'DisplayName','Noisy measurement of v^2')
% xlabel('t [s]')
% ylabel('x [m]')
% grid on
% legend('Location','best')

%% Start of the state estimation

t_eval          =   0:dt:t_sim;

P               =   zeros(10,10,length(t_eval));
x_es            =   zeros(length(t_eval),10);
ab              =   zeros(length(t_eval),3);
a               =   zeros(length(t_eval),3);
w               =   zeros(length(t_eval),3);


x_es(1,:)       =   [r0;v0;q0]';                        %Initial condition 
                                                        %of the estimation 
                                                        %to the real initial condition
                                        
P(:,:,1)        =   Q0;                                 %Initial condition 
                                                        %of the covariance
                                                        %matrix P = Q
                                                        
y_res_GPS       =   zeros(3,length(t_sam_GPS));
y_res_pitot     =   zeros(1,length(t_sam_pitot));
y_res_bar       =   zeros(1,length(t_sam_h));
y_res_mag       =   zeros(3,length(t_sam_mag));
index_GPS=2;
index_pitot=2;
index_bar=2;
index_mag=2;

for i=2:length(t_eval)
    ab(i-1,:)                   =   interp1(t_sam,acc_sam,t_eval(i-1),'previous');  %Evaluation of the 
                                                    %acceleration at estimation 
                                                    %instant with
                                                    %piece-wise constant measurements
     w (i-1,:)                  =   interp1(t_sam,w_sam,t_eval(i-1),'previous');  %Evaluation of the 
                                                %angular velocity at estimation 
                                                %instant with
                                                %piece-wise constant measurements
                                               

    [x_es(i,:),P(:,:,i),a(i,:)] = kalmanFilterPrediction(t_eval(i),x_es(i-1,:),dt,P(:,:,i-1),ab(i-1,:),w(i-1,:),Q0);
    
    if t_eval(i)>=t_sam_GPS(index_GPS)  %Comparison to see the there's a new measurement
       [x_es(i,:),P(:,:,i),y_res_GPS(:,index_GPS)]     = correctionGPS(x_es(i,:),P(:,:,i),x_sam(index_GPS),...
                            y_sam(index_GPS),h_sam2(index_GPS),sigma_GPS,4,1);
        index_GPS   =  index_GPS + 1;
    end
    
    if t_eval(i)>=t_sam_pitot(index_pitot) %Comparison to see the there's a new measurement
       [x_es(i,:),P(:,:,i),y_res_pitot(index_pitot)] = correctionPitot(x_es(i,:),P(:,:,i),v_mod_sam(index_pitot),sigma_pitot);
       index_pitot  =  index_pitot + 1;     
    end
    
    if t_eval(i)>=t_sam_h(index_bar) %Comparison to see the there's a new measurement
       [x_es(i,:),P(:,:,i),y_res_bar(index_bar)]     = correctionBarometer(x_es(i,:),P(:,:,i),h_sam(index_bar),sigma_h);
        index_bar   =  index_bar + 1;     
    end
     
    if t_eval(i)>=t_sam_mag(index_mag) %Comparison to see the there's a new measurement
       [x_es(i,:),P(:,:,i),y_res_mag(:,index_mag)]     = correctionMagnetometer(x_es(i,:),P(:,:,i),mag_sam(index_mag,:),sigma_mag);
       index_mag    =  index_mag + 1;  
    end
    
end

%% Plots


figure
plot(t_eval,x_es(:,4),'LineWidth',2,'DisplayName','\fontsize{16}Estimated vx')
hold on
plot(t_eval,x_es(:,5),'LineWidth',2,'DisplayName','\fontsize{16}Estimated vy')
plot(t_eval,-x_es(:,6),'LineWidth',2,'DisplayName','\fontsize{16}Estimated vz')
plot(t,v(:,1),t,v(:,2),t,-v(:,3),'LineWidth',2,'DisplayName','\fontsize{16}Real v')
title('\fontsize{16}Comparison of real velocities and estimated ones')
xlabel('\fontsize{16}t [s]')
ylabel('\fontsize{16}v [m/s]')
grid on
legend('Location','best')

% figure
% plot(t_sam2,z,'DisplayName','Estimated modulus')
% hold on
% plot(t_sam2,v_mod_sam,'DisplayName','Real modulus')
% title('Comparison of real velocity modulus and estimated one')
% xlabel('t [s]')
% ylabel('v^2 [m^2/s^2]')
% grid on
% legend('Location','best')



% figure
% plot(x,h,'DisplayName','Real trajectory')
% hold on
% plot(x_es(:,1),-x_es(:,3),'DisplayName','Estimated trajectory')
% title('Comparison of real trajectory and estimated one')
% xlabel('x [m]')
% ylabel('h [m]')
% grid on
% legend('Location','best')


figure
plot(t,-h,'LineWidth',2,'DisplayName','\fontsize{16}Real altitude')
hold on
plot(t_eval,-x_es(:,3),'LineWidth',2,'DisplayName','\fontsize{16}Estimated h')
title('\fontsize{16}Comparison of real altitude and estimated one')
xlabel('\fontsize{16}t [s]')
ylabel('\fontsize{16}h [m]')
grid on
legend('Location','best')

figure
plot(t,x,'LineWidth',2,'DisplayName','\fontsize{16}Real x')
hold on
plot(t_eval,x_es(:,1),'LineWidth',2,'DisplayName','\fontsize{16}Estimated x')
title('\fontsize{16}Comparison of real x position and estimated one')
xlabel('\fontsize{16}t [s]')
ylabel('\fontsize{16}x [m]')
grid on
legend('Location','best')


figure
plot(t_eval,x_es(:,7),'LineWidth',2,'DisplayName','\fontsize{16}Estimated q1')
hold on
plot(t_eval,x_es(:,8),'LineWidth',2,'DisplayName','\fontsize{16}Estimated q2')
plot(t_eval,x_es(:,9),'LineWidth',2,'DisplayName','\fontsize{16}Estimated q3')
plot(t_eval,x_es(:,10),'LineWidth',2,'DisplayName','\fontsize{16}Estimated q4')
plot(t,q,'LineWidth',2,'DisplayName','\fontsize{16}Real q')
title('\fontsize{16}Comparison of real quaternions and estimated ones')
xlabel('\fontsize{16}t [s]')
ylabel('\fontsize{16}q [-]')
grid on
legend('Location','best')
grid on
legend('Location','best')


figure
plot(t,acc(:,1),t,acc(:,2),t,-acc(:,3),'LineWidth',2,'DisplayName','\fontsize{16}real acceleration (without gravity)')
hold on
plot(t_eval,a(:,1),t_eval,a(:,2),t_eval,-a(:,3),'LineWidth',2,'DisplayName','\fontsize{16}Rotated acceleration (without gravity)')
title('\fontsize{16}Comparison of real acceleration and the one fed to the estimator')
grid on
legend('Location','best')


figure
plot(t_sam_GPS,y_res_GPS(1,:),'DisplayName','GPS residual x axis')
hold on
plot(t_sam_GPS,y_res_GPS(2,:),'DisplayName','GPS residual y axis')
plot(t_sam_GPS,y_res_GPS(3,:),'DisplayName','GPS residual z axis')
title('Plot of the residual of the position')
xlabel('t [s]')
ylabel('p [m]')
grid on
legend('Location','best')

figure
plot(t_sam_mag,y_res_mag,'DisplayName','Magnetometes residual')
title('Plot of the residual of the magnetometers')
xlabel('t [s]')
ylabel('angle [°]')
grid on
legend('Location','best')

figure
plot(t_sam_h,y_res_bar,'DisplayName','Altimeter residual')
title('Plot of the residual of the altitude')
xlabel('t [s]')
ylabel('h [m]')
grid on
legend('Location','best')

figure
plot(t_sam_pitot,y_res_pitot,'DisplayName','V^2 residual')
title('Plot of the residual of the Velocity')
xlabel('t [s]')
ylabel('V2norm [m^2/s^2]')
grid on
legend('Location','best')

residual_mag      =     mean(y_res_mag,2);
residual_Pitot    =     mean(y_res_pitot,2);
residual_Altitude =     mean(y_res_bar,2);
residual_GPS      =     mean(y_res_GPS,2);

disp("Sample mean of the Magnetometers residual")
disp(residual_mag)
disp("Sample mean of the Altitude residual")
disp(residual_Altitude)
disp("Sample mean of the V^2 residual")
disp(residual_Pitot)
disp("Sample mean of the GPS residual")
disp(residual_GPS)






