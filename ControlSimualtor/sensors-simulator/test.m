close all
clear

% initial sensors
initSensors;

% simulation
t=0:0.01:200;
p=900+100*sin(0.01*2*pi*t);
temp=50;

t=0:0.01:200;
a1=1600+50*sin(0.01*2*pi*t);
a2=1600+50*sin(0.01*2*pi*t-pi*2/3);
a3=1600+50*sin(0.01*2*pi*t-pi*4/3);

for k = 1:length(p)
    ps(k)=MS580301BA01.sens(p(k),temp);
    
    [a1s(k),a2s(k),a3s(k)]=ACCEL_LSM9DS1.sens(a1(k),a2(k),a3(k),temp);
end

figure;
plot(t,p,t,ps);
title('barometer')

figure;
subplot 311
title('accelerometer')
plot(t,a1,t,a1s);
subplot 312
plot(t,a2,t,a2s);
subplot 313
plot(t,a3,t,a3s);

% Calculate magnetic field in magnetometer for sepcific quaternion state
magneticField = 475.937; % magnetic field density in Milan
[mx, my, mz]=createMagneticData(magneticField,cos(pi/4), 0, 0, sin(pi/4));