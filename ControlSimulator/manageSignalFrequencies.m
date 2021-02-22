function [sensorData] = manageSignalFrequencies(magneticFieldApprox, flagAscent, settings, Y, T, x, uw, vw, ww, uncert)

%{

manageSignalFrequencies - This function interpolates the rocket State to
                          output the measures at different frequencies

INTPUTS:
            - freq, acquisition frequencies of the sensors;
            - Y, State Matrix containing the rocket states specified in     
                 ascent.m @ time step given in T
            - T time step vector

OUTPUTS:
            - sensorData, struct containing the interpolated sensor data

Author: Adriano Filippo Inno
Skyward Experimental Rocketry | AFD Dept
email: adriano.filippo.inno@skywarder.eu
Release date: 30/11/2020

%}

% acc = accelerationOnlyAscent(t, Y, settings, c, uw, vw, ww, uncert)

freq = settings.frequencies;

%% accelerometer
if freq.accelerometerFrequency > freq.controlFrequency
    N = freq.accelerometerFrequency/freq.controlFrequency;
    sensorData.accelerometer.time = linspace(T(1), T(end), N);
    if settings.ballisticFligth || (not(settings.ballisticFligth) && flagAscent)
        for i = 1:N
            iTimeAcc = sensorData.accelerometer.time(i);
            if all(iTimeAcc ~= T)
                [index0] = find(iTimeAcc < T);
                index1 = index0(1);
                index0 = index1 - 1;
                Y1 = Y(index1, :);
                Y0 = Y(index0, :);
                T1 = T(index1);
                T0 = T(index0);
                % linear interpolation between the 2 states
                m = (Y1 - Y0)./(T1 - T0);
                q = Y1 - m*T1;
                Yinterp = m*iTimeAcc + q;
            else
                Yinterp = Y(iTimeAcc == T, :);
                
            end
            sensorData.accelerometer.measures(i, :) = accelerationOnlyAscent...
                (iTimeAcc, Yinterp, settings, x, uw, vw, ww, uncert);
        end
    else
        sensorData.accelerometer.measures(1:N, 1:3) = repmat(zeros(1, 3), N, 1);
    end
else
    sensorData.accelerometer.measures(1, :) = accelerationOnlyAscent...
            (T(end), Y(end, :), settings, x, uw, vw, ww, uncert);
    sensorData.accelerometer.time = T(end);
end

%% gyro
if freq.gyroFrequency > freq.controlFrequency
    N = freq.gyroFrequency/freq.controlFrequency;
    if N ~= round(N)
        error('the sensor frequency must be a multiple of the control frequency');
    end
    sensorData.gyro.time = linspace(T(1), T(end), N);
    if settings.ballisticFligth || (not(settings.ballisticFligth) && flagAscent)
        for i = 1:N
            iTimeGyro = sensorData.gyro.time(i);
            if all(iTimeGyro ~= T)
                [index0] = find(iTimeGyro < T);
                index1 = index0(1);
                index0 = index1 - 1;
                Y1 = Y(index1, 7:9);
                Y0 = Y(index0, 7:9);
                T1 = T(index1);
                T0 = T(index0);
                % linear interpolation between the 2 states
                m = (Y1 - Y0)./(T1 - T0);
                q = Y1 - m*T1;
                sensorData.gyro.measures(i, :) = m*iTimeGyro + q;
            else
                sensorData.gyro.measures(i, :) = Y(iTimeGyro == T, 7:9);
            end
        end
    else
        sensorData.gyro.measures(1:N, 1:3) = repmat(zeros(1, 3), N, 1);
    end
else
    sensorData.gyro.measures(1, :) = Y(end, 7:9);
    sensorData.gyro.time = T(end);
end

%% magnetometer
if freq.magnetometerFrequency > freq.controlFrequency
    N = freq.magnetometerFrequency/freq.controlFrequency;
    if N ~= round(N)
        error('the sensor frequency must be a multiple of the control frequency');
    end
    sensorData.magnetometer.time = linspace(T(1), T(end), N);
    Q = zeros(N, 4);
    z = zeros(1, N);
    if settings.ballisticFligth || (not(settings.ballisticFligth) && flagAscent)
        for i = 1:N
            iTimeMagnetometer = sensorData.magnetometer.time(i);
            if all(iTimeMagnetometer ~= T)
                [index0] = find(iTimeMagnetometer < T);
                index1 = index0(1);
                index0 = index1 - 1;
                Y1 = Y(index1, [3, 10:13]);
                Y0 = Y(index0, [3, 10:13]);
                T1 = T(index1);
                T0 = T(index0);
                % linear interpolation between the 2 states
                m = (Y1 - Y0)./(T1 - T0);
                q = Y1 - m*T1;
                Yinterp = m*iTimeMagnetometer + q;
                z(i) = - Yinterp(1);
                Q(i, :) = Yinterp(2:end);
            else
                z(i) = - Y(iTimeMagnetometer == T, 3);
                Q(i, :) = Y(iTimeMagnetometer == T, 10:13);
            end
        end
        magnFieldInertial = magneticFieldApprox(z + settings.z0)';
        sensorData.magnetometer.measure = quatrotate(Q, magnFieldInertial);
    else
        sensorData.magnetometer.measure = repmat(magneticFieldApprox(0)', N, 1);
        a = sensorData.magnetometer.measure
    end
else
    z = -Y(end, 3);
    Q = Y(end, 10:13);
    magnFieldInertial = magneticFieldApprox(z + settings.z0)';
    sensorData.magnetometer.measure = quatrotate(Q, magnFieldInertial);
    sensorData.magnetometer.time = T(end);
end

%% gps
if freq.gpsFrequency > freq.controlFrequency
    N = freq.gpsFrequency/freq.controlFrequency;
    if N ~= round(N)
        error('the sensor frequency must be a multiple of the control frequency');
    end
    sensorData.gps.time = linspace(T(1), T(end), N);
    if settings.ballisticFligth || (not(settings.ballisticFligth) && flagAscent)
        for i = 1:N
            iTimegps = sensorData.gps.time(i);
            if all(iTimegps ~= T)
                [index0] = find(iTimegps < T);
                index1 = index0(1);
                index0 = index1 - 1;
                Y1 = Y(index1, [1:6, 10:13]);
                Y0 = Y(index0, [1:6, 10:13]);
                T1 = T(index1);
                T0 = T(index0);
                % linear interpolation between the 2 states
                m = (Y1 - Y0)./(T1 - T0);
                q = Y1 - m*T1;
                Yinterp = m*iTimegps + q;
                sensorData.gps.positionMeasures(i, :) = Yinterp(1:3);
                sensorData.gps.velocityMeasures(i, :) = quatrotate(quatconj(Yinterp(7:10)), Yinterp(4:6));
            else
                Q = Y(iTimegps == T, 10:13);
                V = Y(iTimegps == T, 4:6);
                sensorData.gps.positionMeasures(i, :) = Y(iTimegps == T, 1:3);
                sensorData.gps.velocityMeasures(i, :) = quatrotate(quatconj(Q), V);
            end
        end
    else
        
        for i = 1:N
            iTimegps = sensorData.gps.time(i);
            if all(iTimegps ~= T)
                [index0] = find(iTimegps < T);
                index1 = index0(1);
                index0 = index1 - 1;
                Y1 = Y(index1, 1:6);
                Y0 = Y(index0, 1:6);
                T1 = T(index1);
                T0 = T(index0);
                % linear interpolation between the 2 states
                m = (Y1 - Y0)./(T1 - T0);
                q = Y1 - m*T1;
                Yinterp = m*iTimegps + q;
                sensorData.gps.positionMeasures(i, :) = Yinterp(1:3);
                sensorData.gps.velocityMeasures(i, :) = Yinterp(4:6);
            else
                sensorData.gps.positionMeasures(i, :) = Y(iTimegps == T, 1:3);
                sensorData.gps.velocityMeasures(i, :) = Y(iTimegps == T, 4:6);
            end
        end

    end
else
    if settings.ballisticFligth || (not(settings.ballisticFligth) && flagAscent)
        Q = Y(end, 10:13);
        V = Y(end, 4:6);
        sensorData.gps.positionMeasures(1, :) = Y(end, 1:3);
        sensorData.gps.velocityMeasures(1, :) = quatrotate(quatconj(Q), V);
        sensorData.gps.time = T(end);
    else
        sensorData.gps.positionMeasures(1, :) = Y(end, 1:3);
        sensorData.gps.velocityMeasures(1, :) = Y(end, 4:6);
        sensorData.gps.time = T(end);
    end
end

sensorData.gps.positionMeasures(:, 3) = -sensorData.gps.positionMeasures(:, 3);
sensorData.gps.velocityMeasures(:, 3) = -sensorData.gps.velocityMeasures(:, 3);



%% barometer 
if freq.barometerFrequency > freq.controlFrequency
    N = freq.barometerFrequency/freq.controlFrequency;
    z = zeros(N, 1);
    if N ~= round(N)
        error('the sensor frequency must be a multiple of the control frequency');
    end
    sensorData.barometer.time = linspace(T(1), T(end), N);
    for i = 1:N
        iTimeBarometer = sensorData.barometer.time(i);
        if all(iTimeBarometer ~= T)
            [index0] = find(iTimeBarometer < T);
            index1 = index0(1);
            index0 = index1 - 1;
            Y1 = -Y(index1, 3);
            Y0 = -Y(index0, 3);
            T1 = T(index1);
            T0 = T(index0);
            % linear interpolation between the 2 states
            m = (Y1 - Y0)./(T1 - T0);
            q = Y1 - m*T1;
            z(i) = m*iTimeBarometer + q;

        else
            z(i) = -Y(iTimeBarometer == T, 3);
        end
    end
    
else
    sensorData.barometer.time = T(end);
    z = -Y(end, 3);
end

[Temp, ~, P, ~] = atmosisa(z);
sensorData.barometer.measures = P;
sensorData.barometer.temperature = Temp;
