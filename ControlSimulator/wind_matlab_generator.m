function [uw, vw, ww] = wind_matlab_generator(settings, z, t, Hour, Day)
%{

wind_generator - Function that generates wind components in NED reference frame, based on hwm07 model

INPUT:      - settings, structure of rocket data;
            - z, local altitude;
            - t, time sample;
            - Hour, hour of the day of the needed simulation;
            - Day, day of the month of the needed simulation.

OUTPUTS:    - uw, wind component along x;
            - vw, wind component along y;
            - ww, wind component along z.

Author: Gabriele Poiana
Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
email: gabriele.poiana@skywarder.eu
January 2016; Last revision: 17.I.2016

%}

h = -z + settings.z0;
if h < 0
    h = 0;
end

if nargin == 3
    if settings.wind.HourMin == settings.wind.HourMax && settings.wind.HourMin == settings.wind.HourMax
        Day = settings.wind.DayMin;
        Hour = settings.wind.HourMin;
    end
end

Seconds = Hour*3600;

%% HORIZONTAL WIND

[uw,vw] = atmoshwm(settings.lat0,settings.lon0,h,'day',Day,...
    'seconds',Seconds+t,'model','quiet','version','14');    % NED reference
ww = settings.wind.ww;


end


