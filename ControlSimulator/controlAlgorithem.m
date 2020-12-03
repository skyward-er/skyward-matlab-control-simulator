function [alphaOut] = controlAlgorithem(z,Vz,Vu,settings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    persistent alpha alphaOld
    if isempty(alpha)
        alpha=settings.control.alphaDefault; % initialize
        alphaOld=settings.control.alphaDefault;
    end


    %% Control
    if Vu<settings.control.start1state && Vu>settings.control.start2state % start control when velocity under 233 m/s
        %% Controller
        % get index of actual hight z
        if z>=settings.control.hightInterval(end)
            iz = length(settings.control.hightInterval);
        elseif z<=settings.control.hightInterval(1)
            iz = 1;
        else
            iz = interp1(settings.control.hightInterval,1:length(settings.control.hightInterval),z,'nearest');
        end

        % linear interpolate between the values for all_alpha at the actual hight
        % index
        if Vz>=settings.control.all_Vz(iz,end)
            alpha = settings.control.all_alpha(end);
        elseif Vz<=settings.control.all_Vz(iz,1) 
            alpha = settings.control.all_alpha(1);
        else
            alpha = interp1(settings.control.all_Vz(iz,:),settings.control.all_alpha,Vz);
        end
    % do not change alpha at end of the flight (not controllable)
    elseif Vu<=settings.control.start2state
        alpha=alphaOld;
    else
        alpha = settings.control.all_alpha(1);
    end
    
    %% Limit changing speed of alpha
    dAlphaMax = (pi/180*60)/0.13*settings.control.dt; % 60 degree/0.13 sec from servo
    dAlphaMax = dAlphaMax*0.5; % multiplied by a value < 1
    difAlphas=alpha-alphaOld;
    difAlphas=max(difAlphas,-dAlphaMax);
    difAlphas=min(difAlphas,dAlphaMax);
    alpha=alphaOld+difAlphas;
    alphaOld=alpha;
    
    % alpha varies from 0 to 0.89 rad
    alpha=max(alpha,0);
    alpha=min(alpha,0.89);
    
    alphaOut=alpha;
end

