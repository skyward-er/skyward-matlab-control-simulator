%{

STRUCT TO SINGLES - a nested structs or a matrix composed of singles into a 
                    row array of singles (doesn't send fields named "time")

Author: Emilio Corigliano
Skyward Experimental Rocketry | ELC Dept | electronics@skywarder.eu
email: emilio.corigliano@skywarder.eu
Release date: 10/03/2021

%}

function [data] = structToSingles(sensorData)
    if(isfield(sensorData, 'time'))
        sensorData = rmfield(sensorData, 'time');
    end
    
    if(class(sensorData)=="struct")
        data = [];
        c = struct2cell(sensorData);
        for i= 1:length(c)
            if(class(c{i})=="struct")
                data = [data structToSingles(c{i})]; %#ok<*AGROW>
            elseif(class(c{i})=="double")
                data = [data reshape(c{i}',[],1)'];
            end
        end
    elseif(class(sensorData)=="double")
        data = reshape(sensorData',[],1)';
    end
end