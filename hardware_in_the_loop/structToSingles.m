%{

STRUCT TO SINGLES - a nested structs or a matrix composed of singles into a 
                    row array of singles (doesn't send fields named "time")

Author: Emilio Corigliano
Skyward Experimental Rocketry | ELC Dept | electronics@skywarder.eu
email: emilio.corigliano@skywarder.eu
Release date: 10/03/2021

%}

function [data] = structToSingles(inputData)
    if(isfield(inputData, 'time'))
        inputData = rmfield(inputData, 'time');
    end
    
    if(class(inputData)=="struct")
        data = [];
        c = struct2cell(inputData);
        for i= 1:length(c)
            if(class(c{i})=="struct")
                data = [data structToSingles(c{i})]; %#ok<*AGROW>
            else
                data = [data reshape(c{i}',[],1)'];
            end
        end
    else
        data = reshape(inputData',[],1)';
    end
    
end