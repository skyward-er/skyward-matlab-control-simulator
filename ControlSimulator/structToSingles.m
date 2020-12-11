function [data] = structToSingles(sensorData)
% converts a struct, a nested structs or a matrix composed of singles into
% a row array of singles
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