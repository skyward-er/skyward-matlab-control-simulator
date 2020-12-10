classdef GPS < Sensor3D
    %GPS GPS Sensor
    %   Uses Sensor3D methods but tronsformes lat and lon degree data in
    %   meter before and afterwards transformes it back to degree
    
    methods
        function obj = GPS()
            %GPS Construct an instance of this class
            %   Calls the constructor of Sensor3D
            
            obj=obj@Sensor3D();
        end
        
        function [outLat,outLon,outAlt] = sens(obj,lat,lon,alt,temp)
            %METHOD1 Summary of this method goes here
            %   Fransforming lat and lat in meter and then call sens method
            %   of Sensor3D https://en.wikipedia.org/wiki/Geographic_coordinate_system
            %
            %  Inputs:
            %  lat: latitude sensor data
            %  lon: latitude sensor data
            %  alt: altidture sensor data
            %  
            %  Outputs:
            %  outLat,outLon,outAlt: sensor data with nois,
            %  offsets, etc.
            
            mLatTrans=111.32e3;
            mLonTrans=(40075e3 * cos(lat)/360);
            
            mLat=mLatTrans*lat;
            mLon=mLonTrans*lon;
            
            [mLat,mLon,alt] = sens@Sensor3D(obj,mLat,mLon,alt,temp);
            
            outLat=mLatTrans^-1*mLat;
            outLon=mLonTrans^-1*mLon;
            outAlt=alt;
        end
    end
end

