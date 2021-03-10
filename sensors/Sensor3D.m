classdef Sensor3D < Sensor
    
    % Author: Jan Hammelman
    % Skyward Experimental Rocketry | ELC-SCS Dept | electronics@skywarder.eu
    % email: jan.hammelmann@skywarder.eu,alessandro.delduca@skywarder.eu
    % Release date: 01/03/2021
    
    %SENSOR_3D For all 3D sensors
    %   Extends Sensor with methods, which work with 3D sensors

    properties (Access='public')
        offsetX; % offset in x direction
        offsetY; % offset in y direction
        offsetZ; % offset in z direction
        
        walkDiffusionCoef; % diffusion coefficient for the random walk
        
        transMatrix % transformation matrix
    end
    
    properties (Dependent=false)
        stateWalkX; % State of random walk in x direction
        stateWalkY; % State of random walk in y direction
        stateWalkZ; % State of random walk in z direction
    end
    
    methods (Access='public')
        function obj = Sensor3D()
            %SENSOR_3D Construct an instance of this class
            %   Sets the random walk states to zero
            
            obj.stateWalkX=0;
            obj.stateWalkY=0;
            obj.stateWalkZ=0;
        end
        
        function [outputArgX,outputArgY,outputArgZ] = sens(obj,inputArgX,inputArgY,inputArgZ,temp)
            %SENS Method to use the sensor.
            %   Gets the simulation data and extract the unideal sensor
            %   output
            %
            %  Inputs:
            %  inputArgX,inputArgY,inputArgZ: sensor data
            %  temp: temperature of the sensor
            %  
            %  Outputs:
            %  outputArgX,outputArgY, outputArgZ: sensor data with nois,
            %  offsets, etc.
            
            [inputArgX,inputArgY,inputArgZ] = tranformAxis(obj,inputArgX,inputArgY,inputArgZ);
            
            inputArgX=obj.add2DOffset(inputArgX,temp);
            inputArgY=obj.add2DOffset(inputArgY,temp);
            inputArgZ=obj.add2DOffset(inputArgZ,temp);
            
            [inputArgX,inputArgY,inputArgZ] = addOffset3D(obj,inputArgX,inputArgY,inputArgZ);
            
            inputArgX=obj.whiteNoise(inputArgX);
            inputArgY=obj.whiteNoise(inputArgY);
            inputArgZ=obj.whiteNoise(inputArgZ);
            
            inputArgX=obj.addOffset(inputArgX);
            inputArgY=obj.addOffset(inputArgY);
            inputArgZ=obj.addOffset(inputArgZ);
            
            inputArgX=obj.addTempOffset(inputArgX,temp);
            inputArgY=obj.addTempOffset(inputArgY,temp);
            inputArgZ=obj.addTempOffset(inputArgZ,temp);
            
            [inputArgX,inputArgY,inputArgZ] = obj.randomWalk(inputArgX,inputArgY,inputArgZ);
            
            inputArgX=obj.quantization(inputArgX);
            inputArgY=obj.quantization(inputArgY);
            inputArgZ=obj.quantization(inputArgZ);
            
            inputArgX=obj.saturation(inputArgX);
            inputArgY=obj.saturation(inputArgY);
            inputArgZ=obj.saturation(inputArgZ);
            
            outputArgX = inputArgX;
            outputArgY = inputArgY;
            outputArgZ = inputArgZ;
        end
    end
    
    methods (Access='protected')
        function [outputArgX,outputArgY,outputArgZ] = addOffset3D(obj,inputArgX,inputArgY,inputArgZ)
            %ADD_OFFSET Adds an 3D offset to the signal
            %   adds the properties offsetX, offsetY and offsetZ to the
            %   input signal components
            %
            %  Necessary properties:
            %  offsetX, offsetY, offsetZ: Offset for all directions
            %
            %  Inputs:
            %  inputArgX,inputArgY,inputArgZ: sensor data
            %  
            %  Outputs:
            %  outputArgX,outputArgY, outputArgZ: sensor data with offset
            
            if (~isempty(obj.offsetX))
                inputArgX=inputArgX+ones(size(inputArgX)).*obj.offsetX;
            end
            if (~isempty(obj.offsetY))
                inputArgY=inputArgY+ones(size(inputArgY)).*obj.offsetY;
            end
            if (~isempty(obj.offsetZ))
                inputArgZ=inputArgZ+ones(size(inputArgZ)).*obj.offsetZ;
            end
            outputArgX=inputArgX;
            outputArgY=inputArgY;
            outputArgZ=inputArgZ;
        end
        
        function [outputArgX,outputArgY, outputArgZ] = randomWalk(obj,inputArgX,inputArgY,inputArgZ)
            %RANDOM_WALK Adds a random walk
            %  Adds a brown motion with the diffusion coefficient walkDiffusionCoef for each step
            %
            %  Necessary properties:
            %  walkDiffusionCoef: diffusion coefficient for the random walk
            %  dt: Sampling time
            %
            %  Inputs:
            %  inputArgX,inputArgY,inputArgZ: sensor data
            %  
            %  Outputs:
            %  outputArgX,outputArgY, outputArgZ: sensor data with random
            %  walk

            if (~isempty(obj.walkDiffusionCoef))
                % distance
                s = sqrt ( 2.0 * 3 * obj.walkDiffusionCoef * obj.dt ) * randn ( 1, 1);

                % direction
                x=randn(1,1);
                y=randn(1,1);
                z=randn(1,1);

                % normalize
                v=s/sqrt(x^2+y^2+z^2);
                dx=x*v;
                dy=y*v;
                dz=z*v;

                % sum random walk
                obj.stateWalkX=obj.stateWalkX+dx;
                obj.stateWalkY=obj.stateWalkY+dy;
                obj.stateWalkZ=obj.stateWalkZ+dz;
                
            end
            
            outputArgX=inputArgX+obj.stateWalkX;
            outputArgY=inputArgY+obj.stateWalkY;
            outputArgZ=inputArgZ+obj.stateWalkZ;
        end
        
        function [outputArgX,outputArgY, outputArgZ] = tranformAxis(obj,inputArgX,inputArgY,inputArgZ)
            %RANDOM_WALK Adds a random walk
            %  Adds a brown motion with the diffusion coefficient walkDiffusionCoef for each step
            %
            %  Necessary properties:
            %  transMatrix: transformation matrix
            %
            %  Inputs:
            %  inputArgX,inputArgY,inputArgZ: sensor data
            %  
            %  Outputs:
            %  outputArgX,outputArgY, outputArgZ: sensor data with
            %  transformation xout=transMatrix*xin
            
            if (~isempty(obj.transMatrix))
                transOut=obj.transMatrix*[inputArgX;inputArgY;inputArgZ];
                inputArgX=transOut(1);
                inputArgY=transOut(2);
                inputArgZ=transOut(3);
            end
            
            outputArgX=inputArgX;
            outputArgY=inputArgY;
            outputArgZ=inputArgZ;
        end
            
    end
end

