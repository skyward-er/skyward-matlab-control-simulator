classdef Sensor < handle
    %SENSOR Super call for all sensor
    %   Inherit from handle to have state variables which can be changed
    %   inside a method
    
    properties (Access='public')
        minMeasurementRange; % Max limit of sensor
        maxMeasurementRange; % Min limit of sensor
        
        resolution; % resolution of the sensor
        
        noiseVariance; % Varianze for the gaussian white noise
        
        offset; % Offset in all directions
        
        tempOffset; % Coefficent for temperature depending offset
        
        dt; % Sampling time
        
        error2dOffset; % first column: inputArg, second column: relativeArg, third column: error
    end
    
    
    methods (Access='public')
        function obj = Sensor()
            %SENSOR Construct an instance of this class
            %   Sensor is the Superclass of all sensors
            
        end
        
        function outputArg = sens(obj,inputArg,temp)
            %SENS Method to use the sensor.
            %   Gets the simulation data and extract the unideal sensor
            %   output
            %
            %  Inputs:
            %  inputArg sensor data
            %  temp: temperature of the sensor
            %  
            %  Outputs:
            %  outputArg sensor data with nois,
            %  offsets, etc.
            
            inputArg=obj.add2DOffset(inputArg,temp);
            inputArg=obj.whiteNoise(inputArg);
            inputArg=obj.addOffset(inputArg);
            inputArg=obj.addTempOffset(inputArg,temp);
            inputArg=obj.quantization(inputArg);
            inputArg=obj.saturation(inputArg); 
            outputArg = inputArg;
        end
      
    end
    
    
    methods (Access='protected')
        function outputArg = saturation(obj,inputArg)
            %SATURATION Includes saturation to the sensor model
            %   inputArg is limited to minMeasurementRange at the lower end and maxMeasurementRange at the upper end
            %
            %  Necessary properties:
            %  minMeasurementRange: Max limit of sensor
            %  maxMeasurementRange: Min limit of sensor
            %
            %  Inputs:
            %  inputArg: sensor data
            %  
            %  Outputs:
            %  outputArg: sensor data with saturation
            
            % checks if sensor data is higher than max possible value
            if (~isempty(obj.maxMeasurementRange))
                inputArg(inputArg>obj.maxMeasurementRange)=obj.maxMeasurementRange;
            end
            
            % checks if sensor data is lower than min possible value
            if (~isempty(obj.minMeasurementRange))
                inputArg(inputArg<obj.minMeasurementRange)=obj.minMeasurementRange;
            end
            
            outputArg = inputArg;
        end
        
        
        function outputArg = quantization(obj,inputArg)
            %QUATIZATION Quantize the sensor data
            %   Quantizes the signal with the resolution properie of the object
            %
            %  Necessary properties:
            %  resolution: resolution of the sensor
            %
            %  Inputs:
            %  inputArg: sensor data
            %  
            %  Outputs:
            %  outputArg: quantized sensor data
            
            if (~isempty(obj.resolution))
                inputArg = obj.resolution*round(inputArg/obj.resolution);
            end
            outputArg = inputArg;
        end
        
        
        function outputArg = whiteNoise(obj,inputArg)
            %WHITE_NOISE Includes gaussian white noise to the sensor data
            %   Adds gaussian white noise with variance noiseVariance
            %
            %  Necessary properties:
            %  noiseVariance: Varianze for the gaussian white noise
            %
            %  Inputs:
            %  inputArg: sensor data
            %  
            %  Outputs:
            %  outputArg: sensor data with white noise
            
            if (~isempty(obj.noiseVariance))
                inputArg=inputArg+ones(size(inputArg)).*sqrt(obj.noiseVariance).*randn(1,1);
            end
            outputArg = inputArg;
        end
        
        
        function outputArg = addOffset(obj,inputArg)
            %ADD_OFFSET Adds an offset to the signal
            %   adds the propertie offset to the input signal
            %
            %  Necessary properties:
            %  offset: Offset in all directions
            %
            %  Inputs:
            %  inputArg: sensor data
            %  
            %  Outputs:
            %  outputArg: sensor data plus offset
            
            if (~isempty(obj.offset))
                inputArg=inputArg+ones(size(inputArg)).*obj.offset;
            end
            outputArg = inputArg;
        end
        
        
        function outputArg = addTempOffset(obj,inputArg,temp)
            %ADD_Temp_OFFSET Adds an temperature dependen offset to the signal
            %   adds the propertie tempOffset multiplied by the temperature input temp to the input signal
            %
            %  Necessary properties:
            %  tempOffset: Coefficent for temperature depending offset
            %
            %  Inputs:
            %  inputArg: sensor data
            %  temp: temperature of the sensor
            %  
            %  Outputs:
            %  outputArg: sensor data plus temerature depending offset
            
            if (~isempty(obj.tempOffset))
                inputArg=inputArg+ones(size(inputArg)).*temp*obj.tempOffset;
            end
            outputArg = inputArg;
        end
        
        
        function outputArg = add2DOffset(obj,inputArg,temp)
            %ADD_2D_OFFSET Adds a offset that depends on the inputArg and temp
            %  The offset is linear interpolated with the data in matrice error2dOffset
            %
            %  Necessary properties:
            %  error2dOffset: data for offset with first column: inputArg, second column: relativeArg, third column: error
            %
            %  Inputs:
            %  inputArg: sensor data
            %  temp: temperature of the sensor
            %  
            %  Outputs:
            %  outputArg: sensor data plus 2D depending offset
            
            if (~isempty(obj.error2dOffset))
                inputArg=inputArg+...
                    ones(size(inputArg))*griddata(obj.error2dOffset(:,1),obj.error2dOffset(:,2),obj.error2dOffset(:,3),inputArg,temp);
            end
            outputArg = inputArg;
        end
        
    end
end

