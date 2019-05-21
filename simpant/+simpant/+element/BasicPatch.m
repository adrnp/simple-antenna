classdef BasicPatch
    %
    %
    % this is based on equations from [1] page 812
    % [1] http://eceweb1.rutgers.edu/~orfanidi/ewa/ch18.pdf
    
    
    properties
        Frequency       % the frequency at which the antenna is operating
        Beamwidth       % the design beamwidth of the antenna in [deg]
        MinimumValue    % the minimum value of the pattern in [dB]
        
        PatterndB       % the pattern in [dB] for all the points on the mesh
        Pattern         % the normalized patter for all the points on the mesh
        
    end
    
    % TODO: make these private properties
    properties
        PhiMesh
        ThetaMesh
    end
    
    properties (Dependent)
        Gain
    end
    
    methods
        function obj = BasicPatch(freq, beamwidth, minValue)
            % handle some of the defaults
            if nargin < 2
                beamwidth = 80;
            end
            
            if nargin < 3
                minValue = -40;
            end
            
            % create the object
            obj.Frequency = freq;
            obj.Beamwidth = beamwidth;
            obj.MinimumValue = minValue;
            
            % build the mesh for defining the points for the antenna
            [obj.ThetaMesh, obj.PhiMesh] = meshgrid(0:180, 0:360);
            
            % build the entire pattern
            obj.Pattern = obj.response(obj.PhiMesh, obj.ThetaMesh);
            obj.PatterndB = obj.responsedB(obj.PhiMesh, obj.ThetaMesh);
        end
        
        function g = get.Gain(obj)
            g = 32383/(obj.Beamwidth*obj.Beamwidth);
        end
        
    end
    
    % TODO: maybe these should shift to their own files if this gets too
    % big
    methods
        function Enorm = response(obj, phi, theta)
            % this response will be the normalized E field magnitude
            
            % TODO: the computation is only a function of theta, but needs
            % to get computed for each of the phi angles desired...
            
            % convert to radian
            phir = deg2rad(phi);
            thetar = deg2rad(theta);
            
            % use the equations from [1] to compute the normalized electric
            % field
            bw = obj.Beamwidth;
            bwMult = 50.76/bw/pi;
            nux = bwMult.*sin(thetar).*cos(phir);
            nuy = bwMult.*sin(thetar).*sin(phir);
            Enorm = (1+cos(thetar))./2.*abs(sinc(pi.*nux).*sinc(pi.*nuy));
        end
        
        function PdB = responsedB(obj, phi, theta)
            % get the normalized E field
            Enorm = obj.response(phi, theta);
            
            % convert to linear power (P = |E|^2)
            Plinear = Enorm.^2;
            
            % convert to db
            PdB = 10*log10(Plinear);
            PdB(PdB < obj.MinimumValue) = obj.MinimumValue;
        end
        
        
        function plot(obj)
            % plot the pattern
            
            % NOTE: need to convert to spherical coordinates
            plot3DGainPattern(obj.ThetaMesh, obj.PhiMesh, obj.PatterndB, obj.MinimumValue);
            
        end
    end
    
end