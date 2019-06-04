classdef Helical
    % antennas.element.Helical   representation of a theoretical helical
    % antenna
    %
    % uses the equation found here:
    % http://www.antenna-theory.com/antennas/travelling/helix.php
    
    % user defined properties that define the antenna
    properties
        Frequency   % the desired operating frequency of the antenna in [GHz]
        N           % the number of turns (larger value means more turns)
        Alpha       % the pitch angle - controls how far the helix grows in [deg]
        MinimumValue    % the minimum value of the pattern in [dB]
        Efficiency  % efficiency of the antenna (as a decimal - e.g 50% is 0.5)
        
        PatterndB   % the pattern in [dB] for all the points on the mesh
        Pattern     % the normalized pattern for all the points on the mesh
    end
    
    % dependent properties based on user set values
    properties (Dependent)
        C   % the circumference of a turn on the helix (a function of the frequency) in [m]
        S   % the vertical separation between turns in [m]
        
        Gain        % the max gain of the antenna in [dB]
        Beamwidth   % the half power beam width in [deg]
    end
    
    % TODO: make these private properties
    properties
        PhiMesh
        ThetaMesh
    end
    
    methods
        function obj = Helical(freq, N, alpha, minValue)
            % Helical   constructor for a theoretical helical antenna
            %   ant = antennas.element.Helical(freq) builds a typical
            %   helical antenna tuned to the desired frequency in GHz.
            %   This antenna will have 10 turns and an alpha value of 13
            %   degrees.
            %
            %   ant = antennas.element.Helical(freq, N) also specifies the
            %   number of turns (which effectively controls the beamwidth).
            %
            %   ant = antennas.element.Helical(freq, N, alpha) also
            %   specifies the alpha value for the antenna.
            %
            %   ant = antennas.element.Helical(freq, N, alpha, minValue)
            %   specifies the minimum value of the dB pattern.
            
            if nargin < 1
                error('not enough inputs, must at least specify frequency');
            end
            
            % set the defaults
            if nargin < 2
                N = 10;
            end
            if nargin < 3
                alpha = 13;
            end
            if nargin < 4
                minValue = -40;
            end
            
            % set the properties
            obj.Frequency = freq;
            obj.N = N;
            obj.Alpha = alpha;
            obj.MinimumValue = minValue;
            obj.Efficiency = 0.7;
            
            % build the mesh for defining the points for the antenna
            [obj.ThetaMesh, obj.PhiMesh] = meshgrid(0:180, 0:360);
            
            % build the entire pattern
            obj.Pattern = obj.response(obj.PhiMesh, obj.ThetaMesh);
            obj.PatterndB = obj.responsedB(obj.PhiMesh, obj.ThetaMesh);
        end
        
        function c = get.C(obj)
            c = 3e8/(obj.Frequency * 1e9);  % C is the wavelength
        end
        
        function s = get.S(obj)
            s = tan(deg2rad(obj.Alpha))*obj.C;
        end
        
        function gdb = get.Gain(obj)
            lambda = 3e8/(obj.Frequency*1e9);   % the wavelength
            g = obj.Efficiency * 12*obj.C^2*obj.N*obj.S/lambda^3;
            gdb = 10*log10(g);
        end
        
        function bw = get.Beamwidth(obj)
            lambda = 3e8/(obj.Frequency*1e9);   % the wavelength
            bw = 52*lambda/(obj.C * sqrt(obj.N*obj.S/lambda));
        end
        
    end
    
    % TODO: maybe these should shift to their own files if this gets too
    % big
    methods
        function Enorm = response(obj, phi, theta)
            % this response will be the normalized E field magnitude
            
            % convert to radian
            thetar = deg2rad(theta);
            
            lambda = 3e8/(obj.Frequency*1e9);   % the wavelength
            k = 2*pi/lambda;                    % the wave constant thing
            Omega = k*obj.S*(cos(thetar) - 1) - pi*(2 + 1/obj.N);
            
            % generate the normalized theta and phi components on the elctric field in
            % a linear scale
            Etheta = sin(pi/(2*obj.N)).*cos(thetar).*sin(obj.N.*Omega./2)./sin(Omega./2);
            Ephi = Etheta;
            
            % need to combine the elements to get the full electric field (?)
            % Though Etheta and Ephi are out of phase from each other....
            % so only really makes sense to combine their magnitudes (which would be
            % their value squared...)
            % should think of Etheta as the real component and Ephi as the imaginary
            % component of a complex number...
            Enorm = sqrt((Etheta.^2 + Ephi.^2)/2);  % dividing by 2 because each of the elements are normalized already
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
            simpant.tools.plot3DGainPattern(obj.ThetaMesh, obj.PhiMesh, obj.PatterndB, obj.MinimumValue);
            
        end
        
        function plotSlice(obj)
            % want to just plot a polar plot of the slice
            slice_phi0 = obj.PatterndB(1,:);
            slice_phi180 = obj.PatterndB(181, end-1:-1:2);
            simpant.tools.plotGainPattern(0:359, [slice_phi0 slice_phi180]);
        end
    end
    
    
    
end