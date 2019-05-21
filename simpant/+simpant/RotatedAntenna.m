classdef RotatedAntenna
    % RotatedAntenna    a class to help with the analysis of mounting a
    % single antenna and physically rotating it to get a gain pattern.
    
    properties
        Element     % the antenna that is being rotated
        MountAngle  % the theta angle of the antenna mount in [deg]
        
        MeasurementNoise    % the standard deviation of the noise on the measurements for a rotation in [dB]
        Sensitivity         % the sensitivity of the sensor reading the incoming gain values (effective min gain) [dB]
        
        NumElements = 1
        
        SensorLosses = 5;  % [dB]
    end
    
    methods
        function obj = RotatedAntenna(element, mountAngle, varargin)
            
            % the first two inputs are required
            if nargin < 2
                error('missing inputs');
            end
            
            % set required values
            obj.Element = element;
            obj.MountAngle = mountAngle;
            
            % parse additional inputs
            parser = inputParser;
            parser.addParameter('MeasurementNoise', 0);
            parser.addParameter('Sensitivity', -65);
            parser.addParameter('SensorLosses', 5);
            parser.parse(varargin{:});
            res = parser.Results;
            
            % set the additional values
            obj.MeasurementNoise = res.MeasurementNoise;
            obj.Sensitivity = res.Sensitivity;
            obj.SensorLosses = res.SensorLosses;
            
            % set the element's min value to the sensitivity
            obj.Element.MinimumValue = obj.Sensitivity;
        end
        
        function gp = rotate(obj, source, signalStrength, thetaMount)
            % rotate    generate the gain pattern for a given source
            % direction.  Simulate the result of what a gain pattern would
            % look like from physically rotating the antenna given a source
            % in the source direction.  NOTE: the result is a normalized
            % response which should be scaled appropriately for distance,
            % etc (?).
            %   gp = rotate(source) rotates the angle through a full 360
            %   degrees of azimuth getting one measurement per angle.
            %   Source is a spherical unit vector (as a column vector) in
            %   the direction of the source.  Source can also be a 3xN
            %   matrix to represent N different simultaneous sources.
            %
            %   gp = rotate(source, signalStrength) specifies the signal
            %   strength of the source's signal at the antenna in [dB].
            %   TODO: determine where the antenna gain will be added...
            %
            %   TODO: maybe want to be able to also pass in distances to
            %   the sources or something that can be used to simulate
            %   either different strength or distance... basically to make
            %   some signals weaker than others.
            %
            
            % need to know how many sources working with
            [~, Nsources] = size(source);
            
            if nargin < 4
                thetaMount = obj.MountAngle;
            end
            
            % get the properties of the rotation being done
            phiRotate = 0:360;
            mountRot = roty(-thetaMount);
            
            % need to convert the signal strength to E to be able to add it
            % to each pattern
            % TODO: need to figure out if this is right....
            if nargin < 3
                signalE = ones(Nsources, 1);
            else
                signalE = sqrt(10.^((signalStrength + obj.Element.Gain)/10));
            end
            
            % need to get the Enorm value due to each source
            Enorm = zeros(Nsources, length(phiRotate));
            for i = 1:Nsources
                % handle the first part of the coordinate transformation to put
                % the source phi direction into the rotated frame
                rotatedSource = wrapTo360(source(:,i) - [phiRotate; zeros(1,length(phiRotate)); zeros(1,length(phiRotate))]);

                % for each rotated source, need to convert it to a phi/theta
                % direction in the antenna frame to get the proper strength
                % from the pattern
                antennaFrameSource = xyz2ptr(mountRot*ptr2xyz(rotatedSource));
                antennaFramePhi = antennaFrameSource(1,:);
                antennaFrameTheta = antennaFrameSource(2,:);
                
                % add to the Enorm matrix
                Enorm(i,:) = obj.Element.response(antennaFramePhi, antennaFrameTheta);
                
                % include the E from the source power
                % TODO: not sure if this is the right way to do it...
                Enorm(i,:) = Enorm(i,:) * signalE(i);
            end
            
            % combine all the normalized electric field responses
            if Nsources > 1
                Enorm = sum(Enorm);
            end
            Plinear = Enorm.^2;
            PdB = 10*log10(Plinear) + normrnd(zeros(1, length(Plinear)), obj.MeasurementNoise);
            PdB(PdB < obj.Element.MinimumValue) = obj.Element.MinimumValue;
            
            % create the gain pattern with this normalized data
            gp = GainPattern(mountRot, phiRotate, PdB);
        end
        
        function rot = fullRotation(obj, antennaPos, sources, thetaMount)
            % fullRotation  generate all the data that would result from
            % the full rotation of this antenna.
            %   rot = ant.fullRotation(antennaPos, sources) computes the
            %   rotation data (gain pattern and metadata) for the result of
            %   a full rotation of this antenna from the given position,
            %   antennaPos, as a result of sources, defined as
            %   InterferenceSource objects

            % default to the mount angle, but for simplicity (since it only
            % affects the rotation math) also allow the specification of a
            % different mount angle here.
            if nargin < 4
                thetaMount = obj.MountAngle;
            end
            
            % get the number of sources
            Nsources = length(sources);
            
            %
            % compute source vectors
            %
            [sphi, stheta] = sources.angleFrom(antennaPos);
            [dn, de] = sources.distanceFrom(antennaPos);
            d = sqrt(dn.^2 + de.^2);
            sv = [sphi; stheta; ones(1,Nsources)];
            
            %
            % compute arrival power of the sources
            %
            sp = zeros(1, Nsources);
            for i = 1:Nsources
                fspl = freespacePathLoss(obj.Element.Frequency, sqrt((d(i)/1000).^2 + (antennaPos(3)/1000).^2));
                fspl(fspl == -Inf) = 0;
                sp(i) = 10*log10(sources(i).Power*1000) + sources(i).Gain - fspl - obj.SensorLosses;
            end
            
            %
            % compute the pattern and the rotation data
            %
            %
            
            % all the parameters above creates just a single situation, so rotate and
            % get the pattern for it
            gp = obj.rotate(sv, sp, thetaMount);

            % create the rotation object
            rot = GPRotation(gp, length(gp.Phi), ...
                              'TrueBearing', sv(1,:)', ...
                              'TrueTheta', sv(2,:)', ...
                              'Position', antennaPos, ...
                              'HeightAboveSource', antennaPos(3), ...
                              'Distance2DFromSource', d');
        end
        
        
    end
    
    
    
end