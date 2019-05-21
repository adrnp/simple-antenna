classdef ArrayAntenna < handle
    % ArrayAntenna is a class to simulate the behavior or an antenna
    % array
    
    properties
        Element             % the antenna model of the element
        ElementPosition     % a vector of the position of each element, as a multiplier of the wavelength
        NumElements
        
        % TODO: maybe set the offsets directly by mode, but for now, let
        % the user do that manually by setting the offsets directly
        ElementPhaseOffsets     % the built in phase offset of each antenna (this helps change the steering mode)
        
        MeasurementNoise    % the standard deviation of the noise on the measurements for a rotation in [dB]
        Sensitivity         % the sensitivity of the sensor reading the incoming gain values (effective min gain) [dB]
    end
    
    
    methods
        
        function obj = ArrayAntenna(element, positions, varargin)
            % ArrayAntenna  constructor for an antenna array with the
            % pattern of each element given by elementAntenna and the
            % positions of the elements within the array given by
            % positions.  Note that positions should be a 3xN matrix with
            % each column containing the [x y z] position of the element in
            % the array as a multiple of the wavelength of the operating
            % frequency of the antenna.
            
            % the first two inputs are required
            if nargin < 2
                error('missing inputs');
            end
            
            % get the number of elements in the array
            [~, nel] = size(positions);
            
            % build the array elements
            obj.Element = element;
            obj.ElementPosition = positions;
            obj.NumElements = nel;
            
            % parse additional inputs
            parser = inputParser;
            parser.addParameter('MeasurementNoise', 0);
            parser.addParameter('Sensitivity', -65);
            parser.addParameter('PhaseOffsets', zeros(1, nel));
            parser.parse(varargin{:});
            res = parser.Results;
            
            % set the additional values
            obj.MeasurementNoise = res.MeasurementNoise;
            obj.Sensitivity = res.Sensitivity;
            obj.ElementPhaseOffsets = res.PhaseOffsets;
            
            % set the element's min value to the sensitivity
            obj.Element.MinimumValue = obj.Sensitivity;
        end
        
        function pattern = fullPattern(obj, phi, theta)
            % fullPattern   compute the resulting gain pattern of the full
            % array steered to the given phi, theta values.  Note that
            % these values default to 0.
            %
            % returns the linear |E| pattern
            
            % handle the defaults
            if nargin < 2
                phi = 0;
            end
            if nargin < 3
                theta = 0;
            end
            
            % compute the array factor
            af = obj.computeArrayFactor(phi, theta);
            
            % TODO: verify that the AF is multiplied by the normalized E
            % field...
            
            % compute the pattern and convert to dB, limiting the values by
            % the limit specified in the element pattern.
            pattern = obj.Element.Pattern .* af;
        end
        
        function patterndB = fullPatterndB(obj, phi, theta)
            % fullPatterndB   compute the resulting gain pattern of the full
            % array steered to the given phi, theta values.  Note that
            % these values default to 0.
            %
            % returns the pattern in [dB] limited by the minimum value
            % specified in the antenna element itself.
            
            % handle the defaults
            if nargin < 2
                phi = 0;
            end
            if nargin < 3
                theta = 0;
            end
            
            % get the linear pattern -> |E|
            pattern = obj.fullPattern(phi, theta);
            
            % convert to linear power (P = |E|^2)
            pattern = pattern.^2;
            
            % convert to dB
            patterndB = 10*log10(pattern);
            patterndB(patterndB < obj.Sensitivity) = obj.Sensitivity;
        end
        
        function setPhaseOffsets(obj, offsets)
            % setPhaseOffsets   set the phase offsets of each of the
            % antennas in the vector to be able to change the operating
            % mode of the antenna (e.g. null or beam steering).
            
            % catch a mistake I might make of passing the data with the
            % wrong dimensions
            if ~isrow(offsets)
                offsets = offsets';
            end
            
            % simply set the offsets
            obj.ElementPhaseOffsets = offsets;
        end
        
        
        function gp = rotate(obj, source, signalStrength, theta)
            % rotate    return the gain pattern that is a result of
            % rotating through all the angles of phi for a given angle of
            % theta.  This is basically a circular slice of data in the
            % horizontal plane.
            %
            % Source is a 2x1 vector of the [phi, theta] angle that points
            % to the source.
            
            % need to know how many sources working with
            [~, Nsources] = size(source);
            
            % the angle through which the pattern will be steered
            phiRotate = 0:1:360;
            thetaSteer = theta;
            
            % need to convert the signal strength to E to be able to add it
            % to each pattern
            % TODO: need to figure out if this is right....
            if nargin < 4
                signalE = ones(Nsources, 1);
            else
                % TODO: need to also account of the effective gain due to
                % the array factor -> this become somewhat in the weeds,
                % can work around it for now by just adjusting the
                % effective signal strength of the incoming signal
                signalE = sqrt(10.^((signalStrength + obj.Element.Gain)/10));
            end
            
            % compute |E| as a result of each of the sources
            Enorm = zeros(Nsources, length(phiRotate));
            for i = 1:Nsources
                sphi = source(1,i);
                stheta = source(2,i);
                
                % get normalized |E|_element in the direction of the ith source
                % (called R(phi, theta) in some of the literature)
                R = obj.Element.response(sphi, stheta);
                
                % include the effect of the signal source strength
%                 R = R * signalE(i);
                R = R;% + signalE(i);
                
                % compute the geometry term of the array factor in the
                % direction of the signal source
                sphir = deg2rad(sphi);
                sthetar = deg2rad(stheta);
                
                dx = obj.ElementPosition(1,:);
                dy = obj.ElementPosition(2,:);
                
                geometryTerms = exp(1i*2*pi*sin(sthetar).*(dx.*cos(sphir) + dy.*sin(sphir)));
                
                % compute |E| for each steered direction
                for j = 1:length(phiRotate)
                    
                    % compute the phase shift needed to apply the desired
                    % weighting to the elements to steer the beam in the
                    % desired direction
                    ps = obj.computePhaseShifts(phiRotate(j), thetaSteer) + obj.ElementPhaseOffsets;
                    psrad = deg2rad(ps);
                    
                    % compute |AF| (normalized array factor magnitude) ->
                    % combines the effects of the geometry of the pattern
                    % and the desired steering weights
                    AF = (1/obj.NumElements) * abs(sum(geometryTerms.*exp(1i*psrad)));
%                     AF = abs(sum(geometryTerms.*exp(1i*psrad)));
                    
                    % |E| for this direction is |R|*|AF|
%                     Enorm(i,j) = R*AF;
                    Enorm(i,j) = R * AF;
                end
            end
            
            % combine all the normalized electric field responses
            if Nsources > 1
                Enorm = sum(Enorm);
            end
            Plinear = Enorm.^2;
            PdB = 10*log10(Plinear) + normrnd(zeros(1, length(Plinear)), obj.MeasurementNoise);
            PdB(PdB < obj.Element.MinimumValue) = obj.Element.MinimumValue;
            
            % create the gain pattern with this normalized data
            gp = GainPattern(thetaSteer, phiRotate, PdB);
        end
        
        function rot = fullRotation(obj, antennaPos, sources, thetaSteer)
            % fullRotation  generate all the data that would result from
            % the full rotation of this antenna.
            %   rot = ant.fullRotation(antennaPos, sources) computes the
            %   rotation data (gain pattern and metadata) for the result of
            %   a full rotation of this antenna from the given position,
            %   antennaPos, as a result of sources, defined as
            %   InterferenceSource objects

            
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
                sp(i) = 10*log10(sources(i).Power*1000) + sources(i).Gain - fspl;
            end
            
            %
            % compute the pattern and the rotation data
            %
            %
            
            % all the parameters above creates just a single situation, so rotate and
            % get the pattern for it
            gp = obj.rotate(sv, sp, thetaSteer);

            % create the rotation object
            rot = GPRotation(gp, length(gp.Phi), ...
                              'TrueBearing', sv(1,:)', ...
                              'TrueTheta', sv(2,:)', ...
                              'Position', antennaPos, ...
                              'HeightAboveSource', antennaPos(3), ...
                              'Distance2DFromSource', d');
        end
        
        function ps = computePhaseShifts(obj, phi, theta)
            % computePhaseShifts    compute the phase shift (angle part of
            % the complex weight value) needed for each of the antenna
            % elements in the array given the desired phi and theta
            % (default to 0) in [deg].
            
            % handle the defaults
            if nargin < 2
                phi = 0;
            end
            if nargin < 3
                theta = 0;
            end
            
            % convert steering angles to rads
            phir = deg2rad(phi);
            thetar = deg2rad(theta);
            
            % get the position elements needed
            dx = obj.ElementPosition(1,:);
            dy = obj.ElementPosition(2,:);
            
            % create a steering vector of the desired length
            ps = rad2deg(-2*pi*sin(thetar).*(dx*cos(phir) + dy*sin(phir)));
        end
        
        
        function plot(obj, phi, theta)
            
            % handle the defaults
            if nargin < 2
                phi = 0;
            end
            if nargin < 3
                theta = 0;
            end
            
            pattern = obj.fullPatterndB(phi, theta);
            plot3DGainPattern(obj.Element.ThetaMesh, obj.Element.PhiMesh, pattern, obj.Element.MinimumValue);
        end
        
        function plotAF(obj, phi, theta)
            % handle the defaults
            if nargin < 2
                phi = 0;
            end
            if nargin < 3
                theta = 0;
            end
            
            af = obj.computeArrayFactor(phi, theta);
            plot3DGainPattern(obj.Element.ThetaMesh, obj.Element.PhiMesh, af, obj.Element.MinimumValue);
        end
        
        function af = getAF(obj, phi, theta)
            
            % handle the defaults
            if nargin < 2
                phi = 0;
            end
            if nargin < 3
                theta = 0;
            end
            
            af = obj.computeArrayFactor(phi, theta);
        end
        
    end
    
    % some private helper methods
    methods (Access = private)
        
        function af = computeArrayFactor(obj, phi, theta)
            % computeArrayFactor    helper function to compute the array
            % factor of the antenna when steered to the given phi and theta
            % values (default to 0) in [deg].
            
            % handle the defaults
            if nargin < 2
                phi = 0;
            end
            if nargin < 3
                theta = 0;
            end
            
            % get the phase shift needed for each of the antennas
            ps = obj.computePhaseShifts(phi, theta) + obj.ElementPhaseOffsets;
            psrad = deg2rad(ps);  % need it as radians
            
            % get the meshes
            pmesh = deg2rad(obj.Element.PhiMesh);
            tmesh = deg2rad(obj.Element.ThetaMesh);
            
            % get the spacing of the elements
            dx = obj.ElementPosition(1,:);
            dy = obj.ElementPosition(2,:);
            
            % compute the array factor, given the phase shift
            af = zeros(size(obj.Element.PhiMesh));
            for i = 1:obj.NumElements
                geometryTerm = exp(1i*2*pi*sin(tmesh).*(dx(i)*cos(pmesh) + dy(i)*sin(pmesh)));
                weightTerm = exp(1i*psrad(i));
                af = af + geometryTerm.*weightTerm;
            end
            
            % normalize for the number of elements
            af = (1/obj.NumElements).*abs(af);
        end
        
        
        function af = computeLimitedArrayFactor(obj, steeringAngle, phi, theta)
            % compute a limited version of the array factor for only a
            % specific slice of data.  In this case steeringAngle is a 2x1
            % array of the desired phi/theta value to steer and phi and
            % theta are the angles for which to compute the array factor.
            %
            % NOTE: this only works with either phi or theta being a single
            % value and the other one being a vector.  This is a function
            % to support the rotate and sweep functions of the pattern.
            
            % get the phase shift needed for each of the antennas
            ps = computePhaseShifts(obj, steeringAngle(1), steeringAngle(2)) + obj.ElementPhaseOffsets;
            psrad = deg2rad(ps);  % need it as radians
            
            % get the meshes
            pmesh = deg2rad(phi);
            tmesh = deg2rad(theta);
            
            % get the spacing of the elements
            dx = obj.ElementPosition(1,:);
            dy = obj.ElementPosition(2,:);
            
            % compute the array factor, given the phase shift
            af = zeros(size(obj.Element.PhiMesh));
            for i = 1:obj.NumElements
                geometryTerm = exp(1i*2*pi*sin(tmesh).*(dx(i)*cos(pmesh) + dy(i)*sin(pmesh)));
                weightTerm = exp(1i*psrad(i));
                af = af + geometryTerm.*weightTerm;
            end
            
            % normalize for the number of elements
            af = (1/obj.NumElements).*af;
            af = abs(af);
        end
        
    end
    
end