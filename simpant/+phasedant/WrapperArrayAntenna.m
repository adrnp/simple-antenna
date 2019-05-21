classdef WrapperArrayAntenna < handle
% WrapperArrayAntenna   a wrappre around the phased array toolbox antenna
% array to be able to interact with the phased array toolbox in the way I
% want to be able to (e.g. with phi/theta as the angles, etc).

    properties
        Element
        PhasedElement
        ElementPosition
        NumElements
        
        % TODO: maybe set the offsets directly by mode, but for now, let
        % the user do that manually by setting the offsets directly
        ElementPhaseOffsets     % the built in phase offset of each antenna (this helps change the steering mode)
        
        MeasurementNoise    % the standard deviation of the noise on the measurements for a rotation in [dB]
        Sensitivity         % the sensitivity of the sensor reading the incoming gain values (effective min gain) [dB]

        SensorLosses = 5            % general losses in the sensor
        AntennaEfficiency = 0.8     % percentage for efficiency of the antenna
    end
    
    % properties that contain Phased Array Toolbox elements
    properties
        Array
        SteerVec
        ArrayResponse
        
        PtotalW  % total power in [W] -> needed for directivity calculation
    end
    
    
    methods
        function obj = WrapperArrayAntenna(element, positions, varargin)
            % TODO: figure out the inputs, etc
            
            % get the number of elements in the array
            [~, nel] = size(positions);
            
            % build the array elements
            obj.Element = element;
            obj.PhasedElement = simpant.element.convertToMatlabElement(element);
            obj.ElementPosition = positions;
            obj.NumElements = nel;
            
            % parse additional inputs
            parser = inputParser;
            parser.addParameter('MeasurementNoise', 0);
            parser.addParameter('Sensitivity', -67);
            parser.addParameter('PhaseOffsets', zeros(1, nel));
            parser.addParameter('SensorLosses', 5);
            parser.parse(varargin{:});
            res = parser.Results;
            
            % set the additional values
            obj.MeasurementNoise = res.MeasurementNoise;
            obj.Sensitivity = res.Sensitivity;
            obj.ElementPhaseOffsets = res.PhaseOffsets;
            obj.SensorLosses = res.SensorLosses;
            
            % set the element's min value to the sensitivity
            obj.Element.MinimumValue = obj.Sensitivity;
            
            % build the phased array conformal array
            % compute the wavelength for the spacing of the antennas
            c = physconst('LightSpeed');            % speed of light in [m/s]
            lambda = c/(element.Frequency*1e9);     % the wavelength of the tuned frequency in [m]
            
            % create the phased array toolbox elements
            obj.Array = phased.ConformalArray('ElementPosition', obj.ElementPosition.*lambda, ...
                              'Element', obj.PhasedElement);
            
            % need helper to generate the steering vectors
            obj.SteerVec = phased.SteeringVector('SensorArray', obj.Array);

            % need helper to get the response from the antenna
            obj.ArrayResponse = phased.ArrayResponse('SensorArray', obj.Array, 'WeightsInputPort', true);
            
            % precompute the total power in W
            obj.PtotalW = obj.computePtotal();
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
        
        
        function Prad = computePtotal(obj)
            % compute the Ptotal for the antenna -> this is in Watts
            fhz = obj.Element.Frequency*1e9;
            myArrayResp = phased.ArrayResponse('SensorArray',obj.Array);
            az = -180:179;
            el = -90:89;
            Prad = 0;
            for m = 1:numel(el)
                Prad = Prad+sum(abs(step(myArrayResp,fhz,...
                [az;el(m)*ones(1,numel(az))])).^2*cosd(el(m)));
            end
            Prad = Prad*2*pi^2/numel(az)/numel(el);
        end
        
        
    end
    
    
    % computation functions
    methods
        function pattern = fullPattern(obj, phi, theta)
            % TODO: implement a function to compute the resulting gain
            % pattern of the full array steered to the given phi/theta
            % values
            %
            % should return the linear |E| pattern
        end
        
        %%% USING DIRECTIVITY AND EFFICIENCY %%%
        function gp = rotate(obj, source, signalStrength, theta)
            % TODO: return the gain pattern that is a result of rotating
            % through all the angles of phi for a given angle of theta.
            % this is basically a circular cslice of data in the horizontal
            % plane.
            %
            % NOTE: source is a 2x1 vector of the [phi, theta] angle that
            % points to the source
            %
            % NOTE: this should be able to handle multiple sources
            
            % source information (convert to AzEl)
            [~, Nsources] = size(source);
            sourceAzEl = phitheta2azel([source(1,:); source(2,:)]);
            
            % signal strength is the received signal strength in [dBm]
            % going to compute the directivity and account for the
            % efficiency to get the antenna gain in [dBi]
            % the resulting RSS is signal strength + antenna gain
            
            % steering information
            phiRotate = 0:1:360;
            thetaSteer = theta;
            nphi = length(phiRotate);
            
            % compute the steering vector for all the steering angles
            % desired
            fhz = obj.Element.Frequency*1e9;
            svs = obj.SteerVec(fhz, phitheta2azel([phiRotate; thetaSteer*ones(1, nphi)]));
            
            
            % compute the directivity in the direction of the signal source
            % for all the different steering vectors
%             d = directivity(obj.Array, fhz*ones(1, nphi), sourceAzEl, 'Weights', svs);

            % self computing the directivity using
            % https://www.mathworks.com/help/phased/ug/element-and-array-radiation-patterns-and-responses.html
            
            Enorm = abs(obj.ArrayResponse(fhz*ones(1,nphi), sourceAzEl, svs));  % the response
            
            % TODO: not sure if this should be power or mag
            % I'm thinknig it sohuld be mag since I haven't squared it
            % yet...
            signalMag = db2mag(signalStrength);
            Enorm = Enorm .* repmat(signalMag', 1, nphi);
            
            % combine the responses at this level
            % TODO: need to account for the signal strength of the
            % source...
            if Nsources > 1
                Enorm = sum(Enorm);
            end
            
            d_W = 4*pi/obj.PtotalW * Enorm.^2;
            
            % account for efficiency
            d_W = obj.AntennaEfficiency * d_W;
            
            % account for the signal strength of the incoming signal
            % NOTE: this needs to be done in the power relm
%             signalPower = db2pow(signalStrength - 30);  % NOTE: signal strength is in dBm and this fn takes dB to convert to W
%             d_W = d_W + repmat(signalPower', 1, nphi);
            
            % convert to dBi
            d_dbi = pow2db(d_W);

            % account for any sensor losses, and adjust for sensitivity
            rss = d_dbi - obj.SensorLosses;
            rss(rss < obj.Sensitivity) = obj.Sensitivity;

            % create the gain pattern with this normalized data
            gp = GainPattern(thetaSteer, phiRotate, rss);
        end
        
        
        function gp = getNominalPattern(obj)
            % get the nominal pattern for this array -> this is the
            % response to an unsteered pattern, gain is irrelevant, it's
            % more about the shape
            
            nomArrayResponse = phased.ArrayResponse('SensorArray', obj.Array);
            theta = [0:180 179:-1:0];
            phi = [180*zeros(1, length(0:180)) zeros(1, length(179:-1:0))];
            azel = phitheta2azel([phi; theta]);
            fhz = obj.Element.Frequency*1e9;
            r = nomArrayResponse(fhz, azel);
            rss = pow2db(abs(r).^2/obj.NumElements);
            rss(rss < obj.Sensitivity) = obj.Sensitivity;
            gp = GainPattern(0, 0:360, rss');
        end
        
        function rot = fullRotation(obj, antennaPos, sources, thetaSteer)
            % TODO: implement a function to generate all the data that
            % would result from the full rotation of this antenna at the
            % given position with sources given by the InterferenceSource
            % objects (sources) with the antenna steered to the theta
            % angle, thetaSteer

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
        
        function [] = plot(obj)
            % plot  plot the array pattern as a normalized power dB pattern
            % using the phased array toolbox plotting function

            pattern(obj.Array, obj.Element.Frequency*1e9, -180:180, -90:90, ...
                    'PropagationSpeed', physconst('LightSpeed'), ...
                    'CoordinateSystem', 'polar', ...
                    'Type', 'powerdb', ...
                    'Normalize', true);
        end
        
        
    end
    
    
end