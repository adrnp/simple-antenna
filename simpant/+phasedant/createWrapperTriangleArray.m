function array = createWrapperTriangleArray(element, varargin)
% createWrapperTriangleArray     creates a 3 element array using the matlab
% hpased array toolbox.

% need to build the element position matrix
dy = 1/4;
dz = sqrt(3)/4;

% TODO: this is not the same location of the origin for the antenna as I
% used in the simpant custom toolbox element stuff

e1 = [0 0 0]';
% e2 = [0 -dy -dz]';  % try this as [0 dy dz]
e2 = [0 dy dz]';  % try this as [0 dy dz]
% e3 = [0 dy -dz]';   % try this as [0 2*dy 0]
e3 = [0 2*dy 0]';   % try this as [0 2*dy 0]

elementPos = [e1 e2 e3];
array = phasedant.WrapperArrayAntenna(element, elementPos, varargin{:});


% % convert the element to a phased array element element
% phasedElement = element.convertToMatlabElement();
% 
% 
% % compute the wavelength for the spacing of the antennas
% c = physconst('LightSpeed');                % speed of light in [m/s]
% lambda = c/(phasedElement.Frequency*1e9);   % the wavelength of the tuned frequency in [m]
% 
% 
% % need to build the element position matrix
% dy = 1/4 * lambda;
% dz = sqrt(3)/4 * lambda;
% 
% % TODO: this is not the same location of the origin for the antenna as I
% % used in the simpant custom toolbox element stuff
% 
% e1 = [0 0 0]';
% e2 = [0 -dy -dz]';  % try this as [0 dy dz]
% e3 = [0 dy -dz]';   % try this as [0 2*dy 0]
% 
% elementPos = [e1 e2 e3];
% array = phased.ConformalArray('ElementPosition', elementPos, ...
%                               'Element', phasedElement);