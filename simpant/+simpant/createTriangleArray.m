function array = createTriangleArray(element, varargin)
% createTriangleArray   create an array with 3 elements laid out in the
% triangle pattern.  This effectively recreate the BeaSt antenna.
%   array = simpant.createTriangleArray(element) creates a model of the
%   beast antenna with each of the antenna elements being defined by
%   element, which is of type simpant.element.*
%

% use the hex grid spacing
dx = 1/4;
dy = sqrt(3)/4;
% dy = 1/2;

% manually create the position matrix
e1 = [0 0 0]';
e2 = [dx dy 0]';
e3 = [2*dx 0 0]';
pos = [e1 e2 e3];

% create the antenna array using the array antenna class
array = simpant.ArrayAntenna(element, pos, varargin{:});