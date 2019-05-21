function array = createHexArray(Nmin, Nmax, element, varargin)
% createHexArray    create an ArrayAntenna with elements arranged in a
% hexagonal pattern.
%   array = simpant.createHexArray(Nmin, Nmax, element) create a hexagonal
%   pattern with the lower row containing Nmin elements and the middle row
%   containing Nmax elements.  The array will have each element be an
%   antenna defined by element (of type simpant.element.*)
%
%   NOTE: good choices are (2,3) and (3,5) to create 7 and 19 element
%   arrays with good patterns.


% use the hex grid spacing
dx = 1/2;
d2 = sqrt(3)/4;

% create the position of all the elements for this
rows = [Nmin:Nmax Nmax-1:-1:Nmin];
N = sum(rows);      % Total number of elements
stop = cumsum(rows);
start = stop-rows+1;
pos = zeros(3,N);
count = 0;
for m = Nmin-Nmax:Nmax-Nmin
    count = count+1;
    idx = start(count):stop(count);
    pos(1,idx) = (-(rows(count)-1)/2:(rows(count)-1)/2)*dx;
    pos(2,idx) = m*d2;
end

% create the antenna array using the array antenna class
array = simpant.ArrayAntenna(element, pos, varargin{:});