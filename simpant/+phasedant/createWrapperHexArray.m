function array = createWrapperHexArray(Nmin, Nmax, element, varargin)

% use the hex grid spacing
dy = 1/2;
dz = sqrt(3)/4;

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
    pos(2,idx) = (-(rows(count)-1)/2:(rows(count)-1)/2)*dy;
    pos(3,idx) = m*dz;
end

array = phasedant.WrapperArrayAntenna(element, pos, varargin{:});