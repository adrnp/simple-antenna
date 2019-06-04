function [] = plotGainPattern(angle, gain, rLimits)
% PLOTGAINPATTERN  plot a gain pattern given the angle of the gain and the
% gain measurement (in dB)
%
%   PLOTGAINPATTERN(ANGLE, GAIN) plots the gain values for each of the
%   angles of a proper polar axis.  ANGLE is in degrees and GAIN in is dB
%
% Examples:
%

% to get the bounds for r lim, we need to get the min value and max value
minval = min(gain);
maxval = max(gain);

if nargin < 3
    rLimits = [minval maxval];
end

%
% make the correct color data...
%

% need to do some funky sign changes and manipulation
startVal = -ceil(maxval);
endVal = -floor(minval);
numvals = length(startVal:0.1:endVal);

rawcd = [uint8(parula(numvals)*255) uint8(zeros(numvals, 1))];
rawcd = flipud(rawcd);
cd = uint8(zeros(length(gain), 4));
for i = 1:length(gain)
    idx = round(-gain(i)*10 - startVal*10 + 1);
    if idx < 1
        idx = 1;
    elseif idx > length(rawcd)
        idx = length(rawcd);
    end
    cd(i,:) = rawcd(idx,:);
end
cd = cd';

% simply do a polar plot
p = polarplot(deg2rad(angle), gain, 'LineWidth', 1.5);
drawnow;
set(p.Edge, 'ColorBinding', 'interpolated', 'ColorData', cd);

% make the axis correct -> this is the whole point of making this function
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
pax.RLim = rLimits;