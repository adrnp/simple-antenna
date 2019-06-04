function [] = plotRotationPattern(rot, method, polarBounds)
% PLOTROTATIONPATTERN   plot the rotation pattern information, this
% includes the gain pattern itself, the true bearing, and some of the
% bearing information.
%
%   PLOTROTATIONPATTERN(ROT) plots the results of the 3dB method with
%   default bounds for the given GPRotation
%
%   PLOTROTATIONPATTERN(ROT, METHOD) specifies the method to use for
%   bearing determination
%
%   PLOTROTATIONPATTERN(ROT, METHOD, POLARBOUNDS) also specifies the polar
%   bounds for the plot

% get the some of the elements for plotting the gain pattern
phi = rot.GP.Phi;
ss = rot.GP.Strength;

if nargin < 2
    method = '3db';
end

if nargin < 3
    minval = min(ss);
    maxval = max(ss);
    polarBounds = [minval maxval];
end

trueBearing = rot.AverageTrueBearing();
% TODO: maybe want to also plat the range of the bearing for duration of
% the rotation...

% calculate the bearing
[b, p, e] = rot.GP.getBearing(method);
Nlobes = length(b);

% plot pattern and true bearing
rot.GP.plot('RLim', polarBounds); hold on;
polarplot(deg2rad([trueBearing trueBearing]), polarBounds, '--');

% plot the lobe bounds
for i = 1:Nlobes
    setData = e.sets{i};
    angs = setData(:,1);
    polarplot(deg2rad(angs), [-70 -40*ones(1, length(angs)-2) -70]);
end

% plot the bearing lines + probability as text
co = get(gca,'colororder');
[nc, ~] = size(co);
for i = 1:Nlobes
    colorInd = 2+i;
    if colorInd > nc
        colorInd = colorInd - nc;
    end
    polarplot(deg2rad([b(i) b(i)]), polarBounds, ...
        'Color', co(colorInd,:));

    text(deg2rad(b(i)), -35, sprintf('%0.2f', p(i)), ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'Color', co(colorInd,:));
end
hold off;