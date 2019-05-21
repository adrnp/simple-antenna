function xyz = ptr2xyz(ptr)
% PTR2XYZ   converts a spherical vector of phi, theta, r in [deg, deg,
% length] to cartesian XYZ coordinates.
%   xyz = PTR2XYZ(ptr) converts from the spherical coordinates to cartesian
%   coordinates.  The vector PTR must be a column vector or a matrix of
%   3xN column vectors to convert with angles defined in degrees.
%
%   NOTE: this uses a different convention than MATLAB does for angles in
%   the spherical coordinate system.  For this coordinate system, phi is
%   defined as the counter-clockwise angle in the XY plane from the X axis
%   (think compass heading) and theta is defined as the angle from the Z
%   axis to the XY plane.  Phi is defined on [0, 360) and theta is deifned
%   on [0, 180).
%
%   SEE ALSO: XYZ2PTR

% extract the parts and convert to radians
phi = deg2rad(ptr(1,:));
theta = deg2rad(ptr(2,:));
r = ptr(3,:);

% do the conversion
x = r.*sin(theta).*cos(phi);
y = r.*sin(theta).*sin(phi);
z = r.*cos(theta);

% combine the information
xyz = [x; y; z];