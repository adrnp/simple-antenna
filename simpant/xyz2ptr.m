function ptr = xyz2ptr(xyz)
% XYZ2PTR   coverts from the cartesian coordinate system to spherical
% coordinates.
%   ptr = XYZ2PTR(xyz) converts the vector, or matrix of column vectors,
%   XYZ to a phi, theta, r spherical vector representation.  The vector XYZ
%   must be a column vector or a matrix of 3xN column vectors to convert.
%   The resulting spherical coordinate vectors are defined with angles in
%   degrees.
%
%   SEE ALSO: PTR2XYZ


% extract the pieces for clarity
x = xyz(1,:);
y = xyz(2,:);
z = xyz(3,:);

% do the conversion
r = sqrt(x.^2 + y.^2 + z.^2);
thetar = acos(z./r);
phir = atan2(y, x);

% convert to degrees and combine
ptr = [wrapTo360(rad2deg(phir)); rad2deg(thetar); r];