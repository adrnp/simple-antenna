function [] = plot3DGainPattern(thetas, phis, pattern, minGain)


if nargin < 4
    minGain = -40;
end

pattern(pattern < minGain) = minGain;

x = -(minGain - pattern).*sind(thetas).*cosd(phis);
y = -(minGain - pattern).*sind(thetas).*sind(phis);
z = -(minGain - pattern).*cosd(thetas);


r = sqrt(x.^2 + y.^2 + z.^2);

surf(x, y, z, r);
shading interp
axis off
daspect([1 1 1])