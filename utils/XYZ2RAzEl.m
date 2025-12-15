function [R, Az, El] = XYZ2RAzEl(X, Y, Z)
% XYZ2RAzEl Convert Cartesian coordinates to spherical coordinates.
%   Inputs:
%       X, Y, Z - Cartesian coordinates.
%
%   Outputs:
%       R  : Radius (Euclidean distance)
%       Az : Azimuth angle in radians [rad], measured from the +X axis toward the +Y axis.
%       El : Elevation angle in radians [rad], measured from the XY-plane toward +Z.
%
%   Note:
%       The coordinate system is right-handed, with Z pointing upward.

    R = sqrt(X .* X + Y .* Y + Z .* Z);
    Az = atan2(Y, X);  % alpha
    El = asin(Z ./ R); % beta
end