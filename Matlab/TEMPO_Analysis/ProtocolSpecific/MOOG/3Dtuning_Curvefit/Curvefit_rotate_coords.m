%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_rotate_coords.m -- This module performs a rotation of spherical coordinates.
%-- Created - pwatkins, 4/21
%-----------------------------------------------------------------------------------------------------------------------
function [rott,rotp] = Curvefit_rotate_coords(t,p,t0,p0)

% (t,p) are the (azimuth,elevation) coordinates to be transformed.
% (t0,p0) are the (azimuth,elevation) rotation angles of the new axes 
% relative to the axes of the given coordinates.
num_points = length(t);  % should be same as length(p)

% NOTE - the definitions of azimuth and elevation used here are
% different from the standard spherical definitions of theta
% as the angle from the z axis in [0,pi] and phi as the angle
% from the x axis in the xy-plane in [0,2*pi).
% using this definition, the sphere is parameterized as follows:
%   x = sin(t)*cos(p)   y = sin(t)*sin(p)   z = cos(t)
% INSTEAD we define theta as the angle from the x axis in the
% xy plane in [-pi,pi) and phi as the elevation angle
% from the xy plane in [-pi/2, pi/2].
% using this definition, the sphere is parameterized as follows:
%   x = cos(p)*cos(t)   y = cos(p)*sin(t)   z = sin(p)
% This definition is consistent with the matlab sph2cart and cart2sph
% functions (see matlab help files for nice pictures).

% clockwise rotation about the y axis
% when looking from the positive y axis towards the origin.
roty = [ cos(p0) 0 -sin(p0);
           0     1    0;
         sin(p0) 0  cos(p0) ];

% clockwise rotation about the z axis
% when looking from the positive z axis towards the origin.
rotz = [  cos(t0) sin(t0) 0;
         -sin(t0) cos(t0) 0;
            0       0     1 ];

% rotation about the z-axis followed by a 
% rotation about the y-axis is roty*rotz.
% counterclockwise rotations are the transpose
% of the above matrices.

% change from spherical to cartesian coordinates
[x y z] = sph2cart(t,p,ones(1,num_points));

% use rotation directions consistent with how we defined spherical
% coordinates above.  positive azimuth angles move counterclockwise
% in the xy plane when viewed from the positive z axis towards the origin.
% positive elevation angles move clockwise in the xz plane when viewed
% from the positive y axis towards the origin.  therefore, we want
% a counterclockwise rotation about the z axis followed by a clockwise
% rotation about the y axis.

% Decided NOT to do this:
% rotate coordinate frame in cartesion coordinates
% using the transpose of the rotation matrix.
% rotating the coordinate frame means finding the coordinates
% of the original points with respect to the rotated axes.
%c_hat = (roty*rotz')'*[x;y;z];

% But INSTEAD to do this:
% I need to come up with a better way to explain this.
% In the meantime, feel free to email me for my current explanation.
c_hat = roty'*rotz*[x;y;z];
        
% change back from cartesion to spherical coordinates.
% rott and rotp are now the respective theta and phi coordinates in the
% rotated coordinate frame.
[rott rotp rotr] = cart2sph(c_hat(1,:), c_hat(2,:), c_hat(3,:));

% take care with the special case at the poles, set theta angle to 0.
e = 1e-10;  % floating point tolerance
select = abs(rotp - pi/2) < e | abs(rotp + pi/2) < e;
rott(select) = 0;

