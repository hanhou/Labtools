% to calculate the angle between two directions in 3D space
% azi -> azimuth
% ele -> elevation
% amp -> amplitude
% 20170324 LBY

function [angleDiff_xy, angleDiff_xz,angleDiff_yz]= angleDiff_2D(azi1,ele1,amp1,azi2,ele2,amp2)

% transform spherical coordinates to cartesian
[x1,y1,z1] = sph2cart(azi1*pi/180,ele1*pi/180,amp1);
[x2,y2,z2] = sph2cart(azi2*pi/180,ele2*pi/180,amp2);

% project vectors to x-y, x-z, y-z planes respectively
% and, calculate the angles between the two vectors

mod1_xy = sqrt(x1.^2+y1.^2);
mod2_xy = sqrt(x2.^2+y2.^2);
if (x1*y2 - x2*y1)>0
    angleDiff_xy = acos((x1.*x2+y1.*y2)/(mod1_xy*mod2_xy))*180/pi;
else
    angleDiff_xy = -acos((x1.*x2+y1.*y2)/(mod1_xy*mod2_xy))*180/pi;
end

mod1_xz = sqrt(x1.^2+z1.^2);
mod2_xz = sqrt(x2.^2+z2.^2);
if (x1*z2 - x2*z1)>0
    angleDiff_xz = acos((x1.*x2+z1.*z2)/(mod1_xz*mod2_xz))*180/pi;
else
    angleDiff_xz = -acos((x1.*x2+z1.*z2)/(mod1_xz*mod2_xz))*180/pi;
end


mod1_yz = sqrt(y1.^2+z1.^2);
mod2_yz = sqrt(y2.^2+z2.^2);
if (y1*z2 - y2*z1)>0
    angleDiff_yz = acos((y1.*y2+z1.*z2)/(mod1_yz*mod2_yz))*180/pi;
else
    angleDiff_yz = -acos((y1.*y2+z1.*z2)/(mod1_yz*mod2_yz))*180/pi;
end

end