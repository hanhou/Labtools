% v = 32;
% disp = 1;
% depth = 6.0759;
% io = 3.5;

v = 1.5;
disp = -(180*3/(pi*1.5));
depth = -1.49999;
io = 3;

newdisp = 180*io*depth/(pi*v*(v+depth))
newdepth = 180*io*v/(-disp*pi*v+180*io)-v