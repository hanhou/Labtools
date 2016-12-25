function [horiz_pars, vert_pars] = EyeCalib3_Fitter(horiz_raw, vert_raw, horiz_desired, vert_desired, lower_bounds, upper_bounds)

%this function takes in the measured horiz/vert eye position for each trial and
%tries to minimize the difference between the observed values and the
%desired values (location of fixation target).  GCD, starting 11/27/00
%NOTE: EyeCalib2_Fitter differs from EyeCalib_Fitter by including a gain term for
%the interaction of the iso and ortho axes.  VR, 7/8/06
%NOTE: the inputs are assumed to be ROW vectors, hence the transpose!

global Eye_Data;

a(1) = 0;  %horiz: offset
a(2) = 1.0; % horiz: gain on horiz axis
a(3) = 0.0; % horiz: gain on vert axis
a(4) = 0.0; % horiz: gain on horiz-vert interaction
a(5) = 0;  %vert: offset
a(6) = 1.0; % vert: gain on vert axis
a(7) = 0.0; % vert: gain on horiz axis
a(8) = 0.0; % vert: gain on horiz-vert interaction

%first, minimize errors for the horizontal dimension
%fill eye data as follows: [iso_raw ortho_raw desired]
Eye_Data = [horiz_raw' vert_raw' horiz_desired' vert_desired'];

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=lower_bounds;
UB=upper_bounds;
OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);
pars = fmincon('EyeCalib3_Func',a,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
horiz_pars = pars(1:4);
vert_pars = pars(5:8);

return;
