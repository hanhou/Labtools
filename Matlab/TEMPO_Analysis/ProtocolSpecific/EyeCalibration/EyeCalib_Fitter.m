function [horiz_pars, vert_pars] = EyeCalib_Fitter(horiz_raw, vert_raw, horiz_desired, vert_desired, H_lower_bounds, H_upper_bounds, V_lower_bounds, V_upper_bounds)
%this function takes in the measured horiz/vert eye position for each trial and
%tries to minimize the difference between the observed values and the
%desired values (location of fixation target).  GCD, starting 11/27/00
%NOTE: input vectors must be ROW VECTORS

global Eye_Data;

a(1) = 0;  %offset
a(2) = 1.0; % gain on iso axis
a(3) = 0.0; % gain on ortho axis

%first, minimize errors for the horizontal dimension
%fill eye data as follows: [iso_raw ortho_raw desired]
Eye_Data = [horiz_raw' vert_raw' horiz_desired'];

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=H_lower_bounds;
UB=H_upper_bounds;
OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);
horiz_pars = fmincon('EyeCalib_Func',a,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);

a(1) = 0;  %offset
a(2) = 1.0; % gain on iso axis
a(3) = 0.0; % gain on ortho axis

%now, minimize errors for the vertical dimension
Eye_Data = [vert_raw' horiz_raw' vert_desired'];

LB=V_lower_bounds;
UB=V_upper_bounds;
vert_pars = fmincon('EyeCalib_Func',a,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);

return;
