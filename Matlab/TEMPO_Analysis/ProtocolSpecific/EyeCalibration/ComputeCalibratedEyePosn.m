function [horiz_calib, vert_calib] = ComputeCalibratedEyePosn(horiz_raw, vert_raw, horiz_pars, vert_pars)
% computes the calibrated eye positions from the raw data and the fit parameters

   horiz_calib = horiz_pars(1) + horiz_pars(2).*horiz_raw + horiz_pars(3).*vert_raw;
   vert_calib = vert_pars(1) + vert_pars(2).*vert_raw + vert_pars(3).*horiz_raw;
   
return;