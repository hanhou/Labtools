%-----------------------------------------------------------------------------------------------------------------------
%-- DDiscrim_GetEyeData.m -- Returns conjugate and disconjugate eye position data, both uncalibrated and calibrated.
%--     Note: Calibration is done here using both fixation and saccade data from the discrimination run.
%--	GCD, 5/21/01
%-----------------------------------------------------------------------------------------------------------------------
function [h_verg, v_verg, h_conj, v_conj, calib_h_verg, calib_v_verg, calib_h_conj, calib_v_conj] = ...
    DDiscrim_GetEyeData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

PLOT_EYE_DATA = 0;

%compute various eye position metrics
h_verg = data.eye_positions(LEYE_H, :) - data.eye_positions(REYE_H, :);
v_verg = data.eye_positions(LEYE_V, :) - data.eye_positions(REYE_V, :);
h_conj = (data.eye_positions(LEYE_H, :) + data.eye_positions(REYE_H, :))/2;
v_conj = (data.eye_positions(LEYE_V, :) + data.eye_positions(REYE_V, :))/2;

%prepare to calibrate eye positions using the fixation and saccade data from this run
num_trials = size(data.event_data, 3);
T1_eye_posn = ComputeMeanEyePos(data, num_trials, IN_T1_WIN_CD, SUCCESS_CD, 3, 10);
T2_eye_posn = ComputeMeanEyePos(data, num_trials, IN_T2_WIN_CD, SUCCESS_CD, 3, 10);

lh = []; lv=[]; rh=[]; rv=[];
lh = [T1_eye_posn(LEYE_H, :) T2_eye_posn(LEYE_H, :) data.eye_positions(LEYE_H, :)]';
rh = [T1_eye_posn(REYE_H, :) T2_eye_posn(REYE_H, :) data.eye_positions(REYE_H, :)]';
lv = [T1_eye_posn(LEYE_V, :) T2_eye_posn(LEYE_V, :) data.eye_positions(LEYE_V, :)]';
rv = [T1_eye_posn(REYE_V, :) T2_eye_posn(REYE_V, :) data.eye_positions(REYE_V, :)]';

fh = [data.targ_params(TARG_XCTR, :, T1) data.targ_params(TARG_XCTR, :, T2) data.targ_params(TARG_XCTR, :, FP)]';
fv = [data.targ_params(TARG_YCTR, :, T1) data.targ_params(TARG_YCTR, :, T2) data.targ_params(TARG_YCTR, :, FP)]';

%remove NaNs from all arrays
NaN_trials = isnan(lh);
lh = lh(~NaN_trials);
rh = rh(~NaN_trials);
lv = lv(~NaN_trials);
rv = rv(~NaN_trials);
fh = fh(~NaN_trials);
fv = fv(~NaN_trials);

%[lh lv rh rv fh fv]
if (PLOT_EYE_DATA)
    figure;
    hold on;
    plot(lh, lv, 'ro');
    plot(rh, rv, 'bo');
    plot(fh, fv, 'g+');
    hold off;
end

HLB=[-10.0; 0.99; -10.0];
HUB=[10.0; 1.01; 10.0];
VLB=[-10.0; 0.1; -0.01];
VUB=[10.0; 10.0; 0.01];
[L_horiz_pars, L_vert_pars] = EyeCalib_Fitter(lh', lv', fh', fv', HLB, HUB, VLB, VUB);
[calib_lh, calib_lv] = ComputeCalibratedEyePosn(lh', lv', L_horiz_pars, L_vert_pars);
[R_horiz_pars, R_vert_pars] = EyeCalib_Fitter(rh', rv', fh', fv', HLB, HUB, VLB, VUB);
[calib_rh, calib_rv] = ComputeCalibratedEyePosn(rh', rv', R_horiz_pars, R_vert_pars);  

if (PLOT_EYE_DATA)
    figure;
    hold on;
    plot(calib_lh, calib_lv, 'ro');
    plot(calib_rh, calib_rv, 'bo');
    plot(fh, fv, 'g+');
    hold off;
end

%now use calibration parameters to compute new eye data to put into analysis
[calib_lh, calib_lv] = ComputeCalibratedEyePosn(data.eye_positions(LEYE_H, :), data.eye_positions(LEYE_V, :), L_horiz_pars, L_vert_pars);
data.eye_positions_calibrated(LEYE_H,:) = calib_lh;
data.eye_positions_calibrated(LEYE_V,:) = calib_lv;
[calib_rh, calib_rv] = ComputeCalibratedEyePosn(data.eye_positions(REYE_H, :), data.eye_positions(REYE_V, :), R_horiz_pars, R_vert_pars);
data.eye_positions_calibrated(REYE_H,:) = calib_rh;
data.eye_positions_calibrated(REYE_V,:) = calib_rv;

%compute various eye position metrics again
calib_h_verg = data.eye_positions_calibrated(LEYE_H, :) - data.eye_positions_calibrated(REYE_H, :);
calib_v_verg = data.eye_positions_calibrated(LEYE_V, :) - data.eye_positions_calibrated(REYE_V, :);
calib_h_conj = (data.eye_positions_calibrated(LEYE_H, :) + data.eye_positions_calibrated(REYE_H, :))/2;
calib_v_conj = (data.eye_positions_calibrated(LEYE_V, :) + data.eye_positions_calibrated(REYE_V, :))/2;

return;