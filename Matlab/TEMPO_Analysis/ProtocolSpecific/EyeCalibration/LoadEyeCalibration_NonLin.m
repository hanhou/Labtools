function [caldata, doneflag] = LoadEyeCalibration_NonLin(data, PATH, FILE);

TEMPO_Defs;

%check if calibration file exists, and read calibration params in from file
i = size(PATH,2) - 1;
while PATH(i) ~='/'	%Analysis directory is one branch below Raw Data Dir
    i = i - 1;
end   
PATHIN = [PATH(1:i) 'Analysis/Eye_Calibration/'];

run_loc = find(FILE == 'r');
file_root_name = FILE(1:run_loc-1);
in_name = [PATHIN file_root_name '_nonlin_eye_calib.mat'];

if exist(in_name)
    buff = sprintf('loading eye calibration data from %s', in_name);
    disp(buff);
    Pars = load(in_name);
    %Note: Pars will be a structure, Pars.M is the loaded matrix
    
    %compute and store the calibrated positions for the left eye
    [clh, clv] = ComputeCalibratedEyePosn_NonLin(data.eye_positions(LEYE_H,:), data.eye_positions(LEYE_V,:), Pars.M(:, LEYE_H)', Pars.M(:, LEYE_V)' );  
    caldata(LEYE_H, :) = clh;
    caldata(LEYE_V, :) = clv;

    %compute and store the calibrated positions for the right eye
    [crh, crv] = ComputeCalibratedEyePosn_NonLin(data.eye_positions(REYE_H,:), data.eye_positions(REYE_V,:), Pars.M(:, REYE_H)', Pars.M(:, REYE_V)' );  
    caldata(REYE_H, :) = crh;
    caldata(REYE_V, :) = crv;
    doneflag = 1;
    disp('calibrated eye signals stored in data.eye_positions_calibrated()');
else
    caldata = [];
    doneflag = 0;
end

return;