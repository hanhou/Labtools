function     ZCFUNC(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
DataMatrixDef;
all_data = loadzcdata(FILE);
[allData{1},allData{2},allData{3}] = SepByStimType(all_data);
%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
SAVPATH = 'E:\Data\Gutou\Raw\';
FILENAME = [FILE(1:end-4) '.mat']; 
data_extract(good_data,FILENAME,SAVPATH);
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

a = [FILE_list.name]