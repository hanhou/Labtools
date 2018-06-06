% models fitting for 3D tuning
% LBY 20170612

%% load data
clear all;
tic;
reps = 200;

% TEMPO_Defs;
% Path_Defs;
% ProtocolDefs;

stimType{1}='Vestibular';
stimType{2}='Visual';
stimType{3}='Combined';

models = {'VO','AO','VA','VJ','AJ','VAJ','PVAJ'};
% models = {'VA'};
models_color = {'r','b','m','y','c','k','k'}; % for figures
% models_color = {'r','b','m','k'}; % for figures

global PSTH3Dmodel PSTH FigPara;

Protocol = 2;
FnameCode = 2;

switch FnameCode
    case 1 % choose files mannully
        [filename pathname] = uigetfile('Z:\Data\TEMPO\BATCH');
        
    case 2 % load all mat files from interested folder
        switch Protocol
            case 1
                pathname = 'Z:\Data\TEMPO\BATCH\QQ_3D_Tuning_T_PCC';
                cd('Z:\Data\TEMPO\BATCH\QQ_PSTH3DModel_T_PCC');
                disp('Load PSTH T data...');
            case 2
                pathname = 'Z:\Data\TEMPO\BATCH\QQ_3D_Tuning_R_PCC';
                cd('Z:\Data\TEMPO\BATCH\QQ_PSTH3DModel_R_PCC');
                disp('Load PSTH R data...');
        end
end

filename = dir(pathname);
for c_inx = 1:length(filename)-2
    % load PSTH data
    temp = load([pathname '\' filename(2+c_inx).name]);
    temp_name = [filename(2+c_inx).name(1:10) '_PSTH'];
    eval([temp_name,'= temp.result.PSTH',';']);
    eval([temp_name,'.stimType = temp.result.unique_stimType',';']);
    eval([temp_name,'.cellName = temp.result.FILE',';']);
    eval([temp_name,'.nBins = temp.result.nBins',';']);
    eval([temp_name,'.pro = temp.result.Protocol',';']);
    eval([temp_name,'.meanSpon = temp.result.meanSpon',';']);
    disp([filename(2+c_inx).name(1:10), ' Pack PSTH data SUCCESSE!']);
    
    % pack data for models
    eval(['temp_PSTH = ', temp_name,'.spk_data_bin_mean_rate;']);
    eval(['temp_spatialData = ', temp_name,'.spk_data_count_mean_rate;']);
    eval(['temp_nBins = ', temp_name,'.nBins;']);
    eval(['temp_stimType = ', temp_name,'.stimType;']);
    eval(['temp_meanSpon = ', temp_name,'.meanSpon;']);
    eval(['temp_pro = ', temp_name,'.pro;']);
    eval(['temp_sigTrue = ', temp_name,'.respon_sigTrue;']);
    eval(['temp_spon = ', temp_name,'.spon_spk_data_bin_mean_rate;']);
    
    % fit models for those have significant responses
    for k = 1: length(temp_stimType)
        
        if temp_sigTrue(k) == 1
            models_fitting(models, models_color, temp_name, temp_pro, temp_stimType(k),temp_meanSpon, temp_PSTH{k}, temp_spatialData{k}, temp_spon,temp_nBins, reps);
        end
        
        eval([temp_name,'_modelData {temp_stimType(k)} = PSTH3Dmodel',';']);
        
        % data saving
        switch Protocol
            case 1
                temp_save_file = [stimType{temp_stimType(k)},'_',temp_name(1:end-5),'_model_T'];
            case 2
                temp_save_file = [stimType{temp_stimType(k)},'_',temp_name(1:end-5),'_model_R'];
        end
        eval(['model_para = ',temp_name,'_modelData{',num2str(temp_stimType(k)),'};']);
        eval(['save '  temp_save_file ' model_para']);
        
    end
    
end
toc;
