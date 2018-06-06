% for 3D tuning with 26 directions in 3D space
% fig.2 plot PSTH for each trial and raster plot across directions
% fig.3 plot mean PSTHs across directions
% fig.4 plot coutour tuning responses (sum response of total 1.5s)
% fig.5 PSTH - spon
% fig.10 3D model
%
% LBY,201612



function Tuning_3D_analysis_LBY_old(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag)

tic;

% contains protocol etc. specific keywords
TEMPO_Defs;
Path_Defs;
ProtocolDefs;

stimType{1}='Vestibular';
stimType{2}='Visual';
stimType{3}='Combined';

global PSTH3Dmodel PSTH;
PSTH3Dmodel = [];PSTH = [];

%% get data
% stimulus type,azi,ele,amp and duration

% temporarily, for htb missing cells
if data.one_time_params(NULL_VALUE) == 0
    data.one_time_params(NULL_VALUE) = -9999;
end

switch Protocol
    case DIRECTION_TUNING_3D % for translation
        trials = 1:size(data.moog_params,2);
        % temp_spon_trials = trials( (data.moog_params(AZIMUTH,trials,MOOG) == data.one_time_params(NULL_VALUE)) ); % all sponteneous trials
        
        % find the trials for analysis
        temp_trials = trials(find( (data.moog_params(AZIMUTH,trials,MOOG) ~= data.one_time_params(NULL_VALUE)) )); % 26*rept trials
        temp_stimType = data.moog_params(STIM_TYPE,temp_trials,MOOG);
        unique_stimType = munique(temp_stimType');
        BegTr = ceil((BegTrial-1)/((length(unique_stimType)*26+1)))*(length(unique_stimType)*26+1)+1;
        EndTr =  floor(EndTrial/((length(unique_stimType)*26+1)))*(length(unique_stimType)*26+1);
        %     BegTr = ceil((BegTrial-1)/((length(unique_stimType)*26)))*(length(unique_stimType)*26)+1; % for no-spon trials
        %     EndTr =  floor(EndTrial/((length(unique_stimType)*26)))*(length(unique_stimType)*26); % for no-spon trials
        select_trials = find( (trials >= BegTr) & (trials <= EndTr) );
        spon_trials = select_trials(data.moog_params(AZIMUTH,select_trials,MOOG) == data.one_time_params(NULL_VALUE));
        real_trials = select_trials(data.moog_params(AZIMUTH,select_trials,MOOG) ~= data.one_time_params(NULL_VALUE));
        
        temp_azimuth = data.moog_params(AZIMUTH,real_trials,MOOG);
        temp_elevation = data.moog_params(ELEVATION,real_trials,MOOG);
        temp_amplitude = data.moog_params(AMPLITUDE,real_trials,MOOG);
        
        
    case ROTATION_TUNING_3D % for rotation
        trials = 1:size(data.moog_params,2);
        
        % find the trials for analysis
        temp_trials = trials(find( (data.moog_params(ROT_AZIMUTH,trials,MOOG) ~= data.one_time_params(NULL_VALUE)) )); % 26*rept trials
        temp_stimType = data.moog_params(STIM_TYPE,temp_trials,MOOG);
        unique_stimType = munique(temp_stimType');
        BegTr = ceil((BegTrial-1)/((length(unique_stimType)*26+1)))*(length(unique_stimType)*26+1)+1;
        EndTr =  floor(EndTrial/((length(unique_stimType)*26+1)))*(length(unique_stimType)*26+1);
        select_trials = find( (trials >= BegTr) & (trials <= EndTr) );
        spon_trials = select_trials(data.moog_params(ROT_AZIMUTH,select_trials,MOOG) == data.one_time_params(NULL_VALUE));
        real_trials = select_trials(data.moog_params(ROT_AZIMUTH,select_trials,MOOG) ~= data.one_time_params(NULL_VALUE));
        
        temp_azimuth = data.moog_params(ROT_AZIMUTH,real_trials,MOOG);
        temp_elevation = data.moog_params(ROT_ELEVATION,real_trials,MOOG);
        temp_amplitude = data.moog_params(ROT_AMPLITUDE,real_trials,MOOG);
        
end

temp_duration = data.moog_params(DURATION,real_trials,MOOG); % in ms
temp_stimType = data.moog_params(STIM_TYPE,real_trials,MOOG);

unique_azimuth = munique(temp_azimuth');
unique_elevation = munique(temp_elevation');
unique_amplitude = munique(temp_amplitude');
unique_duration = munique(temp_duration');
unique_stimType = munique(temp_stimType');

% the least repitition numbers
repNums = (EndTr - BegTr +1)/(length(unique_stimType)*26+1);
% repNums = (EndTr - BegTr +1)/(length(unique_stimType)*26); % for no-spon trials
% time information?
eye_timeWin = 1000/(data.htb_header{EYE_DB}.speed_units/data.htb_header{EYE_DB}.speed/(data.htb_header{EYE_DB}.skip+1)); % in ms
spike_timeWin = 1000/(data.htb_header{SPIKE_DB}.speed_units/data.htb_header{SPIKE_DB}.speed/(data.htb_header{SPIKE_DB}.skip+1)); % in ms
event_timeWin = 1000/(data.htb_header{EVENT_DB}.speed_units/data.htb_header{EVENT_DB}.speed/(data.htb_header{EVENT_DB}.skip+1)); % in ms

event_in_bin = squeeze(data.event_data(:,:,trials));
stimOnT = find(event_in_bin == VSTIM_ON_CD); 
stimOffT = find(event_in_bin == VSTIM_OFF_CD);
FPOnT = find(event_in_bin == FP_ON_CD);
FixInT = find(event_in_bin == IN_FIX_WIN_CD);

% spike data
% 5000*265 (for 5 repetitions)
spike_data = squeeze(data.spike_data(SpikeChan,:,real_trials)); % sipke data in 5000ms for each trial
spike_data( spike_data > 100 ) = 1; % something is absolutely wrong
spon_spk_data = squeeze(data.spike_data(SpikeChan,:,spon_trials));


%% calculate data for drawing and modeling
% k=1,2,3
% j= -90,-45,0,45,90 (up->down)
% i=0 45 90 135 180 225 270 315
% for PSTH, from 500ms before stim on to 500ms after stim off


% for models
timeWin = 300; % in ms
timeStep = 25; % in ms
gau_sig = 100;

% for PSTH analysis
% timeWin = 400; % in ms
% timeStep = 30; % in ms

tOffset1 = 200; % in ms
tOffset2 = 200; % in ms
% align to stim on
nBins = floor((temp_duration(1)+tOffset2+timeStep/2)/timeStep)+floor((tOffset1-timeStep/2)/timeStep); % -500:2000ms
PSTH_onT = stimOnT(1) - floor((tOffset1-timeStep/2)/timeStep)*timeStep;
stimOnBin = (stimOnT(1)-PSTH_onT+timeStep)/timeStep;
stimOffBin = (stimOffT(1)-PSTH_onT+timeStep)/timeStep;

% spk rates of real trials (no spon)
for k = 1:length(unique_stimType)
    pc = 0;
    PSTH.spk_data_bin_rate_aov{k} = [];% for PSTH ANOVA
    PSTH.spk_data_count_rate{k} = [];%
    PSTH.spk_data_bin_mean_rate{k} = []; % mean spike count per time window (for PSTH) % 1*nbins
    PSTH.spk_data_bin_mean_rate_std{k} = [];
    PSTH.spk_data_count_mean_rate{k} = [];
    PSTH.spk_data_bin_mean_rate_ste{k} = [];
    for j = 1:length(unique_elevation)
        for i = 1:length(unique_azimuth)
            spk_data{k,j,i} = []; % raw spike count
            
            PSTH.spk_data_bin_rate{k,j,i} = []; % spike counts per time window (for PSTH across trials)
            
            select = find( (temp_azimuth == unique_azimuth(i)) & (temp_elevation == unique_elevation(j)) & (temp_stimType == unique_stimType(k))) ;
            if sum(select)>0
                spk_data{k,j,i}(:,:) = spike_data(:,select);
                dim = size(spk_data{k,j,i}(:,:));
                pc = pc+1;
                spk_data_count_rate_anova{k}(:,pc) =  sum(spk_data{k,j,i}(stimOnT(1):stimOffT(1),:),1)/(unique_duration(1,1)/1000);
                resp_sse(k,pc) = sum((spk_data_count_rate_anova{k}(:,pc) - mean(spk_data_count_rate_anova{k}(:,pc))).^2);
                resp_trialnum(k,pc)= length(spk_data_count_rate_anova{k}(:,pc));
            else
                spk_data{k,j,i}(:,:) = spk_data{k,j,1}(:,:);
            end
            % for PSTH
            try
            PSTH.spk_data_bin_rate{k,j,i} = PSTH_smooth( nBins, PSTH_onT, timeWin, timeStep, spk_data{k,j,i}(:,:), 2, gau_sig);
            catch 
                keyboard;
            end
            PSTH.spk_data_bin_rate_aov{k}(pc,:,:) = PSTH_smooth( nBins, PSTH_onT, timeWin, timeStep, spk_data{k,j,i}(:,:), 2, gau_sig);
            
            %             for nn = 1:nBins(1,1)
            %                 PSTH.spk_data_bin_rate{k,j,i}(nn,:) = sum(spk_data{k,j,i}(PSTH_onT-timeWin/2+timeStep*(nn-1):PSTH_onT+timeWin/2+timeStep*(nn-1),:),1)/(timeWin/1000);
            %
            %                 try
            %                     PSTH.spk_data_bin_rate_aov{k}(pc,nn,:) = sum(spk_data{k,j,i}(PSTH_onT-timeWin/2+timeStep*(nn-1):PSTH_onT+timeWin/2+timeStep*(nn-1),:),1)/(timeWin/1000);
            %                 catch
            %                     keyboard;
            %                 end
            %             end
            
            PSTH.spk_data_count_rate{k}(j,i,:) = sum(spk_data{k,j,i}(stimOnT(1):stimOffT(1),:),1)/(unique_duration(1,1)/1000); % rates
            PSTH.spk_data_bin_mean_rate{k}(j,i,:) = mean(PSTH.spk_data_bin_rate{k,j,i},2);
            PSTH.spk_data_bin_vector{k} = PSTH.spk_data_bin_mean_rate{k};
            PSTH.spk_data_bin_vector{k}([1,5],2:end,:) = 0;
            PSTH.spk_data_bin_mean_rate_std{k}(j,i,:) = std(PSTH.spk_data_bin_rate{k,j,i},0,2);
            PSTH.spk_data_bin_mean_rate_ste{k}(j,i,:) = PSTH.spk_data_bin_mean_rate_std{k}(j,i,:)/sqrt(size(PSTH.spk_data_bin_rate{k,j,i}(1,:),2));
            PSTH.spk_data_count_mean_rate{k}(j,i) = mean(PSTH.spk_data_count_rate{k}(j,i,:),3);% for countour plot and ANOVA
            PSTH.spk_data_vector{k} = PSTH.spk_data_count_mean_rate{k};
            PSTH.spk_data_vector{k}([1,5],2:end) = 0;
            
        end
    end
    PSTH.time_profile{k} = squeeze(sum(sum(PSTH.spk_data_bin_mean_rate{k}(j,i,:),1),2));
    p_anova_dire(k) = anova1(spk_data_count_rate_anova{k},'','off'); % tuning anova
    maxSpkRealMean(k) = max(PSTH.spk_data_count_mean_rate{k}(:));
    minSpkRealMean(k) = min(PSTH.spk_data_count_mean_rate{k}(:));
    PSTH.maxSpkRealBinMean(k) = max(PSTH.spk_data_bin_mean_rate{k}(:));
    PSTH.minSpkRealBinMean(k) = min(PSTH.spk_data_bin_mean_rate{k}(:));
    PSTH.maxSpkRealBinAll(k) = max(PSTH.spk_data_bin_rate_aov{k}(:));
    PSTH.maxSpkRealBinMeanSte(k) = max(PSTH.spk_data_bin_mean_rate_ste{k}(:));
    
    resp_std(k) = sum(resp_sse(k,:))/(sum(resp_trialnum(k,:))-26);
    DDI(k) = (maxSpkRealMean(k)-minSpkRealMean(k))/(maxSpkRealMean(k)-minSpkRealMean(k)+2*sqrt(resp_std(k)));
    [Azi, Ele, Amp] = vectorsum(squeeze(PSTH.spk_data_vector{k}(:,:)));
    preferDire{k} = [Azi, Ele, Amp];
    
    
end

% sponteneous spk rates
spon_spk_count_rate = sum(spon_spk_data(stimOnT(1):stimOffT(1),:),1)/(unique_duration(1,1)/1000);
% for nn = 1:nBins(1,1)
%     PSTH.spon_spk_data_bin_rate(nn,:) = sum(spon_spk_data(PSTH_onT-timeWin/2+timeStep*(nn-1):PSTH_onT+timeWin/2+timeStep*(nn-1),:),1)/(timeWin/1000);
% end
PSTH.spon_spk_data_bin_rate = PSTH_smooth( nBins, PSTH_onT, timeWin, timeStep, spon_spk_data, 2, gau_sig);
PSTH.spon_spk_data_bin_mean_rate = mean(PSTH.spon_spk_data_bin_rate(:,:),2);
PSTH.spon_spk_data_bin_mean_rate_std = std(PSTH.spon_spk_data_bin_rate(:,:),0,2);
PSTH.spon_spk_data_bin_mean_rate_ste = PSTH.spon_spk_data_bin_mean_rate_std/sqrt(size(PSTH.spon_spk_data_bin_rate(:,:),2));

% for k = 1:length(unique_stimType)
%
%     p_wilcox_respon(k) = ranksum(squeeze(PSTH.spk_data_bin_rate_aov{k}(pc,ii,:))',PSTH.spon_spk_data_bin_rate(ii,:)') % if response
%
% end

% maximum spon value
PSTH.maxSpkSponAll = max(max(PSTH.spon_spk_data_bin_rate)); % for PSTH across trials
PSTH.maxSpkSponBinMean = max(PSTH.spon_spk_data_bin_mean_rate); % for PSTH mean
maxSpkSponMean = max(spon_spk_count_rate); % for contour
PSTH.maxSpkSponBinMeanSte = max(PSTH.spon_spk_data_bin_mean_rate_ste);
% mean spon value
PSTH.meanSpkSponBinMean = mean(PSTH.spon_spk_data_bin_mean_rate);% for PSTH mean & PSTH across trials
meanSpkSponMean = mean(spon_spk_count_rate);% for contour
meanSpon = meanSpkSponMean;


%% Analysis

% find the max Firing Rate (FR) of every time bin

for k = 1:length(unique_stimType)
    preDirAcrossTBin{k} = [];
    angleDiffs{k} = [];
    for nn = 1:nBins(1,1)
        PSTH.spk_data_bin_mean_rate{k}(j,i,nn);
        [Azi, Ele, Amp] = vectorsum(squeeze(PSTH.spk_data_bin_vector{k}(:,:,nn)));
        preDirAcrossTBin{k}(nn,:) = [Azi, Ele, Amp];
        for j = 1:length(unique_elevation)
            for i = 1:length(unique_azimuth)
                try
                    angleDiffs{k}(j,i,nn) = angleDiff(unique_azimuth(i),unique_elevation(j),0.11,preDirAcrossTBin{k}(nn,1),preDirAcrossTBin{k}(nn,2),preDirAcrossTBin{k}(nn,3));
                catch
                    keyboard;
                end
            end
            
        end
    end
end
%% analysis
for k = 1:length(unique_stimType)
    spk_datda_bin_rate_mean_minusSpon{k} = mean(PSTH.spk_data_bin_rate_aov{k},3)-repmat(PSTH.spon_spk_data_bin_mean_rate',size(mean(PSTH.spk_data_bin_rate_aov{k},3),1),1);
    spk_datda_bin_rate_minusSpon{k} = PSTH.spk_data_bin_rate_aov{k}-permute(repmat(PSTH.spon_spk_data_bin_rate,1,1,size(PSTH.spk_data_bin_rate_aov{k},1)),[3 1 2]);
    % significant temporal modulation (response to the stimulus)
    for pc = 1:size(PSTH.spk_data_bin_rate_aov{k},1)
        PSTH.sigBin{k,pc} = []; % find if this bin is sig.
        PSTH.sigTrue{k}(pc) = 0; % PSTH.sigTrue(pc) -> sig. of this direction
        PSTH.respon_sigTrue(k) = 0; % PSTH.sigTrue(pc) -> sig. of this cell
        PSTH.localPeak{k,pc} = []; % find the local peak bin of each direction
        PSTH.localTrough{k,pc} = []; % find the local trough bin of each direction
        PSTH.s{k,pc} = []; % just for check
        for nn = 3 : floor(stimOffBin - stimOnBin)-3
            PSTH.s{k,pc}(nn) = 0;
            % if 5 consecutive bins are sig.,then we say this bin is sig.
            for ii = nn-2:nn+2
                if ranksum(squeeze(PSTH.spk_data_bin_rate_aov{k}(pc,ii,:))',PSTH.spon_spk_data_bin_rate(ii,:)') < 0.05
                    PSTH.s{k,pc}(nn) =  PSTH.s{k,pc}(nn)+1;
                end
            end
            if  PSTH.s{k,pc}(nn) == 5
                PSTH.sigBin{k,pc}= [PSTH.sigBin{k,pc};nn]; % find the bin with significant response to the stim.
            end
            if ~isempty(PSTH.sigBin{k,pc})
                PSTH.sigTrue{k}(pc) = 1; % PSTH.sigTrue(pc) -> sig. of this direction
                for ii = 1:length(PSTH.sigBin{k,pc})
                    pp = 0;tt = 0;
                    if (spk_datda_bin_rate_mean_minusSpon{k}(pc,nn)>=spk_datda_bin_rate_mean_minusSpon{k}(pc,nn-1)) &&(spk_datda_bin_rate_mean_minusSpon{k}(pc,nn)>=spk_datda_bin_rate_mean_minusSpon{k}(pc,nn-1))
                        pp = pp+1;
                        PSTH.localPeak{k,pc}(pp,1)= nn; % indicate the bin num.
                        PSTH.localPeak{k,pc}(pp,2) = spk_datda_bin_rate_mean_minusSpon{k}(pc,nn); % indicate the mean-value of this peak
                    else if (spk_datda_bin_rate_mean_minusSpon{k}(pc,nn)<=spk_datda_bin_rate_mean_minusSpon{k}(pc,nn-1)) &&(spk_datda_bin_rate_mean_minusSpon{k}(pc,nn)<=spk_datda_bin_rate_mean_minusSpon{k}(pc,nn-1))
                            tt = tt+1;
                            PSTH.localTrough{k,pc}(tt,1) = nn; % indicate the bin num.
                            PSTH.localTrough{k,pc}(tt,2) = spk_datda_bin_rate_mean_minusSpon{k}(pc,nn); % indicate the mean-value of this trough
                        end
                    end
                end
            end
            
        end
    end
    if sum(PSTH.sigTrue{k}(:)) ~=0
        PSTH.respon_sigTrue(k) = 1;
    end
    PSTH.sig(k) = sum(PSTH.sigTrue{k}(:));
end

% PSTH.localPeak = cellfun( @(x) sortrows(x),PSTH.localPeak,'UniformOutput', false);
% PSTH.localTrough = cellfun( @(x) sortrows(x),PSTH.localTrough,'UniformOutput', false);

%% plot figures
%
% k=1,2,3
% j= -90,-45,0,45,90 (up->down)
% i=0 45 90 135 180 225 270 315

% 270-225-180-135-90-45-0-315-270 for figures
iAzi = [7 6 5 4 3 2 1 8 7];

% initialize default properties
set(0,'defaultaxesfontsize',24);
colorDefsLBY;

% the lines markers
markers = {
    % markerName % markerTime % marker bin time % color
    'FPOnT',FPOnT(1),(FPOnT(1)-PSTH_onT+timeStep)/timeStep,colorDGray;
    'stim_on',stimOnT(1),(stimOnT(1)-PSTH_onT+timeStep)/timeStep,colorDRed;
    'stim_off',stimOffT(1),(stimOffT(1)-PSTH_onT+timeStep)/timeStep,colorDRed;
    'aMax',stimOnT(1)+731,(stimOnT(1)+719-PSTH_onT+timeStep)/timeStep,colorDBlue;
    'aMin',stimOnT(1)+1099,(stimOnT(1)+1074-PSTH_onT+timeStep)/timeStep,colorDBlue;
    };


% preferDirectionOfTime;
% CosineTuningPlot;
% PSTH_3D_Tuning; % plot PSTHs across sessions;
% Contour_3D_Tuning; % plot countour figures;
% models_figure; % plot figures for 3D models



%% Data Saving

% Reorganized. HH20141124
config.batch_flag = batch_flag;

% Output information for test. HH20160415
if isempty(batch_flag)
    config.batch_flag = 'test.m';
    disp('Saving results to \batch\test\ ');
end

%%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = PackResult(FILE, PATH, SpikeChan, repNums, unique_stimType,Protocol, ... % Obligatory!!
    unique_azimuth, unique_elevation, unique_amplitude, unique_duration,...   % paras' info of trial
    markers,...
    timeWin, timeStep, tOffset1, tOffset2,nBins, ... % PSTH slide window info
    meanSpon, p_anova_dire, DDI,preferDire,PSTH,... % PSTH and mean FR info
    PSTH3Dmodel); % model info

switch Protocol
    case DIRECTION_TUNING_3D
        config.suffix = 'PSTH_T';
    case ROTATION_TUNING_3D
        config.suffix = 'PSTH_R';
end
config.xls_column_begin = 'meanSpon';
config.xls_column_end = 'SigDireNum_vis';

% figures to save
config.save_figures = [];

% Only once
config.sprint_once_marker = {'0.2f'};
config.sprint_once_contents = 'result.meanSpon';

% loop across stim_type
config.sprint_loop_marker = {{'0.0f','0.0f','0.2f','g'};
    {'d','d'}};
config.sprint_loop_contents = {'result.preferDire{k}(1), result.preferDire{k}(2),result.DDI(k),result.p_anova_dire(k)';
    'result.PSTH.respon_sigTrue(k), result.PSTH.sig(k)'};

config.append = 1; % Overwrite or append




SaveResult(config, result);




toc;
end