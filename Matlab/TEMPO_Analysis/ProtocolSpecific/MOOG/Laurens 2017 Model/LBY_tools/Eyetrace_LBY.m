% protocol: 3D tuning
% plot eyetrace for 3D tuning
% LBY 20161221

%-----------------------------------------------------------------------------------------------------------------------
function Eyetrace_LBY(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag)

TEMPO_Defs;
Path_Defs;
ProtocolDefs;
%% get data
% stimulus type,azi,ele,amp and duration

trials = 1:size(data.moog_params,2);
% find the trials for analysis
temp_trials = trials(find( (data.moog_params(AZIMUTH,trials,MOOG) ~= data.one_time_params(NULL_VALUE)) )); % 26*rept trials
temp_stimType = data.moog_params(STIM_TYPE,temp_trials,MOOG);
unique_stimType = munique(temp_stimType');
BegTr = ceil((BegTrial-1)/((length(unique_stimType)*26+1)))*(length(unique_stimType)*26+1)+1;
EndTr =  floor(EndTrial/((length(unique_stimType)*26+1)))*(length(unique_stimType)*26+1);
select_trials = find( (trials >= BegTr) & (trials <= EndTr) );
spon_trials = select_trials(data.moog_params(AZIMUTH,select_trials,MOOG) == data.one_time_params(NULL_VALUE));
real_trials = select_trials(data.moog_params(AZIMUTH,select_trials,MOOG) ~= data.one_time_params(NULL_VALUE));
%
temp_azimuth = data.moog_params(AZIMUTH,real_trials,MOOG);
temp_elevation = data.moog_params(ELEVATION,real_trials,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,real_trials,MOOG);
temp_duration = data.moog_params(DURATION,real_trials,MOOG); % in ms
temp_stimType = data.moog_params(STIM_TYPE,real_trials,MOOG);

unique_azimuth = munique(temp_azimuth');
unique_elevation = munique(temp_elevation');
unique_amplitude = munique(temp_amplitude');
unique_duration = munique(temp_duration');
unique_stimType = munique(temp_stimType');

% the least repitition numbers
repNums = (EndTr - BegTr +1)/(length(unique_stimType)*26+1);

% time information
eye_timeWin = 1000/(data.htb_header{EYE_DB}.speed_units/data.htb_header{EYE_DB}.speed/(data.htb_header{EYE_DB}.skip+1)); % in ms
event_timeWin = 1000/(data.htb_header{EVENT_DB}.speed_units/data.htb_header{EVENT_DB}.speed/(data.htb_header{EVENT_DB}.skip+1)); % in ms

event_in_bin = squeeze(data.event_data(:,:,trials));
stimOnT = find(event_in_bin == VSTIM_ON_CD);
stimOffT = find(event_in_bin == VSTIM_OFF_CD);
FPOnT = find(event_in_bin == FP_ON_CD);
fixInT = find(event_in_bin == IN_FIX_WIN_CD);
FPOnBin=ceil(FPOnT*event_timeWin/eye_timeWin);
stimOffBin=ceil(stimOffT*event_timeWin/eye_timeWin);
stimOnBin = ceil(stimOnT*event_timeWin/eye_timeWin);
fixInBin = ceil(fixInT*event_timeWin/eye_timeWin);

% eye data
real_eye_data_LH = squeeze(data.eye_data(1,:,real_trials));
real_eye_data_LV = squeeze(data.eye_data(2,:,real_trials));
real_eye_data_RH = squeeze(data.eye_data(3,:,real_trials));
real_eye_data_RV = squeeze(data.eye_data(4,:,real_trials));

% spon eye data
spon_eye_data_LH = squeeze(data.eye_data(1,:,spon_trials));
spon_eye_data_LV = squeeze(data.eye_data(2,:,spon_trials));
spon_eye_data_RH = squeeze(data.eye_data(3,:,spon_trials));
spon_eye_data_RV = squeeze(data.eye_data(4,:,spon_trials));

% %% calculate data for drawing
% % eye data of real trials (no spon)
% for k = 1:length(unique_stimType)
%     pc = 0;
%     eye_data_LH_bin{k} = [];% for PSTH ANOVA
%     for j = 1:length(unique_elevation)
%         for i = 1:length(unique_azimuth)
%             eye_data_LH_bin{k,j,i} = [];
%             eye_data_LV_bin{k,j,i} = [];
%             eye_data_RH_bin{k,j,i} = [];
%             eye_data_RV_bin{k,j,i} = [];
%             select = find( (temp_azimuth == unique_azimuth(i)) & (temp_elevation == unique_elevation(j)) & (temp_stimType == unique_stimType(k))) ;
%             if sum(select)>0
%                 eye_data_LH_bin{k,j,i} = real_eye_data_LH(eyeBegBin(1,1):eyeEndBin(1,1),select);
%                 eye_data_LV_bin{k,j,i} = real_eye_data_LV(eyeBegBin(1,1):eyeEndBin(1,1),select);
%                 eye_data_RH_bin{k,j,i} = real_eye_data_RH(eyeBegBin(1,1):eyeEndBin(1,1),select);
%                 eye_data_RV_bin{k,j,i} = real_eye_data_RV(eyeBegBin(1,1):eyeEndBin(1,1),select);
%             end
%         end
%     end
% end


% for drawing pictures
eye_data_LH = spon_eye_data_LH;
eye_data_LV = spon_eye_data_LV;
eye_data_RH = spon_eye_data_RH;
eye_data_RV = spon_eye_data_RV;


%% plot figures

% initialize default properties
set(0,'defaultaxesfontsize',24);
colorDefsLBY;
% the lines markers
markers = {
    % markerName % markerTime % marker bin time % color
    'fix_in',fixInT(1),fixInBin(1),colorLGray;
    'stim_on',stimOnT(1),stimOnBin(1),colorDRed;
    'stim_off',stimOffT(1),stimOffBin(1),colorDRed;
    };

figure(111);
clf;
set(gcf,'Position', [0 0 1500 800], 'Name', 'Eye Trace for sponteneous');
h1 = axes('unit','pixels','pos',[100 120 600 600]);
for r = 1:size(spon_eye_data_LH,2)
    color = ['color',num2str(r)];
    plot(h1,eye_data_LH(1:fixInBin(1,1),r),eye_data_LV(1:fixInBin(1,1),r),'linestyle','--','color',eval(color));
    hold on;
    plot(h1,eye_data_LH(fixInBin(1,1):stimOnBin(1,1),r),eye_data_LV(fixInBin(1,1):stimOnBin(1,1),r),'linestyle',':','color',eval(color));
    hold on;
    plot(h1,eye_data_LH(stimOnBin(1,1):stimOffBin(1,1),r),eye_data_LV(stimOnBin(1,1):stimOffBin(1,1),r),'linestyle','-','color',eval(color));
    hold on;
    plot(h1,eye_data_LH(stimOffBin(1,1),r),eye_data_LV(stimOffBin(1,1),r),'linestyle','none','markersize',5,'color','k');
    hold on;
end
set(h1,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5]); % for 2 degree window
set(h1,'xtick',[-1 0 1],'ytick',[-1 0 1]);% for 2 degree window
box off;
title(h1,'Eye trace for sponteneous trials');

% plot h & v as a func of time
h3 = axes('unit','pixels','pos',[800 120 600 200]);
for r = 1:size(spon_eye_data_LH,2)
    color = ['color',num2str(r)];
    plot(h3,1:fixInBin(1,1),eye_data_LV(1:fixInBin(1,1),r),'linestyle','--','color',eval(color));
    hold on;
    plot(h3,fixInBin(1,1):stimOnBin(1,1),eye_data_LV(fixInBin(1,1):stimOnBin(1,1),r),'linestyle',':','color',eval(color));
    hold on;
    plot(h3,stimOnBin(1,1):stimOffBin(1,1),eye_data_LV(stimOnBin(1,1):stimOffBin(1,1),r),'linestyle','-','color',eval(color));
    hold on;
end
box off;
title(h3,'Vetical');
set(h3,'xlim',[1 stimOffBin(1)+10],'ylim',[-1.5 1.5]); % for 2 degree window
set(h3,'ytick',[-1 0 1]);% for 2 degree window
set(h3,'xtick',[]);

h2 = axes('unit','pixels','pos',[800 440 600 200]);
for r = 1:size(spon_eye_data_LH,2)
    color = ['color',num2str(r)];
    plot(h2,1:fixInBin(1,1),eye_data_LH(1:fixInBin(1,1),r),'linestyle','--','color',eval(color));
    hold on;
    plot(h2,fixInBin(1,1):stimOnBin(1,1),eye_data_LH(fixInBin(1,1):stimOnBin(1,1),r),'linestyle',':','color',eval(color));
    hold on;
    plot(h2,stimOnBin(1,1):stimOffBin(1,1),eye_data_LH(stimOnBin(1,1):stimOffBin(1,1),r),'linestyle','-','color',eval(color));
    hold on;
end
box off;
title(h2,'Horizontal');
set(h2,'xlim',[1 stimOffBin(1)+10],'ylim',[-1.5 1.5]); % for 2 degree window
set(h2,'ytick',[-1 0 1]);% for 2 degree window
set(h2,'xtick',[]);

for n = 1:size(markers,1)
    plot(h2,[markers{n,3} markers{n,3}], [-1.5,1.5], ':','color',markers{n,4},'linewidth',1.5);
    hold on;
    box off;
    plot(h3,[markers{n,3} markers{n,3}], [-1.5,1.5], ':','color',markers{n,4},'linewidth',1.5);
    hold on;
    box off;
end

axes('unit','pixels','pos',[800 60 600 40]);
set(gca,'xlim',[1 stimOffBin(1)+10],'ylim',[0 10]); % to align with h2 & h3
text(fixInBin(1)-30,10,'fix in','fontsize',16);
text(stimOnBin(1)-30,2,'stim on','fontsize',16);
text(stimOffBin(1)-30,2,'stim off','fontsize',16);
axis off;

% text on the figure
axes('unit','pixels','pos',[800 745 700 80]);
xlim([0,100]);
ylim([0,10]);
FileName_Temp = num2str(FILE);
FileName_Temp =  FileName_Temp(1:end-4);
str1 = [FileName_Temp,'\_Ch' num2str(SpikeChan)];
text(0,2,str1,'fontsize',20);
axis off;

% save the figure
str2 = [FileName_Temp,'_Ch' num2str(SpikeChan), '_eyetrace_spon'];
%  str3 = [str1, '_eyetrace'];
saveas(gcf,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_eyetrace\' str2], 'emf')
%  saveas(3,['Z:\LBY\Recording Data\Qiaoqiao\3D_Tuning_eyetrace\' str3], 'emf');


% % %% Data Saving
% %
% % % Reorganized. HH20141124
% % config.batch_flag = batch_flag;
% %
% % % Output information for test. HH20160415
% % if isempty(batch_flag)
% %     config.batch_flag = 'test.m';
% %     disp('Saving results to \batch\test\ ');
% % end
% %
% % %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % if exist('CP','var')  % Full version
% %     result = PackResult(FILE, SpikeChan, repetitionN, unique_stim_type, ... % Obligatory!!
% %                     stim_type_per_trial, heading_per_trial, choice_per_trial,... % Trial info
% %                     align_markers, align_offsets_others, sort_info, smoothFactor ,...
% %                     PREF, PREF_CP_obsolete, PREF_target_location, outcome_mask_enable, ...
% %                     binSize_rate, stepSize_rate, rate_ts, binSize_CP, stepSize_CP, CP_ts ,...
% %                     spike_aligned, spike_hist, CP, PSTH,...
% %                     ChoiceDivergence_ALL, ChoiceDivergence_Difficult, ChoiceDivergence_Easy,ChoicePreference,ChoicePreference_pvalue,...
% %                     ModalityDivergence, ModalityPreference, ModalityPreference_pvalue);
% %     % Figures to save
% %     config.save_figures = [60 + (1:length(sort_info)) ,161, 1899, 1999, 2000, 59];
% %
% % else
% %     result = PackResult(FILE, SpikeChan, repetitionN, unique_stim_type, ... % Obligatory!!
% %                         stim_type_per_trial, heading_per_trial, choice_per_trial,... % Trial info
% %                         align_markers, align_offsets_others, sort_info, smoothFactor ,...
% %                         PREF, PREF_target_location, outcome_mask_enable, ...
% %                         binSize_rate, stepSize_rate, rate_ts,...
% %                         spike_aligned, spike_hist, PSTH,...
% %                         ChoiceDivergence_ALL, ChoiceDivergence_Difficult, ChoiceDivergence_Easy,ChoicePreference,ChoicePreference_pvalue,...
% %                         ModalityDivergence, ModalityPreference, ModalityPreference_pvalue);
% %
% % %     result = PackResult(FILE, SpikeChan, repetitionN, unique_stim_type, ... % Obligatory!!
% % %                         PREF_target_location);
% %
% %     % Figures to save
% %     config.save_figures = [];
% %
% % end
% %
% % config.suffix = 'PSTH';
% % config.xls_column_begin = 'HD_rep';
% % % config.xls_column_end = 'HD_comb_p';
% % config.xls_column_end = 'HD_comb_ChoicePref_p';
% %
% % % Only once
% % config.sprint_once_marker = 'gs';
% % config.sprint_once_contents = 'result.repetitionN,num2str(result.PSTH{2,1,1}.ts)';
% % % Loop across each stim_type
% % % config.sprint_loop_marker = {'gg';
% % %                            'gg';
% % %                            'sss'};
% % % config.sprint_loop_contents = {'result.CP{2,k}.Psy_para(2), result.CP{2,k}.Psy_para(1)';
% % %                         'result.CP{2,k}.CP_grand_center, result.CP{2,k}.CP_grand_sac';
% % %                         'num2str(result.PSTH{2,1,1}.ys((k-1)*2+1,:)), num2str(result.PSTH{2,1,1}.ys((k-1)*2+2,:)), num2str(result.PSTH{2,1,1}.ps(k,:))'};
% %
% % % Replace the meaningless CPs with ChoicePreference (surprisingly, they have the same abbreviation!). HH20160419
% % config.sprint_loop_marker = {'gg';
% %                            'gg';
% %                            'gg'};
% % config.sprint_loop_contents = {'result.CP{2,k}.Psy_para(2), result.CP{2,k}.Psy_para(1)';
% %                         'result.CP{2,k}.CP_grand_center, result.CP{2,k}.CP_grand_sac';
% %                         'result.ChoicePreference(1,k), result.ChoicePreference_pvalue(1,k)'};
% %
% % config.append = 1; % Overwrite or append
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % SaveResult(config, result);
% %
return;