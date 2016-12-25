%%  LIP HD Pooled Data
% Begin @ HH 201406

function function_handles = Group_HD_dt(XlsData)

% try
tmp1 = []; tmp2 = [];

%% Get data
num = XlsData.num;
txt = XlsData.txt;
raw = XlsData.raw;
header = XlsData.header;

stim_type_num = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch address
mat_address = {
    % Major protocol goes here (Address, Suffix)
     'Z:\Data\Tempo\Batch\20160424_HD_dt_m10','PSTH';
%      'Z:\Data\Tempo\Batch\20160603_HD_dt_m10_Smooth20ms','PSTH';
%      'Z:\Data\Tempo\Batch\20150418_LIP_HD_m5_m10_modifiedBatch','PSTH';
    %  'Z:\Data\Tempo\Batch\20150411_LIP_HD_m5','PSTH';
    %  'Z:\Data\Tempo\Batch\20141118_LIP_Decision_WithAUC_ResultPSTHOrderChanged','PSTH';
    
    % Associative protocols
     'Z:\Data\Tempo\Batch\20150725_BP_allAreas_m5_m10','MemSac';
    % 'Z:\Data\Tempo\Batch\20150411_LIP_memSac_m5_m10','MemSac';
    };

% Global Mask for the each protocol (for XLS)
mask_all = {
    strcmp(txt(:,header.Protocol),'HD_dt') & (strcmpi(txt(:,header.Area),'LIP') | strcmpi(txt(:,header.Area),'LIP-V')) ; %& ( num(:,header.Monkey)~=10 | (num(:,header.Session) >= 27 & num(:,header.Session) <= 42));
    % &(num(:,header.HD_rep) >= 8); % This has been moved to sections below (I will find repN from .mat files)
    strcmp(txt(:,header.Protocol),'MemSac') & (strcmpi(txt(:,header.Area),'LIP') | strcmpi(txt(:,header.Area),'LIP-V'));
    }; % Now no constraint on monkeys

% Add flexible monkey mask here (but I've still decided to choose monkey for analysis below). HH20150723
monkey_included_for_loading = [5 10]; 

monkey_mask_for_loading = false(size(num,1),1);
for mm = 1:length(monkey_included_for_loading)
    monkey_mask_for_loading = monkey_mask_for_loading | (num(:,header.Monkey) == monkey_included_for_loading(mm));
end

% Now apply monkey mask
for mm = 1:length(mask_all)
    mask_all{mm} = mask_all{mm} & monkey_mask_for_loading;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for pp = 1:size(mat_address,1)
    % Xls Data
    xls_num{pp} = num(mask_all{pp},:);
    xls_txt{pp} = txt(mask_all{pp},:);
    xls_raw{pp} = raw(mask_all{pp},:);
    
    % Basic information : E.g.  '20140721m05s034h2x04y12d08478u6m5c67r2'
    cell_info_tmp = xls_num{pp}(:,[header.Date:header.Yloc header.Depth header.Chan1]);
    
    % Override MU with SU. HH20150422
    cell_info_tmp(xls_num{pp}(:,header.Units_RealSU)==1 & xls_num{pp}(:,header.Chan1)==1,end) = 5;
    
    cell_info{pp} = strsplit(sprintf('%8dm%02ds%03dh%01dx%02dy%02dd%05du%02d\n',cell_info_tmp'),'\n')';
    cell_info{pp} = strcat(cell_info{pp}(1:end-1),xls_txt{pp}(:,header.FileNo));
end
cd(mat_address{1,1});

%% Establish mat-to-xls relationship and load data into group_result.mat
% %{

global group_result; % This lets me reload the whole m.file without loading group_result, which speeds up my debugging.

%

if isempty(group_result)
    group_result(size(xls_txt{1},1)).cellID = [];  % The unique cellID
    load_group_result();
end

    function load_group_result()
        %
        % Load .mat files and put the data into a large structure array "group_result(i).mat_raw"
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tolerance_enable = 0;
        depthTol = 200; % Depth tolerance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tic;
        progressbar('Load .mat files');
        not_match = 0;
        
        for major_i = 1:length(group_result) % Major protocol loop
            for pp = 1:size(mat_address,1)
                if pp == 1 % Major protocol
                    % Get basic information for cellID
                    file_i = major_i;
                else % Search for corresponding files in xls
                    match_i = find(strncmp(cell_info{1}{major_i},cell_info{pp},32));
                    
                    if length(match_i) == 1 % Exact match
                        file_i = match_i;
                    elseif length(match_i) > 1
                        fprintf('More than one match have been found for %s,\n Major ID = %s\n', xls_txt{1}{major_i,header.FileNo},cell_info{1}{major_i});
                        disp(cell_info{pp}(match_i));
                        %                 file_choose = input('Which one do you want?');
                        %                 file_i = match_i(file_choose);
                        file_i = match_i(1);  % The first appearance by default.
                        
                        %                 keyboard;
                    else % Special cases for inexact match (MU-SU)
                        if ~tolerance_enable
                            fprintf('No exact matched files for %s, ID = %s\n', xls_txt{1}{major_i,header.FileNo},cell_info{1}{major_i});
                            not_match = not_match + 1;
                            file_i = NaN;
                        else
                            % Unit tolerance
                            match_i = find(strncmp(cell_info{1}{major_i},cell_info{pp},29),1);
                            
                            if length(match_i) == 1
                                fprintf('Unit tolerance found:\n   %s (major)\n   %s\n\n', cell_info{1}{major_i}, cell_info{pp}{match_i});
                                file_i = match_i;
                            else
                                % Depth tolerance
                                match_i = find(strncmp(cell_info{1}{major_i},cell_info{pp},23),1);
                                depDiff = abs(xls_num{pp}(match_i,header.Depth)-xls_num{1}(major_i,header.Depth));
                                
                                if any(depDiff <= depthTol) % Within depth tolerance
                                    match_i = match_i(find(depDiff == min(depDiff),1));
                                    fprintf('Depth tolerance found:\n   %s (major)\n   %s\n\n', cell_info{1}{major_i}, cell_info{pp}{match_i});
                                    file_i = match_i;
                                else
                                    fprintf('!! All tolerance failed:\n    %s (major)\n\n',cell_info{1}{major_i});
                                    file_i = NaN;
                                    not_match = not_match + 1;
                                end
                            end
                        end
                        
                    end
                end
                
                % Load .mat for major and associative protocols
                if ~isnan(file_i)
                    try
                        group_result(major_i).cellID{pp} = cell_info{pp}{file_i};
                        
                        mat_file_name = sprintf('%s_%g',xls_txt{pp}{file_i,header.FileNo},xls_num{pp}(file_i,header.Chan1));
                        mat_file_fullname = [mat_address{pp,1} '\' mat_file_name '_' mat_address{pp,2}];
                        
                        raw = load(mat_file_fullname);
                        group_result(major_i).(['mat_raw_' mat_address{pp,2}]) = raw.result;  % Dynamic structure
                    catch
                        fprintf('Error Loading %s\n',[mat_file_name '_' mat_address{pp,2}]);
                        keyboard;
                    end
                else
                        group_result(major_i).(['mat_raw_' mat_address{pp,2}]) = [];  
                end
            end
            
            progressbar(major_i/length(group_result));
            
        end
        fprintf('Loaded %g files (%g of them are incomplete).\n',length(group_result),not_match);
        toc
        
    end

%% Get some values from group_result(i).mat_raw_xxx our to group_result(i) for easier access
%  Preprocess each .mat_raw for different purposes (which have not been done in Batch Processing)

% HH20150414: Modality divergence added.
% HH20150419: Choice Divergence/Preference and Modality Divergence/Preference have been moved to Batch processing
% HH20150419: Different temporal alignment have been included

%%%%%%% Order corresponding to "Sort_Id" in TEMPO_GUI processing %%%%
ALL_CHOICE = 1; CORRECT_ANGLE = 2; CHOICE_DIFFICULT = 3; OUTCOME = 4; WRONG_ANGLE = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(group_result)
    
    group_result(i).repN = group_result(i).mat_raw_PSTH.repetitionN;
    
    unique_stim_type = group_result(i).mat_raw_PSTH.unique_stim_type;  % In case we don't have all three conditions
    group_result(i).length_unique_stim_type = length(unique_stim_type);
    group_result(i).unique_heading = group_result(i).mat_raw_PSTH.CP{1,1}.raw_CP_result{1}.Neu_tuning(:,1);
    
    %     if isempty(group_result(i).mat_raw_PSTH)
    %         continue
    %     end
    
    % Time markers (stim on, stim off, sac on)
    group_result(i).time_marker{1} = [0 mean(group_result(i).mat_raw_PSTH.align_offsets_others{1})];
    group_result(i).time_marker{2} = [mean(group_result(i).mat_raw_PSTH.align_offsets_others{2}) 0];
    
    % Psychophysics
    group_result(i).Psy_para = nan(2,3);
    for stim_type = 1:stim_type_num
        k = find(stim_type == unique_stim_type);
        if ~isempty(k)   % We have this condition
            group_result(i).Psy_para(1,stim_type) = group_result(i).mat_raw_PSTH.CP{1,k}.Psy_para(2); % Threshold
            group_result(i).Psy_para(2,stim_type) = group_result(i).mat_raw_PSTH.CP{1,k}.Psy_para(1); % Bias
        end
    end
    
    for j = 1:2 % For two temporal alignments
        
        % 1). For CP, here I decide to use pref. direction defined in Anne 2014 paper (100-200
        % ms before decision). But in the previous CP_HH.m, I assigned pref. direction for each time bin.
        % So here I check whether the local pref. direction is aligned with the "mode" (Zhong Shu) preferred
        % direction several time bins before decision. If not, I flip the original CP value.
        %  XXX @HH20141203 NO NEED To FLIP ANYMORE (I did this in HD_PSTH_HH)
        
        CP_ts = group_result(i).mat_raw_PSTH.CP{j,1}.ts;
        group_result(i).CP_ts{j} = CP_ts;
        
        
        %         group_result(i).CP = NaN(3,length(CP_ts));
        %         group_result(i).CP_p_value = NaN(3,length(CP_ts));
        
        % *** NOTE: I have done "always output three conditions" in TEMPO_GUI for all data except CP and PSTH. @HH20150419
        group_result(i).CP{j} = nan(3,length(CP_ts));
        group_result(i).CP_p_value{j} = nan(3,length(CP_ts));
        
        for stim_type = 1:stim_type_num  % Always output three conditions
            
            k = find(stim_type == unique_stim_type);
            
            if ~isempty(k)   % We have this condition
                %             for tt = 1:length(CP_ts)
                %                 pref(tt) = group_result(i).mat_raw_PSTH.CP{j,k}.raw_CP_result{tt}.pref;
                %             end
                
                %             pref_Anne = mode(pref((CP_ts > -500) & (CP_ts <= -200)));  % Each neuron is assigned by only one pref. direction
                %             need_flip = (pref ~= pref_Anne);
                
                group_result(i).CP{j}(stim_type,:) = group_result(i).mat_raw_PSTH.CP{j,k}.CP_grand;
                
                % --------  No need to flip. HH20141203 -------
                %             group_result(i).CP(stim_type,need_flip) = 1 - group_result(i).CP(k,need_flip);
                
                % CP p value (permutation 1000)
                group_result(i).CP_p_value{j}(stim_type,:) = group_result(i).mat_raw_PSTH.CP{j,k}.CP_p;
            end
        end
        
        % 2). Rate p_value is straightforward
        
        rate_ts_temp = group_result(i).mat_raw_PSTH.PSTH{j,ALL_CHOICE,1}.ts;
        group_result(i).rate_ts{j} =  rate_ts_temp;
        
        group_result(i).rate_p_value{j} = NaN(3,length(rate_ts_temp));
        
        for stim_type = 1:stim_type_num  % Always output three conditions
            
            k = find(stim_type == unique_stim_type);
            if ~isempty(k)   % We have this condition
                group_result(i).rate_p_value{j}(stim_type,:) =  group_result(i).mat_raw_PSTH.PSTH{j,ALL_CHOICE,1}.ps(k,:);
            end
            
        end
    end
    
    %     figure(1799); clf;
    %     set(0,'defaultAxesColorOrder',[0 0 1; 1 0 0; 0 0.8 0;]);
    %     plot(rate_ts, sqrt(-log(result{i}.rate_p_value'))); SetFigure(); axis tight;
    %     title('Sqrt(-log(Rate p value))');
    
    group_result(i).PREF_PSTH = group_result(i).mat_raw_PSTH.PREF;
    
    % 3). Next I calculate the "Choice divergence": AUC between different
    % choices under different conditions and difficulty levels.
    
    group_result(i).ChoicePreference = group_result(i).mat_raw_PSTH.ChoicePreference;
    
    % In TEMPO_GUI, choice preference uses the cell's PREF as its preferred direction
    % (because the hemisphere is unknown unless accessible to Result.xls)
    % Now I transform it to be related to "Contralateral" (Anne 2014)
    group_result(i).if_contralateral = xls_num{1}(i,header.Hemisphere) ~= group_result(i).PREF_PSTH;
    group_result(i).ChoicePreference = group_result(i).ChoicePreference * sign(group_result(i).if_contralateral - 0.5);
    
    group_result(i).ChoicePreference_pvalue =  group_result(i).mat_raw_PSTH.ChoicePreference_pvalue;
    
    %{
    % XXX @ HH20150419 Moved to Batch processing
    
%     group_result(i).PREF_ANNE_CP = group_result(i).mat_raw_PSTH.PREF_CP;
    
    group_result(i).ChoiceDivergence_ALL = NaN(3,length(rate_ts));
    group_result(i).ChoiceDivergence_Difficult = NaN(3,length(rate_ts));
    group_result(i).ChoiceDivergence_Easy = NaN(3,length(rate_ts));
    
    for stim_type = 1:3  % Always output three conditions
        
        k = find(stim_type == unique_stim_type);
        
        if ~isempty(k)   % We have this condition
            
            % Re-determine the preferred direction according to Anne's method (see above)
            mean_rates = group_result(i).mat_raw_PSTH.PSTH{j,ALL_CHOICE,1}.ys(2*k-1:2*k,:);      % Note for sortInd = 1 (ALL_CHOICE), all conditions are put together
            first_larger = mean(mean_rates(1, rate_ts > -500 & rate_ts <= -200)) > mean(mean_rates(2,rate_ts > -500 & rate_ts <= -200));
            
            % Calculate auROC for each time bin
            for tt = 1:length(rate_ts)
                group_result(i).ChoiceDivergence_ALL(stim_type,tt) = rocN(group_result(i).mat_raw_PSTH.PSTH{j,ALL_CHOICE,1}.raw{2*k-1}(:,tt),...
                    group_result(i).mat_raw_PSTH.PSTH{j,ALL_CHOICE,1}.raw{2*k}(:,tt)) - 0.5;      % Default: first larger
                
                group_result(i).ChoiceDivergence_Difficult(stim_type,tt) = rocN(group_result(i).mat_raw_PSTH.PSTH{j,CHOICE_DIFFICULT,k}.raw{1}(:,tt),...
                    group_result(i).mat_raw_PSTH.PSTH{j,CHOICE_DIFFICULT,k}.raw{2}(:,tt)) - 0.5;
                
                group_result(i).ChoiceDivergence_Easy(stim_type,tt) = rocN(group_result(i).mat_raw_PSTH.PSTH{j,CHOICE_DIFFICULT,k}.raw{3}(:,tt),...
                    group_result(i).mat_raw_PSTH.PSTH{j,CHOICE_DIFFICULT,k}.raw{4}(:,tt)) - 0.5;
            end
            
            if ~first_larger % If not, flip
                group_result(i).ChoiceDivergence_ALL(stim_type,:) = - group_result(i).ChoiceDivergence_ALL(k,:);
                group_result(i).ChoiceDivergence_Difficult(stim_type,:) = - group_result(i).ChoiceDivergence_Difficult(k,:);
                group_result(i).ChoiceDivergence_Easy(stim_type,:) = - group_result(i).ChoiceDivergence_Easy(k,:);
            end
            
        end  % if ~isempty(k)
        
    end
    
    %}
    
    %     figure(1899); clf;
    %     plot(rate_ts, group_result(i).ChoiceDivergence_ALL','Linewidth',2); SetFigure(); axis tight;
    %     title([mat_file_name{i} '_ Choice Divergence']);
    %     print(1899,'-dbitmap',[mat_file_fullname{i} '_ChoiceDivergenceAll.bmp']);
    %
    %     figure(1999); clf;
    %     plot(rate_ts, group_result(i).ChoiceDivergence_Easy'- group_result(i).ChoiceDivergence_Difficult','Linewidth',2); SetFigure(); axis tight;
    %     title([mat_file_name{i} '_ Easy - Difficult']);
    %     print(1999,'-dbitmap',[mat_file_fullname{i} '_ChoiceDivergenceEasyMinusDifficult.bmp']);
    
    
    % 4). Here comes the "modality divergence" (1-2,1-3,2-3)
    %     AUC between different modalities (regardless of choice)     HH20150415
    
    group_result(i).ModalityPreference = group_result(i).mat_raw_PSTH.ModalityPreference;
    group_result(i).ModalityPreference_pvalue =  group_result(i).mat_raw_PSTH.ModalityPreference_pvalue;
    
    %{
    % XXX @ HH20150419 Moved to Batch processing
  
    group_result(i).ModalityDivergence = NaN(3,length(rate_ts));
    
    modality_pair = {[1 2],[1 3],[2 3]};  % The later is set to be "Pref modality"
    
    if length(group_result(i).mat_raw_PSTH.unique_stim_type) == 3 % If we have all three modalities
        
        % Reorganize the modality data (which has not been done in batch processing before)
        
        raw_for_md = group_result(i).mat_raw_PSTH.PSTH{j,ALL_CHOICE,1}.raw;
        
        for k = 1:3
            modality_PSTH{k} = [ raw_for_md{(k-1)*2+2} ; raw_for_md{k*2} ];
        end
        
        % Now calculate auROC
        
        for mp = 1:length(modality_pair)
            
            pref_mod = modality_pair{mp}(2);
            null_mod = modality_pair{mp}(1);
            
            % Calculate auROC for each time bin
            for tt = 1:length(rate_ts)
                group_result(i).ModalityDivergence (mp,tt) = ...
                    rocN(modality_PSTH{pref_mod}(:,tt), modality_PSTH{null_mod}(:,tt)) - 0.5;      % Default: first larger
            end
            
            
        end
    end
    
    %}
    
    
    % 5) Mem-sac dynamics (from .mat file) goes here
    if ~isempty(group_result(i).mat_raw_MemSac)
        group_result(i).MemSac_p = group_result(i).mat_raw_MemSac.p;
        group_result(i).MemSac_DDI = group_result(i).mat_raw_MemSac.DDI;
        group_result(i).MemSac_vectSum = group_result(i).mat_raw_MemSac.vectSum;
        group_result(i).MemSac_AI = group_result(i).mat_raw_MemSac.activityIndex;
        group_result(i).MemSac_ts = group_result(i).mat_raw_MemSac.t_centers{3};
        group_result(i).MemSac_PSTH = reshape([group_result(i).mat_raw_MemSac.result_PSTH_anne_mean{3,:}],length(group_result(i).MemSac_ts),[]);
    end
    
    %     catch
    %         fprintf('Error in preprocessing %g\n',i);
    %     end
    
end

% ------ Load Gaussian velocity -----
% Measured by accelerometer at 109. Should retest it on 103. HH20150422
% I've done that, they are almost the same.
temp = load('Gaussian_vel_real_sigma3p5.mat');
Gauss_vel = temp.Gaussian_vel_real_sigma3p5;
Gauss_vel(:,2) = Gauss_vel(:,2)/max(Gauss_vel(:,2));

%{
% XXX @HH20150419.
% Saving and Loading group_result.mat seem too slow, so I choose to
% recalculate it everytime.

disp('Saving to group_result...'); tic
save('group_result.mat','group_result','-v7.3');
disp('Done.'); toc

%}

%}

% keyboard

%% Reorganize data into matrices which are easier to plot and compare
%{
if ~exist('group_result')
    disp('Loading group_result ...'); tic
    load group_result;
    disp('Done.');
    toc
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smooth_factor_for_divergence = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(group_result);
rate_ts = group_result(end).rate_ts;
CP_ts = group_result(end).CP_ts;
memsac_ts = []; %group_result(end).MemSac_ts;

% Initialization

% j-sensitive
for j = 1:2
    % For normalized PSTH plotting (weighted by dynamic ranges for each neuron)
    PSTH_all_Norm{j} = NaN(N,length(rate_ts{j}),stim_type_num*2);
    
%     PSTH_correct_angles_Norm{j} = NaN(N,length(rate_ts{j}),1+length(group_result(end).unique_heading),stim_type_num);  % Two zero headings
    PSTH_correct_angles_Norm{j} = NaN(N,length(rate_ts{j}),length(group_result(end).unique_heading),stim_type_num);  % HH20160426

    PSTH_outcomes_Norm{j} = NaN(N,length(rate_ts{j}),4,stim_type_num);

%     PSTH_wrong_angles_Norm{j} = NaN(N,length(rate_ts{j}),length(group_result(end).unique_heading)-1,stim_type_num); % @HH20150523
    PSTH_wrong_angles_Norm{j} = NaN(N,length(rate_ts{j}),length(group_result(end).unique_heading),stim_type_num); % % HH20160426
    
    % Also pack raw data for different ways of weighted sum PSTH plotting (weighted by SVM / targeted dimensionality reduction, etc.)
    PSTH_all_raw{j} = NaN(N,length(rate_ts{j}),stim_type_num*2);
%     PSTH_correct_angles_raw{j} = NaN(N,length(rate_ts{j}),1+length(group_result(end).unique_heading),stim_type_num);
    PSTH_correct_angles_raw{j} = NaN(N,length(rate_ts{j}),length(group_result(end).unique_heading),stim_type_num);% HH20160426
    PSTH_outcomes_raw{j} = NaN(N,length(rate_ts{j}),4,stim_type_num);
%     PSTH_wrong_angles_raw{j} = NaN(N,length(rate_ts{j}),length(group_result(end).unique_heading)-1,stim_type_num); % @HH20150523
    PSTH_wrong_angles_raw{j} = NaN(N,length(rate_ts{j}),length(group_result(end).unique_heading),stim_type_num); % HH20160426
    
    CP{j} = NaN(N,length(CP_ts{j}),3);
    
    ChoiceDiv_All{j} = NaN(N,length(rate_ts{j}),3);
    ChoiceDiv_Easy{j} = NaN(N,length(rate_ts{j}),3);
    ChoiceDiv_Difficult{j} = NaN(N,length(rate_ts{j}),3);
    ChoiceDiv_EasyMinusDifficult{j} = NaN(N,length(rate_ts{j}),3);
    
    ChoiceDiv_ModDiffer{j} = NaN(N,length(rate_ts{j}),3);  % 3-1, 3-2, 1-2
    ModDiv_All{j} = NaN(N,length(rate_ts{j}),3); % HH20140415   % 2-1, 3-1, 3-2
end

% j-insensitive
MemSac_DDI = NaN(N,6);
MemSac_PREF_Null_DI = NaN(N,6);
MemSac_actual_DI = NaN(N,1); % DI of actuall target locations. @HH20150524

MemSac_PREFmNULL_PSTH = []; %NaN(N,length(group_result(2).MemSac_ts));

MemSac_PSTH_AngDiff = NaN(N,6);

% Reorganize data
for i = 1:N
    
    % j-sensitive
    for j = 1:2
        
        % Calculate normalized PSTH (from stim onset to sac onset)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        norm_range =  group_result(i).time_marker{j}(1) <= rate_ts{j} & rate_ts{j} <= group_result(i).time_marker{j}(3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Normalization according to dynamic range of each neuron
        % (In the norm_range, lowest in all conditions = 0, highest = 1)
        PSTH_all_Norm_this = group_result(i).mat_raw_PSTH.PSTH{j,ALL_CHOICE,1}.ys;
        
        PSTH_all_raw_this = PSTH_all_Norm_this;
        
        offset = min(min(PSTH_all_Norm_this(:,norm_range)));
        gain = max(max(PSTH_all_Norm_this(:,norm_range))) - offset;   % Dynamic range (e.g., Meister 2013 JNS)
        %         gain = mean(mean(PSTH_all_Norm_this(:,-300 < rate_ts{1} & rate_ts{1} < -100)));  % Background activity before stimulus onset
        
        PSTH_all_Norm_this = PSTH_all_Norm_this - offset;
        PSTH_all_Norm_this = PSTH_all_Norm_this / gain;
        
        % Save weights of classical normalization method for comparison
        % with weights from other methods (targed dim reduction, etc.)
        weights_normalized_PSTH(i) = 1/gain;
        
        for stim_type = 1:stim_type_num % Stim_type check
            
            % I have done "always output three conditions" in TEMPO_GUI for all data except CP and PSTH. @HH20150419
            k = find(stim_type == group_result(i).mat_raw_PSTH.unique_stim_type);
            if ~isempty(k)   % We have this condition
                
                % ----------- Pack PSTH_all_Norm -------------
                PSTH_all_Norm{j}(i,:,k*2-1) = PSTH_all_Norm_this(k*2-1,:);
                PSTH_all_Norm{j}(i,:,k*2) = PSTH_all_Norm_this(k*2,:);
                PSTH_all_raw{j}(i,:,k*2-1) = PSTH_all_raw_this(k*2-1,:);
                PSTH_all_raw{j}(i,:,k*2) = PSTH_all_raw_this(k*2,:);
                
                % ---------- Normalize and pack PSTH_angles_Norm ---------
                PSTH_correct_angles_norm_this = group_result(i).mat_raw_PSTH.PSTH{j,CORRECT_ANGLE,k}.ys;
                PSTH_correct_angles_raw_this = PSTH_correct_angles_norm_this;
                PSTH_correct_angles_norm_this = PSTH_correct_angles_norm_this - offset;
                PSTH_correct_angles_norm_this = PSTH_correct_angles_norm_this / gain;
                
                if size(PSTH_correct_angles_norm_this,1) == size(PSTH_correct_angles_Norm{j},3)
                    PSTH_correct_angles_Norm{j}(i,:,:,k) = PSTH_correct_angles_norm_this';
                    PSTH_correct_angles_raw{j}(i,:,:,k) = PSTH_correct_angles_raw_this';
                elseif size(PSTH_correct_angles_norm_this,1) == size(PSTH_correct_angles_Norm{j},3) - 2 % Without zero heading
                    PSTH_correct_angles_Norm{j}(i,:,3:end,k) = PSTH_correct_angles_norm_this';
                    PSTH_correct_angles_raw{j}(i,:,3:end,k) = PSTH_correct_angles_raw_this';
                else
                    disp('No match PSTH_angles_norm...');
                end
                
                % --------- Normalize and pack PSTH_outcome_Norm ---------
                % {
                %    'outcome_per_trial',{'=='}, [CORRECT ERR_WRONG_CHOICE], .................
                %    'choice_per_trial', {'=='}, [PREF NULL], [], {'-'; '--'}, ..............
                % }
                
                PSTH_outcome_norm_this = group_result(i).mat_raw_PSTH.PSTH{j,OUTCOME,k}.ys;
                PSTH_outcome_raw_this = PSTH_outcome_norm_this;
                
                PSTH_outcome_norm_this = PSTH_outcome_norm_this - offset;
                PSTH_outcome_norm_this = PSTH_outcome_norm_this / gain;
                
                PSTH_outcomes_Norm{j}(i,:,:,k) = PSTH_outcome_norm_this';
                PSTH_outcomes_raw{j}(i,:,:,k) = PSTH_outcome_raw_this';
                
                
                % -------- Normalize and pack PSTH_wrong_angles_Norm ------
                % Note in batch file (HeadingDis_cum_PSTH_HH.m), sort_ind{5} has been defined as :
                %     {
                %        'outcome_per_trial',{'=='}, [ERR_WRONG_CHOICE], ...........
                %        'heading_per_trial''',{'=='}, unique_heading(unique_heading~=0), ................
                %     }
                % That is to say, unlike the others, PSTH_wrong_angles_Norm does not have any information about PREF and
                % NULL, but just headings. Therefore, I need to check the PREF of each cell in order to organize
                % PSTH_wrong_angles like this:
                %      [smallest heading & PREF, smallest heading & NULL, ... , largest heading & PREF, largest heading & NULL]
                % which is similar to PSTH_correct_angles except that it does not have 0 headings.
                
                % Maybe someday I will move this part to the batch file.  @HH20150524
                
                PSTH_wrong_angles_norm_this = group_result(i).mat_raw_PSTH.PSTH{j,WRONG_ANGLE,k}.ys;
                order = nan(size(group_result(i).mat_raw_PSTH.sort_info{WRONG_ANGLE}{2}{2,3}));
                
                if group_result(i).mat_raw_PSTH.PREF == 1 % If PREF is leftward, then the positive headings go first (because they are WRONG trials)
                    order(1:2:end) = length(order)/2+1:1:length(order);
                    order(2:2:end) = length(order)/2:-1:1;
                else % vice versa
                    order(1:2:end) = length(order)/2:-1:1;
                    order(2:2:end) = length(order)/2+1:1:length(order);
                end
                
                PSTH_wrong_angles_norm_this = PSTH_wrong_angles_norm_this(order,:);
                PSTH_wrong_angles_raw_this = PSTH_wrong_angles_norm_this;
                
                PSTH_wrong_angles_norm_this = PSTH_wrong_angles_norm_this - offset;
                PSTH_wrong_angles_norm_this = PSTH_wrong_angles_norm_this / gain;
                
                PSTH_wrong_angles_Norm{j}(i,:,:,k) = PSTH_wrong_angles_norm_this';
                PSTH_wrong_angles_raw{j}(i,:,:,k) = PSTH_wrong_angles_raw_this';
                
                
            end
            
            % Stim type check already done
            ChoiceDiv_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_ALL{j}(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_Easy{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Easy{j}(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_Difficult{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Difficult{j}(stim_type,:),smooth_factor_for_divergence);
            
            ChoiceDiv_EasyMinusDifficult{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Easy{j}(stim_type,:) ...
                - group_result(i).mat_raw_PSTH.ChoiceDivergence_Difficult{j}(stim_type,:),smooth_factor_for_divergence);
            
            
            CP{j}(i,:,stim_type) = group_result(i).CP{j}(stim_type,:);
            CP_p{j}(i,:,stim_type) = group_result(i).CP_p_value{j}(stim_type,:);
            
            if stim_type <=3
                ModDiv_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ModalityDivergence{j}(stim_type,:),smooth_factor_for_divergence);
            end
        end
        
        ChoiceDiv_ModDiffer{j}(i,:,1) = smooth(ChoiceDiv_All{j}(i,:,3) - ChoiceDiv_All{j}(i,:,1),smooth_factor_for_divergence);
        ChoiceDiv_ModDiffer{j}(i,:,2) = smooth(ChoiceDiv_All{j}(i,:,3) - ChoiceDiv_All{j}(i,:,2),smooth_factor_for_divergence);
        ChoiceDiv_ModDiffer{j}(i,:,3) = smooth(ChoiceDiv_All{j}(i,:,1) - ChoiceDiv_All{j}(i,:,2),smooth_factor_for_divergence);
    end
    
    % j-insensitive
    
    % Mem-sac stuff (Note some are NaNs)
    if isfield(group_result(i),'MemSac_DDI') & ~isempty(group_result(i).MemSac_DDI)
        % -------- Global memsac indicator -------
        MemSac_DDI(i,:) = group_result(i).MemSac_DDI;
        
        % -------- Local memsac indicator --------
        % 1. Left and Right (Pref/null of HD task)
        MemSac_PREFmNULL_PSTH(i,:) = (group_result(i).MemSac_PSTH(:,1) - group_result(i).MemSac_PSTH(:,5)) *  sign((group_result(i).PREF_PSTH == 2)-0.5);  % (Right - Left)* Right is pref
        MemSac_PSTH_AngDiff(i,:) = mod(group_result(i).MemSac_vectSum - ((group_result(i).PREF_PSTH == 2)*0 + (group_result(i).PREF_PSTH == 1)*180),360);
        
        % 2. Pref_null DI. @HH20150524
        MemSac_PREF_Null_DI(i,:) = group_result(i).mat_raw_MemSac.PREF_NULL_DI;
        
        % 3. Actual DI of PREf/NULL of HD task.  @HH20150524
        
        % Interpolate original Memsac traces into higher spatial resolution
        
        % Moved to batch file
        %         MemSac_interp_locations = 0:5:360; % Resoluation = 1 degree
        %         MemSac_interp_PSTH = nan(size(MemSac_original_PSTH,1),length(MemSac_interp_locations));
        %
        %         for tttt = 1:size(MemSac_original_PSTH,1)
        %             MemSac_interp_PSTH(tttt,:) = interp1(MemSac_original_locations, MemSac_original_PSTH(tttt,:), MemSac_interp_locations);
        %         end
        
        % Save data for individual cell plotting
        group_result(i).MemSac_interp_locations = group_result(i).mat_raw_MemSac.MemSac_interp_locations;
        group_result(i).MemSac_interp_PSTH = group_result(i).mat_raw_MemSac.MemSac_interp_PSTH{3};
        
        % Find temporal periods
        MemSac_actual_DI_period = [3 4];  % I use memory and pre period to calculate the actual DI
        
        MemSac_temporal_Slice = group_result(i).mat_raw_MemSac.temporal_Slice;
        MemSac_align_offsets = group_result(i).mat_raw_MemSac.align_offsets;
        MemSac_align_markers = group_result(i).mat_raw_MemSac.align_markers;  
        MemSac_ts = group_result(i).MemSac_ts;
        MemSac_actual_DI_time_ind = false(size(MemSac_ts));
        
        for p_ind = 1:length(MemSac_actual_DI_period)
            
            ppp = MemSac_actual_DI_period(p_ind);
            
            if MemSac_temporal_Slice{ppp,3} == 7  % Precise windows (because MemSac_PSTH is aligned to saccade onset)
                add_ind = MemSac_temporal_Slice{ppp,1} <= MemSac_ts & MemSac_ts <= MemSac_temporal_Slice{ppp,2};
            else  % Windows that is not so precise
                % Mean shift
                meanShift = mean(MemSac_align_offsets(:,MemSac_align_markers == MemSac_temporal_Slice{ppp,3})-MemSac_align_offsets(:,MemSac_align_markers==7),1);
                add_ind = meanShift + MemSac_temporal_Slice{ppp,1} <= MemSac_ts & MemSac_ts <= meanShift + MemSac_temporal_Slice{ppp,2};
            end
            
            MemSac_actual_DI_time_ind = MemSac_actual_DI_time_ind | add_ind ;
        end
        
        % Merge the gap
        MemSac_actual_DI_time_ind(find(MemSac_actual_DI_time_ind,1):find(MemSac_actual_DI_time_ind,1,'last')) = true;
        
        %         figure(); plot(group_result(i).MemSac_ts,MemSac_interp_PSTH','-')
        %         hold on; plot(MemSac_ts(MemSac_actual_DI_time_ind),50 * ones(size(MemSac_ts(MemSac_actual_DI_time_ind))),'k.','linew',2);
        
        % Find actual DI for this cell
        PREF_target_location = group_result(i).mat_raw_PSTH.PREF_target_location; % Actual location of PREF target
        
        [~,pref] = min(abs(PREF_target_location - group_result(i).MemSac_interp_locations));
        pref = mod(pref-1,length(group_result(i).MemSac_interp_locations)-1)+1;
        null = mod(pref + (length(group_result(i).MemSac_interp_locations)-1)/2 -1, length(group_result(i).MemSac_interp_locations)-1)+1; % The opposite position
        
        MemSac_PREF_mean = mean(group_result(i).MemSac_interp_PSTH(MemSac_actual_DI_time_ind,pref),1);
        MemSac_NULL_mean = mean(group_result(i).MemSac_interp_PSTH(MemSac_actual_DI_time_ind,null),1);
        
        MemSac_actual_DI(i) = (MemSac_PREF_mean - MemSac_NULL_mean) / (MemSac_PREF_mean + MemSac_NULL_mean);
        
    end
    
    MemSac_PSTH_AngDiff(MemSac_PSTH_AngDiff > 180) = 360 - MemSac_PSTH_AngDiff(MemSac_PSTH_AngDiff > 180); % Ang diffs are less than 180
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- Choose Memsac indicator (affect all Memsac measures) -------------
% 1 = Background, 2 = LS, 3 = Mem, 4 = Pre, 5 = Co, 6 = Post

MemSac_indicator = mean(MemSac_DDI(:,[3 4]),2); % Global DDI
MemSac_indicator_txt = 'MemSac\_DDI ([3 4])';

% MemSac_indicator = MemSac_actual_DI; % Actual target locations DI
% MemSac_indicator_txt = 'MemSac\_actual\_DI ([3 4])';

% MemSac_indicator = mean(MemSac_PREF_Null_DI(:,[3 4]),2); % PREF and NULL DI (this is BAD correlated with choice weights)
% MemSac_indicator_txt = 'MemSac\_PREF\_Null\_DI';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Psychometric parameters
Psy_all = [group_result.Psy_para];
Psy_thres = reshape(Psy_all(1,:),3,[])';

bayes_pred = @(x,y) sqrt(1./(1./x.^2 + 1./y.^2));
Psy_thres_pred = bayes_pred(Psy_thres(:,1),Psy_thres(:,2));
Psy_pred_ratio = Psy_thres(:,3)./Psy_thres_pred;
Psy_pred_ratio_vestibular_visual = Psy_thres(:,1:2)./repmat(Psy_thres_pred,1,2);


select_all = []; select_sus = []; select_bottom_line = []; select_tcells = []; select_no_tcells = [];
select_psy_good = []; select_psy_bad = []; find_bottom_line = [];

Choice_pref_all = []; Choice_pref_p_value_all = []; Modality_pref_all = []; Modality_pref_p_value_all = [];
select_cpref_mpref = [];

cell_selection;

%% Cell Selection and Cell Counter
    function cell_selection()  % Cell Selection and Cell Counter

        % --------  @ HH20150413 -------- 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Limitations on Repetition number & Target first
        select_all_all_monkey = ([group_result.repN]' >= 8) & ([group_result.length_unique_stim_type]' >=3 ) & (xls_num{1}(:,header.HD_TargFirst)~=0);
        
        % + SUs
        select_sus_all_monkey = select_all_all_monkey & (xls_num{1}(:,header.Units_RealSU) > 0) & (xls_num{1}(:,header.Chan1) < 20);
        
        % Bottom line for most figures
        select_bottom_line_all_monkey = select_sus_all_monkey;
        
        % + T(ypical) Cells
%         select_tcells_all_monkey = select_sus_all_monkey & (xls_num{1}(:,header.HD_MemSac) >= 0.8);    % Manually assigned
        select_tcells_all_monkey = select_sus_all_monkey & (xls_num{1}(:,header.HD_comb_ChoicePref_p) <0.05);    % HH20160426
        % select_tcells = select_sus & mean(MemSac_DDI(:,[3 4]),2) >= 0.55;
        
        % Non-T(ypical) Cells
%         select_no_tcells_all_monkey = select_sus_all_monkey & (xls_num{1}(:,header.HD_MemSac) <0.8);   % Manually assigned
        select_no_tcells_all_monkey = select_sus_all_monkey & (xls_num{1}(:,header.HD_comb_ChoicePref_p) >=0.05);   % HH20160426
        % select_no_tcells = select_sus & mean(MemSac_DDI(:,[3 4]),2) < 0.55;
        
        % Psychometric good (prediction ratio <= 1.3)
        select_psy_good = Psy_pred_ratio <= 1.3;
        select_psy_bad = Psy_pred_ratio > 1.3;
        
        % select_memsac_good = xls_num{1}(:,header.HD_MemSac) == 1; % MemSac_DDI(:,4) >= 0.55;
        % select_memsac_bad = xls_num{1}(:,header.HD_MemSac) < 1; % MemSac_DDI(:,4) < 0.55;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -------- Count cell numbers for each monkey. HH20150723 -------- 
        n_monkey = length(monkey_included_for_loading);
        cell_nums = zeros(n_monkey + 1,3); % All units / SUs / T Cells
        for mm = 1:n_monkey
            select_monkey{mm} = xls_num{1}(:,header.Monkey) == monkey_included_for_loading(mm);
            cell_nums(mm,1) = sum(select_all_all_monkey & select_monkey{mm});
            cell_nums(mm,2) = sum(select_sus_all_monkey & select_monkey{mm});
            cell_nums(mm,3) = sum(select_tcells_all_monkey & select_monkey{mm});
        end

        % -------- Update actual dataset for analysis. HH20150723 -------- 
        monkey_included_for_analysis = monkey_included_for_loading(logical([get(findall(gcbf,'tag','Polo_data'),'value') get(findall(gcbf,'tag','Messi_data'),'value')]));
        monkey_mask_for_analysis = false(length(group_result),1);
        for mm = 1:length(monkey_included_for_analysis)
            monkey_mask_for_analysis = monkey_mask_for_analysis | (xls_num{1}(:,header.Monkey) == monkey_included_for_analysis(mm));
        end
        
        % -------- Affect all analysis below --------
        select_all = select_all_all_monkey & monkey_mask_for_analysis;
        select_sus = select_sus_all_monkey & monkey_mask_for_analysis;
        select_bottom_line = select_bottom_line_all_monkey & monkey_mask_for_analysis;   find_bottom_line = find(select_bottom_line);
        select_tcells = select_tcells_all_monkey & monkey_mask_for_analysis;
        select_no_tcells = select_no_tcells_all_monkey & monkey_mask_for_analysis;
        
        cell_nums(end,:) = [sum(select_all) sum(select_sus) sum(select_tcells)];
        
        % -------- Update cell counter ---------
        h_all = findall(gcbf,'tag','num_all_units');
        set(h_all,'string',sprintf('%5d%5d%5d\n',cell_nums'),'fontsize',20);
        %         h_su = findall(gcbf,'tag','num_su');
        %         set(h_su,'string',num2str(sum(select_sus)));
        %         h_tcell = findall(gcbf,'tag','num_t_cells');
        %         set(h_tcell,'string',num2str(sum(select_tcells)));
        
        % -------- Update/refresh some related datasets that are influenced by cell_selection ---------
        % For Choice and modality preference
        select_cpref_mpref = select_bottom_line;
        Choice_pref_all = reshape([group_result(select_cpref_mpref).ChoicePreference]',stim_type_num,[],2);  % Stim, Cell No, Pre/Post
        Choice_pref_p_value_all = reshape([group_result(select_cpref_mpref).ChoicePreference_pvalue]',stim_type_num,[],2);
        Modality_pref_all = reshape([group_result(select_cpref_mpref).ModalityPreference]',3,[],2);
        Modality_pref_p_value_all = reshape([group_result(select_cpref_mpref).ModalityPreference_pvalue]',3,[],2);
        
        select_for_SVM = select_bottom_line;
        select_for_PCA_B = select_bottom_line;
        
        PCA_A = [];  % Reset PCA_A
        PCA_B_projPC = []; % Reset PCA_B
        thres_choice = []; % Reset SVM training
        weights_PCA_B_PC = []; weights_svm_choice_mean = []; weights_TDR_PCA_SVM_mean = []; % Reset TDR
    end

%% Final Preparation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===========  Common parameters   ============
% Overall time marker (Stim on, Stim off, Sac on)
tt = reshape(cell2mat([group_result.time_marker]),6,[])';

for j = 1:2
    time_markers{j} = [mean(tt(:,(j-1)*3+1:(j-1)*3+3)); std(tt(:,(j-1)*3+1:(j-1)*3+3))];
end

p_critical = 0.01;

% =========== Data for common use ============
function_handles = [];

% For PCA_A and hotgram
j_PCA_A = 1;
PCA_A = []; PCA_A_PC = []; sort_time_interval1 = []; sort_time_interval2 = []; sort_time_interval3 = [];
A_memSac = []; A_choicediv = []; A_moddiv = []; A_CP = [];
enlarge_factor = 30; % Enlarge memsac DDIs

% For PCA_B
j_PCA_B = 1;
select_for_PCA_B = select_bottom_line;
PCA_B_time_range = min(rate_ts{j_PCA_B})+100 <= rate_ts{j_PCA_B} & rate_ts{j_PCA_B} <= time_markers{j_PCA_B}(1,3);  % Before saccade
PCA_B_times = rate_ts{j_PCA_B}(PCA_B_time_range);
denoised_dim = 9;
PCA_B = []; weights_PCA_B_PC = []; PCA_B_projPC = []; PCA_B_explained = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ================ Miscellaneous ===================

colors = [0 0 1; 1 0 0; 0 0.8 0.4; 1 0.8 0 ; 0 1 1];
set(0,'defaultAxesColorOrder',colors);
modality_diff_colors = colors;    modality_diff_colors(3,:) = [0 0 0];

transparent = 0;
figN = 2499;

marker_for_time_markers{1} = {'-','-','--'};
marker_for_time_markers{2} = {'--','--','-'};


%% ====================================== Function Handles =============================================%

% Today I reconstruct the group analysis codes in a more logical and effcient way:
% 1. I divided each figure into small nested functions in Group_HD
% 2. Group_HD returns a cell called "function_handles" back to GROUP_GUI (with all the data available for each nested function)
% 3. Now I can choose which figure(s) to plot or which function to debug in GROUP_GUI
% 4. I feel so happy.
% @HH20150425


function_handles = {
    'Rate Metrics', {
    'Correct only, all choices (conventional)', @f1p1;
    '   Different weighting methods',@f1p1p5;
    '   Single cell analysis', @f1p1p6;
    'Different headings' , @f1p2;
    'Correct / Wrong Trials', @f1p3;
    };
    
    'ROC Metrics',{
    'CP, CDiv, and MDiv',@f2p1;
    'Multisensory Enhancement of CDiv',@f2p2;
    'Easy and Difficult',@f2p3;
    };
    
    'Correlations', {
    'Mem-sac vs. CD/CP', @f3p1;
    'Choice Preference vs. Modality Preference', @f3p2;
    'Psychophysics vs. CD', @f3p3;
    };
    
    'PCA_A analysis (Eigen-feature)',{
    'Hot-gram', @f4p1;
    'Cluster and Trajectory', @f4p2;
    };
    
    'PCA_B analysis (Eigen-neuron)',{
    'Weights and correlations', @f5p1;
    '1-D Trajectory',@f5p2;
    '2-D Trajectory',@f5p3;
    };
    
    'Linear SVM decoder',{
    'Training SVM', @f6p0;
    'Weights', @f6p1;
    'Performance (overall)', @f6p2;
    'SVM weighted sum',@f6p3;
    }
    
    'Linear regression',{
    'Comb = w1 Vest + w2 Vis (Fig.5 in Gu 2008)',@f7p1;
    };
    
    'Targeted Dimensionality Reduction',{
    'PCA + SVM: weights', @f8p1;
    'PCA + SVM: all correct + all angles',@f8p2;
    }
    
    'Others',{
    'Cell Counter',@f9p1;
    'Target first vs. Target last',@f9p2;
    'Test',@f9p9;
    };
    
    'NoShow',{@cell_selection};
    
    };

%% ====================================== Function Definitions =============================================%

    function f1p1(debug, methods_of_select)      % Rate 1. Correct only, all choices
        if debug  ; dbstack;   keyboard;      end
        
        % Reuse elsewhere than directly from Group_GUI
        if nargin < 2
            methods_of_select = {
                select_bottom_line, 'All cells';
                select_tcells, 'Typical cells';
                select_no_tcells, 'Non-typical cells'};
        end
        
        %% ------- Averaged norm PSTH --------
        
        set(figure(999),'name','Average norm PSTH (Correct only, all choices)','pos',[27 63 1449 892]); clf
        h_subplot = tight_subplot(2,size(methods_of_select,1),[0.11 0.05],[0.05 0.1],0.07);
        
        %         for j = 1:2
        j = 1;
        
        for ms = 1:size(methods_of_select,1)
            
            % --- Pref and Null ---
            SeriesComparison({PSTH_all_Norm{1}(methods_of_select{ms,1},:,:), PSTH_all_Norm{2}(methods_of_select{ms,1},:,:)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',mat2cell(reshape(repmat(colors,1,2)',3,[])',ones(stim_type_num*2,1)),'LineStyles',{'-','--','-','--','-','--'},...
                'ErrorBar',4,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot((ms-1)*2+1),...
                'CompareIndex',[1:2:2*stim_type_num;2:2:2*stim_type_num],...
                'CompareColor',mat2cell(colors,ones(stim_type_num,1)),'Border',[1800 -350]);
            
%             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end

            xlabel('Time (ms)');
            legend off;
            title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}>0))]);
            axis tight;
            
            axis tight;
            xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 0.7]);
            
            % Gaussian vel
%             plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/6,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
            
            % --- Difference (Pref - Null) ---
            PSTH_all_Norm_PrefminusNull{1} = PSTH_all_Norm{1}(methods_of_select{ms,1},:,1:2:end)...
                                           - PSTH_all_Norm{1}(methods_of_select{ms,1},:,2:2:end);
            PSTH_all_Norm_PrefminusNull{2} = PSTH_all_Norm{2}(methods_of_select{ms,1},:,1:2:end)...
                                           - PSTH_all_Norm{2}(methods_of_select{ms,1},:,2:2:end);
            
            SeriesComparison({PSTH_all_Norm_PrefminusNull{1}, PSTH_all_Norm_PrefminusNull{2}},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',mat2cell(colors,ones(stim_type_num,1)),'LineStyles',{'-'},...
                'ErrorBar',4,'Xlabel',[],'Ylabel','\Delta Norm firing rate (pref-null)','axes',h_subplot((ms-1)*2+2),...
                'CompareIndex',[1:stim_type_num,3,4,5;1:stim_type_num,4,5,3],...
                'CompareColor',[mat2cell(colors,ones(stim_type_num,1));'b';'r';'g'],'Border',[1800 -350]);
            
            %             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            xlabel('Time (ms)');
            legend off;
            
            axis tight;
            xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 0.7]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/6,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            plot(xlim,[0 0],'k--');
            
            legend off;

        end
        %         end
        SetFigure(15);
        
        
         %% ------- Averaged raw PSTH --------
        set(figure(998),'name','Average raw PSTH (Correct only, all choices)','pos',[27 63 1449 892]); clf
        h_subplot = tight_subplot(2,size(methods_of_select,1),[0.11 0.05],[0.05 0.1],0.07);
                
        %         for j = 1:2
        j = 1;
        for ms = 1:size(methods_of_select,1)
            
            % --- Pref and Null ---
            
            SeriesComparison({PSTH_all_raw{1}(methods_of_select{ms,1},:,:), PSTH_all_raw{2}(methods_of_select{ms,1},:,:)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',mat2cell(reshape(repmat(colors,1,2)',3,[])',ones(stim_type_num*2,1)),'LineStyles',{'-','--','-','--','-','--'},...
                'ErrorBar',4,'Xlabel',[],'Ylabel','Raw firing rate','axes',h_subplot((ms-1)*2+1),...
                'CompareIndex',[1:2:2*stim_type_num;2:2:2*stim_type_num],...
                'CompareColor',mat2cell(colors,ones(stim_type_num,1)),'Border',[1800 -350]);
            
%             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end

            xlabel('Time (ms)');
            legend off;
            title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}>0))]);
            axis tight;
            
            axis tight;
            xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 0.7]);
            
            % Gaussian vel
%             plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/6,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
            
            % --- Difference (Pref - Null) ---
            PSTH_all_raw_PrefminusNull{1} = PSTH_all_raw{1}(methods_of_select{ms,1},:,1:2:end)...
                                           - PSTH_all_raw{1}(methods_of_select{ms,1},:,2:2:end);
            PSTH_all_raw_PrefminusNull{2} = PSTH_all_raw{2}(methods_of_select{ms,1},:,1:2:end)...
                                           - PSTH_all_raw{2}(methods_of_select{ms,1},:,2:2:end);
            
            SeriesComparison({PSTH_all_raw_PrefminusNull{1}, PSTH_all_raw_PrefminusNull{2}},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',mat2cell(colors,ones(stim_type_num,1)),'LineStyles',{'-'},...
                'ErrorBar',4,'Xlabel',[],'Ylabel','\Delta Raw firing rate (pref-null)','axes',h_subplot((ms-1)*2+2),...
                'CompareIndex',[1:stim_type_num,4;1:stim_type_num,5],...
                'CompareColor',[mat2cell(colors,ones(stim_type_num,1));'k'],'Border',[1800 -350]);
            
            %             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            xlabel('Time (ms)');
            legend off;
            
            axis tight;
            xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 0.7]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/6,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            plot(xlim,[0 0],'k--');
            
            legend off;
        end
        %         end
        SetFigure(15);
    end

    function f1p1p5(debug)
        if debug  ; dbstack;   keyboard;      end
        
        % ------- Weighted by different methods ------
        % Straight mean
        weights_straight_mean = ones(sum(select_bottom_line),1);
        Weighted_sum_PSTH( weights_straight_mean,{'Weighted by straight mean'},select_bottom_line);
        title('Weighted by same weights (straight mean)');
%         Weights_Property_Correlation( weights_straight_mean /sum(weights_straight_mean) ,...
%             {'Weighted by straight mean','Weighted by straight mean'},select_bottom_line);
        
        % Norm by dynamic range
        weights_norm = weights_normalized_PSTH(select_bottom_line)';
        Weighted_sum_PSTH( weights_norm,{'Weighted by 1/dynamic range'},select_bottom_line);
        title('Weighted by 1/dynamic range');
%         Weights_Property_Correlation(weights_norm/sum(weights_norm),...
%             {'Weighted by 1/dynamic range','Weighted by 1/dynamic range'},select_bottom_line);
        
        % Norm by dynamic range with cell selection
        weights_norm_tcells = weights_norm;
        weights_norm_tcells(~select_tcells(select_bottom_line)) = 0;
        Weighted_sum_PSTH( weights_norm_tcells,{'Weighted by 1/dynamic range with cell selection'},select_bottom_line);
        title('Weighted by 1/dynamic range with cell selection');
%         Weights_Property_Correlation(weights_norm_tcells/sum(weights_norm_tcells),...
%             {'Weighted by 1/dynamic range with cell selection','Weighted by 1/dynamic range with cell selection'},select_bottom_line);
        
    end
    function f1p1p6(debug)
        if debug  ; dbstack;   keyboard;      end
        %%
        % --- Difference (Pref - Null) ---
        % All angles
        PSTH_all_raw_PrefminusNull{1} = PSTH_all_raw{1}(select_bottom_line,:,1:2:end)...
            - PSTH_all_raw{1}(select_bottom_line,:,2:2:end);
        PSTH_all_raw_PrefminusNull{2} = PSTH_all_raw{2}(select_bottom_line,:,1:2:end)...
            - PSTH_all_raw{2}(select_bottom_line,:,2:2:end);

        % 8-degree      
%         PSTH_all_raw_PrefminusNull{1} = squeeze(PSTH_correct_angles_raw{1}(select_bottom_line,:,7,:)...
%             - PSTH_correct_angles_raw{1}(select_bottom_line,:,8,:));
%         PSTH_all_raw_PrefminusNull{2} = squeeze(PSTH_correct_angles_raw{2}(select_bottom_line,:,7,:)...
%             - PSTH_correct_angles_raw{2}(select_bottom_line,:,8,:));
        
        [~,plot_order] = sort(mean(abs(Choice_pref_all(:,:,1))),'descend');
        
        %%{
        % === All traces ===
        
        set(figure(799),'name','Average norm PSTH (Correct only, all choices)','pos',[27 63 1449 892]); clf
        [~,h_subplot] = tight_subplot(4,5,[0.03 0.02]);
        
        for nn = 1:length(plot_order)
            nnn = plot_order(nn);
            
            axes(h_subplot(nn));
            SeriesComparison({PSTH_all_raw_PrefminusNull{1}(nnn,:,:), PSTH_all_raw_PrefminusNull{2}(nnn,:,:)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',mat2cell(colors,ones(stim_type_num,1)),'LineStyles',{'-'},...
                'ErrorBar',0,'Xlabel',[],'Ylabel','Raw firing rate','axes',h_subplot(nn),...
                'Border',[1800 -350]);
            
            title(sprintf('Cell %g (p-values: %g, %g, %g, %g, %g)',nn,Choice_pref_p_value_all(:,nnn,1)'));
            legend off;
            xlim([0 1750]); axis tight;
            ylabel('');
            if nn ~= 15
                set(gca,'xticklabel','');
            end
            plot(xlim,[0 0],'k--');
            
        end
        delete(h_subplot(end-1:end));
        %}
        
        % === Curve fitting ===
        set(figure(798),'name','Average norm PSTH (Correct only, all choices)','pos',[27 63 1449 892]); clf
        [~,h_subplot] = tight_subplot(4,5,[0.03 0.02]);
        
        % Define sigmoid function for fitting
        func_sigmoid = fittype('aa*1/(1+exp(-(x-tt)/bb))','coeff',{'aa','bb','tt'});
        options = fitoptions(func_sigmoid);
        options.StartPoint = [10 100 600];
        
        aas = nan(length(plot_order),stim_type_num);
        bbs = aas; tts = aas;
        
        for nn = 1:length(plot_order)
            nnn = plot_order(nn);
            
            axes(h_subplot(nn));

            this_ys = squeeze(PSTH_all_raw_PrefminusNull{1}(nnn,rate_ts{1}<=1750,:));
            plot(rate_ts{1}(rate_ts{1}<=1750),this_ys,'linew',2); hold on;
            
            for kk = 1:stim_type_num
                % Fitting range
                this_ys = PSTH_all_raw_PrefminusNull{1}(nnn,:,kk);
                this_ys_smooth = smooth(this_ys,5);
                
                if kk == 2
%                     fit_t_range = 0 <= rate_ts{1} & rate_ts{1} <= 1500 ;
                    t_end = rate_ts{1}(this_ys_smooth == max(this_ys_smooth(0 <= rate_ts{1} & rate_ts{1}<=1500)));
                else
%                     fit_t_range = 0 <= rate_ts{1} & rate_ts{1} <= 1300 ;
                    t_end = rate_ts{1}(this_ys_smooth == max(this_ys_smooth(0 <= rate_ts{1} & rate_ts{1}<=1200)));
                end
                
                fit_t_range = 0 <= rate_ts{1} & rate_ts{1} <= t_end + 200 ;
 
                this_t = rate_ts{1}(fit_t_range);
                this_ys = this_ys(fit_t_range);  

                % Fit and plot
                try
                    [fitresult{nn,kk}, gof{nn,kk}] = fit( this_t(:) , this_ys(:), func_sigmoid, options );
                    plot(this_t, fitresult{nn,kk}(this_t),'color',colors(kk,:),'linewidth',3);
                    
                    aas(nn,kk) = fitresult{nn,kk}.aa;
                    bbs(nn,kk) = fitresult{nn,kk}.bb;
                    tts(nn,kk) = fitresult{nn,kk}.tt;
                    rrs(nn,kk) = gof{nn,kk}.adjrsquare;
                catch
                end
            end
            
            title(sprintf('Cell %g (p-values: %g, %g, %g, %g, %g)',nn,Choice_pref_p_value_all(:,nnn,1)'));
            legend off;
            xlim([0 1750]); axis tight;
            ylabel('');
            if nn ~= 15
                set(gca,'xticklabel','');
            end
            plot(xlim,[0 0],'k--');
            
        end
        delete(h_subplot(end-1:end));
        
        %% === Statistics
        select_cells = [1:6 8:11 14];  exclude_cells = setdiff(1:18,select_cells);
        result_bb = BarComparison(bbs(select_cells,:),'figN',556,'Colors',mat2cell(colors,ones(stim_type_num,1)));
        title('bbs'); ylim([1 200]); plot(1:5,bbs(exclude_cells,:),'o:','color',[0.7 0.7 0.7]);
        
        result_tt = BarComparison(tts(select_cells,:),'figN',557,'Colors',mat2cell(colors,ones(stim_type_num,1)));
        title('tts'); ylim([1 1400]); plot(1:5,tts(exclude_cells,:),'o:','color',[0.7 0.7 0.7]);
        
        result_aa = BarComparison(aas(select_cells,:),'figN',558,'Colors',mat2cell(colors,ones(stim_type_num,1)));
        title('aas'); ylim([1 40]); plot(1:5,aas(exclude_cells,:),'o:','color',[0.7 0.7 0.7]);
        
        %% === Redraw traces
%         t = 1:1500; figure(559); clf;
%         for nn = 1:18
%             for kk = 1:5
%                 subplot(1,5,kk); hold on;
%                 plot(t,fitresult{nn,kk}(t));
%             end
%         end

        %% === Replot average
        %% ------- Averaged raw PSTH --------
        
        methods_of_select = {
            select_cells, 'non "X" cells'
            exclude_cells, '"X" cells'
            [select_cells exclude_cells],'All cells'
            };
        
        f1p1(0, methods_of_select);
        
    end

    function f1p2(debug)      % Rate 2. Different headings
        %% Different headings
        
        if debug
            dbstack;
            keyboard;
        end
        
        methods_of_select = {
            select_bottom_line, 'All cells';
            select_tcells, 'Typical cells';
            select_no_tcells, 'Non-typical cells'};
        
        colors_angles = colormap(gray);
        colors_angles = colors_angles(round(linspace(length(colors_angles)-15,1,5)),:);
        colors_angles = reshape(repmat(colors_angles,1,2)',3,[])';
        colors_angles = mat2cell(colors_angles,ones(10,1));
        
        for ms = 1:3
            set(figure(999-ms),'name',['Average PSTH (correct only, all headings), ' methods_of_select{ms,2}],'pos',[27 57 919 898]); clf
            h_subplot = tight_subplot(stim_type_num,2,[0.05 0.1],[0.05 0.05],[0.05 0.05]);
            
            for k = 1:stim_type_num
                
                colors_angles = colormap(gray); 
                colors_angles = ones(size(colors_angles,1),3) - colors_angles .* repmat([1 1 1]-colors(k,:),size(colors_angles,1),1);
                colors_angles = colors_angles(round(linspace(20,length(colors_angles),5)),:);
                colors_angles = reshape(repmat(colors_angles,1,2)',3,[])';
                colors_angles = mat2cell(colors_angles,ones(10,1));
                
                % --- Ramping with different angles ---
                for j = 1:2
                    yyy{j} = PSTH_correct_angles_Norm{j}(methods_of_select{ms,1},:,:,k);
%                     yyy{j} = PSTH_correct_angles_raw{j}(methods_of_select{ms,1},:,:,k);
                    ttt{j} = rate_ts{j};
                    
                    yyy_diff{k}{j} =  yyy{j}(:,:,1:2:end) - yyy{j}(:,:,2:2:end);
                end
                
                
                h = SeriesComparison(yyy,{ttt{1} ttt{2} time_markers},...
                    'Colors',colors_angles,'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k));
                
                if k < stim_type_num ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                if k == 1
                    title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  ylim([0.1 max(ylim)]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
                
                % ----  Difference ---
                
                h = SeriesComparison(yyy_diff{k},{ttt{1} ttt{2} time_markers},...
                    'Colors',colors_angles(1:2:end),'LineStyles',{'-'},...
                    'ErrorBar',2,'Xlabel',[],'Ylabel','Diff','axes',h_subplot(k+(2-1)*stim_type_num));
                
                if k < stim_type_num ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]); ylim([-.1 max(ylim)]); plot(xlim,[0 0],'k--');
                
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
            end
            SetFigure(15); drawnow;
            
            % --- Multisensory enhancement of different angles --- HH20160126
            set(figure(995-ms),'name',['Enhancement for different angles (correct only, all headings), ' methods_of_select{ms,2}],'pos',[856 273 818 679]); clf
            h_subplot = tight_subplot(2,2,[0.05 0.1],0.1,[0.12 0.03]); 
            unique_abs_heading = unique(abs(group_result(end).mat_raw_PSTH.heading_per_trial));
            
            for aa = 1:size(yyy_diff{1}{1},3) % Different angles
                yyy_diff_this_angle_all_stim_type = [];
                for j = 1:2
                    for k = 1:stim_type_num
                        yyy_diff_this_angle_all_stim_type{j}(:,:,k) = yyy_diff{k}{j}(:,:,aa); % Reorganize data
                    end
                end
                
                h = SeriesComparison(yyy_diff_this_angle_all_stim_type, {ttt{1} ttt{2} time_markers},...
                    'Colors',mat2cell(colors,ones(stim_type_num,1)),'LineStyles',{'-'},...
                    'ErrorBar',6,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(aa),'Transparent',0.5,...
                    'CompareIndex',[1:stim_type_num,3,4,5;1:stim_type_num,4,5,3],...
                    'CompareColor',[mat2cell(colors,ones(stim_type_num,1));'b';'r';'g'],...
                    'Border',[1800 -350]);
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]); ylim([-.1 max(ylim)]); plot(xlim,[0 0],'k--');
                
               
                title(['|heading| = ',num2str(unique_abs_heading(aa))]);
                if k == 1
                    title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                if aa ~= 4 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end

                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                legend off;    SetFigure(15); drawnow
            end
            
        end
        
    end
    function f1p3(debug)      % Rate 3. Correct / Wrong Trials
        if debug
            dbstack;
            keyboard;
        end
        %%
        methods_of_select = {
            select_bottom_line, 'All cells';
            select_tcells, 'Typical cells';
            select_no_tcells, 'Non-typical cells'};
        
        unique_heading = group_result(end).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
        unique_heading_for_correct_wrong = unique_heading(unique_heading > 0);
        
        for ms = 1:3
            set(figure(999-ms),'name',['Average PSTH (Correct vs Wrong), ' methods_of_select{ms,2}],'pos',[18 67 1645 898]); clf
            h_subplot = tight_subplot(stim_type_num,1 + length(unique_heading_for_correct_wrong),[0.05 0.02],[0.05 0.05]);
            
            j = 1;
            for k = 1:stim_type_num
                
                % ------ All correct and wrong ------ %
                h = SeriesComparison({PSTH_outcomes_Norm{1}(methods_of_select{ms,1},:,:,k) PSTH_outcomes_Norm{2}(methods_of_select{ms,1},:,:,k)},...
                    {rate_ts{1} rate_ts{2} time_markers},...
                    'Colors',{'k','k','m','m'},'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                % legend off;
                
                if k == 1
                    title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                for tt = 1:3
                    plot([1 1] * time_markers{j}(1,tt),[0 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
                end
                
                xlim([rate_ts{j}(10) rate_ts{j}(end-10)]);  ylim([.0 .7]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                % legend off;
                
                % ------ Correct and wrong for each heading ------ %    @HH20150523
                
                for hh = 1:length(unique_heading_for_correct_wrong)
                    
                    % Reorganize to fit "correct +, correct -, wrong +, wrong -" for each heading
                    correct_pref_ind = 2*sum(unique_heading==0) + hh*2-1 ; % Index in PSTH_correct_angles_norm (exclude first two 0 heading)
                    correct_null_ind = 2*sum(unique_heading==0) + hh*2 ;  % Index in PSTH_correct_angles_norm
                    wrong_pref_ind =  hh*2-1 ; % Index in PSTH_wrong_angles_norm
                    wrong_null_ind =  hh*2 ; % Index in PSTH_wrong_angles_norm
                    
                    for jjj = 1:2
                        PSTH_correct_wrong_this{jjj} = ...
                            cat(3, PSTH_correct_angles_Norm{jjj}(methods_of_select{ms,1},:,[correct_pref_ind correct_null_ind],k),...
                            PSTH_wrong_angles_Norm{jjj}(methods_of_select{ms,1},:,[wrong_pref_ind wrong_null_ind],k));
                    end
                    
                    SeriesComparison({PSTH_correct_wrong_this{1} PSTH_correct_wrong_this{2}},...
                        {rate_ts{1} rate_ts{2} time_markers},...
                        'Colors',{'k','k','m','m'},'LineStyles',{'-','--'},...
                        'ErrorBar',0,'Xlabel',[],'Ylabel',[],'axes',h_subplot(hh * stim_type_num + k));
                    
                    if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                    set(gca,'ytick',[]);
                    
                    % legend off;
                    
                    if k == 1
                        title([methods_of_select{ms,2} ', ' num2str(unique_heading_for_correct_wrong(hh))  ' n = ' num2str(sum(methods_of_select{ms,1}))]);
                    end
                    
                    for tt = 1:3
                        plot([1 1] * time_markers{j}(1,tt),[0 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
                    end
                    
                    xlim([rate_ts{j}(10) rate_ts{j}(end-10)]);  ylim([.0 .7]);
                    
                    % Gaussian vel
                    plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                end
                
            end
            
            SetFigure(10);
        end
    end

    function f2p1(debug)      % ROC 1. CP, CDiv, MDiv
        if debug
            dbstack;
            keyboard;
        end
        
        %%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         select_for_CP = select_bottom_line;
%         select_for_div = select_bottom_line;
        
        select_for_CP = select_tcells;
        select_for_div = select_tcells;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% CP: Time-course

        
        set(figure(2099+figN),'name','CP','pos',[287 483 1161 475]); figN = figN+1; clf;
        
        for j = 1:2
            subplot(1,2,j);
            
            for k = 1:stim_type_num
                
                %     selectCells_notNaN =  select_for_CP & (~isnan(CP(:,1,k)));
                ys = mean(CP{j}(select_for_CP,:,k));
                errors = std(CP{j}(select_for_CP,:,k))/sqrt(sum(select_for_CP));
                h = shadedErrorBar(CP_ts{j},ys,errors,{'Color',colors(k,:)},transparent);
                set(h.mainLine,'LineWidth',2);
                hold on;
                
            end
            
            xlabel(['Center of ' num2str(group_result(end).mat_raw_PSTH.binSize_CP) ' ms time window']);
            
            title(sprintf('Grand CP (N = %g)',sum(select_for_CP)));
            xlim([CP_ts{j}(1) CP_ts{j}(end)]);
            ylim([0.45 0.75]);
            
            for tt = 1:3
                plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
            end
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',3,'color',[0.6 0.6 0.6]);
            
            plot(xlim,[0.5 0.5],'k--');
            
        end
        
        SetFigure(17);
        
        %% CDiv: Time-course 
        
        % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) >1) & (xls_num{1}(:,header.HD_TargFirst)~=0);
        
        
        set(figure(2099+figN),'name','CDiv and MDiv','pos',[62 109 1528 584]); figN = figN+1; clf;
        
        for j = 1:2
            
            subplot(1,2,j); hold on;
            
            % Modality divergence (abs(vis-vest))
            
            % The following two lines are problematic. I should flip each ModDiv
            % accoring to its "preferred modality" and then average them directly, not
            % use abs(...)!
            
            % ys_MD = mean(abs(ModDiv_All(select_cd_all,:,1)));
            % err_MD = std(abs(ModDiv_All(select_cd_all,:,1)))/sqrt(sum(select_cd_all));
            
            % Comb v.s. Vestibular
            ModDiv_to_average = ModDiv_All{j}(select_for_div,:,2);
            ModDiv_to_average = ModDiv_to_average.* repmat(sign(sum(ModDiv_to_average,2)),1,size(ModDiv_to_average,2));
            
            ys_MD = mean(ModDiv_to_average);
            err_MD = std(ModDiv_to_average)/sqrt(size(ModDiv_to_average,1));
            
            h = shadedErrorBar(rate_ts{j},ys_MD,err_MD,{'Color',[0.5 0.5 0.5]},transparent);
            set(h.mainLine,'LineWidth',2);
            
            % Comb v.s. Visual
            ModDiv_to_average = ModDiv_All{j}(select_for_div,:,3);
            ModDiv_to_average = ModDiv_to_average.* repmat(sign(sum(ModDiv_to_average,2)),1,size(ModDiv_to_average,2));
            
            ys_MD = mean(ModDiv_to_average);
            err_MD = std(ModDiv_to_average)/sqrt(size(ModDiv_to_average,1));
            
            h = shadedErrorBar(rate_ts{j},ys_MD,err_MD,{'Color',[0.5 0.5 0.5]},transparent);
            set(h.mainLine,'LineWidth',2,'linestyle','-.'); hold on;
            
            
            % Choice divergence
            for k = 1:stim_type_num
                selectCells_notNaN =  select_for_div & (~isnan(ChoiceDiv_All{j}(:,1,k)));
                
                ys_CD{j}(k,:) = mean(ChoiceDiv_All{j}(selectCells_notNaN,:,k));
                errors = std(ChoiceDiv_All{j}(selectCells_notNaN,:,k))/sqrt(sum(selectCells_notNaN));
                h = shadedErrorBar(rate_ts{j},ys_CD{j}(k,:),errors,{'Color',colors(k,:),'linestyle','-'},transparent);
                set(h.mainLine,'LineWidth',2);
            end
            
            % Linear sum of vest and vis
            plot(rate_ts{j}, sum(ys_CD{j}(1:2,:)),'k-');
            
            axis tight;  ylim([-0.1 0.35]);
            plot(xlim,[0 0],'k--');
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',3,'color',[0.6 0.6 0.6]);
            
            title(sprintf('Choice divergence (SU, N = %g)',sum(select_for_div)));
            
            for tt = 1:3
                plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
            end
            
        end
        
        SetFigure(17);
        
        %% For clearer drawing. HH20160213
        
        SeriesComparison({ChoiceDiv_All{1}(select_for_div,:,:) ChoiceDiv_All{2}(select_for_div,:,:)},...
            {rate_ts{1} rate_ts{2} time_markers},...
            'Colors',mat2cell(colors,ones(stim_type_num,1)),'LineStyles',{'-'},...
            'ErrorBar',2,'Xlabel',[],'Ylabel',[],'transparent',0.5);
        xlim([-300 2300]); ylim([-0.1 0.4]); legend off;  plot(xlim,[0 0],'k--'); 
        title(sprintf('Choice divergence (SU, N = %g)',sum(select_for_div)));
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
        
        SetFigure(17);
        
%         %% A or V? HH20151220
%         figure(); plot(rate_ts{1},ys_CD{1}');
%         hold on;
%         plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',3,'color',[0.6 0.6 0.6]);
%         dt_Gauss = Gauss_vel(2,1)-Gauss_vel(1,1);
%         Gauss_acc = diff(Gauss_vel(:,2));
%         plot(Gauss_vel(2:end,1),Gauss_acc);
%         plot(Gauss_vel(1:end,1),[0 ;cumsum(abs(Gauss_acc))]/10,'k','linew',2);
%         plot(Gauss_vel(1:end,1),[cumsum(Gauss_vel(:,2))]/1000,'k','linew',2);
%         
%         
%         %% Linear fitting combined trace
%         ts = rate_ts{1};
%         t_select = rate_ts{1}> 0 & rate_ts{1} <=1500;
%    
%         ramping1 = ys_CD{1}(1,t_select);
%         ramping2 = ys_CD{1}(2,t_select);
%         ramping3 = ys_CD{1}(3,t_select);
%         
%         w = fminsearch(@(w) sum((w(1)*ramping1 + w(2)*ramping2 - ramping3).^2), [.5 .5])
%         
%         figure(); plot(ts(t_select),ramping1,'b',ts(t_select),ramping2,'r',ts(t_select),ramping3,'g',ts(t_select),ramping1*w(1)+ramping2*w(2),'k');
        
        
    end
    function f2p2(debug)      % ROC 2. CDiv: Multisensory Enhancement
        if debug
            dbstack;
            keyboard;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        select_for_div =  select_bottom_line;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% CDiv: Multisensory Enhancement
        % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) > 0) & (xls_num{1}(:,header.HD_TargFirst)~=0);% ...
        % & (~isnan(ChoiceDiv_All(:,1,k))) ;% & MemSac_DDI(:,4)<0.55 ;
        
        set(figure(2099+figN),'name','Multisensory Enhancement of CDiv','pos',[12 276 1272 684]); clf; figN = figN+1;
        
        for j = 1:2
            for k = 1:3
                subplot(2,3,k + (j-1)*3);
                
                select_this_k = select_for_div & (~isnan(ChoiceDiv_ModDiffer{j}(:,1,k)));
                
                ys = mean(ChoiceDiv_ModDiffer{j}(select_this_k,:,k));
                errors = std(ChoiceDiv_ModDiffer{j}(select_this_k,:,k))/sqrt(sum(select_this_k));
                [~,ps] = ttest(ChoiceDiv_ModDiffer{j}(select_this_k,:,k));  % p_value
                
                h = shadedErrorBar(rate_ts{j},ys,errors,{'Color',modality_diff_colors(k,:),'Linestyle','-'},transparent);
                set(h.mainLine,'LineWidth',2); xlim([time_markers{j}(1,1)-200 time_markers{j}(1,3)+200]); ylim([-0.05 0.1]);
                hold on;
                
                plot(rate_ts{j},ps,'k');
                
                if sum(ps < p_critical) > 0
                    plot(rate_ts{j}(ps<p_critical),max(ylim)*.9,'.','color',modality_diff_colors(k,:),'linew',5);
                end
                
                plot(xlim,[0 0],'k--');
                for tt = 1:3
                    plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
                end
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
                
            end
            
        end
        
        title(sprintf('3-1 3-2 1-2 (N = %g)',sum(select_this_k)));
        SetFigure(15);
    end
    function f2p3(debug)      % ROC 3. Easy and Difficult
        if debug
            dbstack;
            keyboard;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        select_for_div =  select_tcells; 
%         select_for_div =  select_bottom_line;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% CDiv: Easy AND difficult
        
        % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) > 0) & (xls_num{1}(:,header.HD_TargFirst)~=0);% ...
        % & (~isnan(ChoiceDiv_All(:,1,k))) ;% & MemSac_DDI(:,4)<0.55;
        
        set(figure(2099+figN),'pos',[12 276 1272 684],'name','Easy & difficult of CDiv'); clf; figN = figN+1;
        
        for j = 1:2
            for k = 1:3
                subplot(2,3,k + (j-1)*3);
                
                select_this_k = select_for_div & (~isnan(ChoiceDiv_All{j}(:,1,k))) ;
                
                ys = mean(ChoiceDiv_Easy{j}(select_this_k,:,k));
                errors = std(ChoiceDiv_Easy{j}(select_this_k,:,k))/sqrt(sum(select_this_k));
                h = shadedErrorBar(rate_ts{j},ys,errors,{'Color',colors(k,:)},transparent);
                set(h.mainLine,'LineWidth',2);  hold on;
                
                ys = mean(ChoiceDiv_Difficult{j}(select_this_k,:,k));
                errors = std(ChoiceDiv_Difficult{j}(select_this_k,:,k))/sqrt(sum(select_this_k));
                h = shadedErrorBar(rate_ts{j},ys,errors,{'Color',colors(k,:)*0.4 + [0.6 0.6 0.6]},transparent);
                set(h.mainLine,'LineWidth',2);  hold on;
                
                xlim([time_markers{j}(1,1)-200 time_markers{j}(1,3)+200]); ylim([-0.1 0.25]); ylim([-0.05 0.25]);
                
                plot(xlim,[0 0],'k--');
                for tt = 1:3
                    plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
                end
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/8 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
                
            end
        end
        
        title(sprintf('Easy - Difficult (N = %g)',sum(select_this_k)));
        SetFigure(15);
        
        %% CDiv: Easy-difficult
        
        % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) > 0) & (xls_num{1}(:,header.HD_TargFirst)~=0);% ...
        % & (~isnan(ChoiceDiv_All(:,1,k))) ;% & MemSac_DDI(:,4)<0.55;
        
        set(figure(2099+figN),'pos',[12 276 1272 684],'name','Easy-difficult of CDiv'); clf; figN = figN+1;
        
        for j = 1:2
            for k = 1:3
                subplot(2,3,k + (j-1)*3);
                
                select_this_k = select_for_div & (~isnan(ChoiceDiv_All{j}(:,1,k)));
                
                ys = mean(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k));
                errors = std(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k))/sqrt(sum(select_this_k));
                [~,ps] = ttest(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k));  % p_value
                
                h = shadedErrorBar(rate_ts{j},ys,errors,{'Color',colors(k,:)},transparent);
                set(h.mainLine,'LineWidth',2); ; hold on;
                xlim([time_markers{j}(1,1)-200 time_markers{j}(1,3)+200]); ylim([-0.1 0.25]); ylim([-0.05 0.09]);
                
                plot(rate_ts{j},ps,'k');
                
                if sum(ps < p_critical)>0
                    plot(rate_ts{j}(ps < p_critical),max(ylim)*.9,'.','color',modality_diff_colors(k,:));
                    ps(ps<p_critical)
                end
                
                plot(xlim,[0 0],'k--');
                for tt = 1:3
                    plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
                end
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
                
            end
        end
        
        title(sprintf('Easy - Difficult (N = %g)',sum(select_this_k)));
        SetFigure(15);
        
        %% Tuning: At 3 different phases. HH20150414
        % Now this part is dirty and quick.
        % Here are several things to be done:
        %  1. Align to stimulus onset AND saccade onset  (done)
        %  2. Deal with different heading angles
        %  3. Put this part elsewhere
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
        select_for_tuning = select_for_div; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        find_for_tuning = find(select_for_tuning);
        
        % Time centers (of CP time windows, width = 500 ms)
%         t_stim_center_earlier = mean(time_markers{j}(1,1:2)) -100; % Stim center + sensory delay
%         t_stim_center = mean(time_markers{j}(1,1:2)) + 150; % Stim center + sensory delay
%         t_pre_center = mean(time_markers{j}(1,3)) - group_result(end).mat_raw_PSTH.binSize_CP/2;   % Pre-sac epoch
%         t_post_center = mean(time_markers{j}(1,3)) + group_result(end).mat_raw_PSTH.binSize_CP/2;  % Post-sac epoch
        t_stim_center = 585; % Stim center + sensory delay
        t_pre_center = 900;   % Pre-sac epoch
        t_post_center = 1500;  % Post-sac epoch
        
%         [~,center_t_ind_earlier] = min(abs(CP_ts{j} - t_stim_center_earlier));
        [~,center_t_ind] = min(abs(CP_ts{j} - t_stim_center));
        [~,pre_t_ind] = min(abs(CP_ts{j} - t_pre_center));
        [~,post_t_ind] = min(abs(CP_ts{j} - t_post_center));
        
        tuning_time_phase = [center_t_ind pre_t_ind post_t_ind];    
        tuning_time_phase_title = {'Vest largest', 'Vis largest', 'Pre-sac'};
        
        % Suppose all the unique_headings are the same!!
        unique_heading = group_result(find_for_tuning(end)).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
        % HH20160214: Patch calculation for LEFT & RIGHT choices of 0 headings
        zero_index = find(unique_heading == 0);
        unique_heading_two_zeros = [unique_heading(1:zero_index-1); -eps; eps ;unique_heading(zero_index+1:end)];  
%         unique_heading_two_zeros = [unique_heading(1:zero_index-1) ;unique_heading(zero_index+1:end)];  % No zero version

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for tuning curve normalization
        
        %         time_for_dynamic_range = find(CP_ts{j} < time_markers{j}(1,2)); % all pre-sac times
        time_for_dynamic_range = [center_t_ind - 5 center_t_ind + 5]; % Around stimulus center
        modalities_share_same_dynamic_range = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Pack tuning data
        for i = 1:sum(select_for_tuning)  % For cells
            
            % Only include three modalities
            if length(group_result(find_for_tuning(i)).mat_raw_PSTH.unique_stim_type) < 3
                continue;
            end
            
            % Align preferred direcions to 'rightward'. So if PREF == 1 for a certain cell, we
            % flip all its tuning curve.
            flip_tuning = group_result(find_for_tuning(i)).PREF_PSTH == 1;
            
            % Find dynamic range
            dynamic_min_own_range = inf * [1 1 1]; % Each modality have its own dynamic range
            dynamic_max_own_range = -inf * [1 1 1];
            
            dynamic_min = inf; % Each modality have its own dynamic range
            dynamic_max = -inf;
            
            % Pack data (all tuning phases)
            
            for pp = 1:length(CP_ts{j})
                for k = 1:stim_type_num
                    this_raw = group_result(find_for_tuning(i)).mat_raw_PSTH.CP{j,k}.raw_CP_result{pp};
                    this_tuning = this_raw.Neu_tuning(:,2);
                    this_tuning_correctonly = this_raw.Neu_tuning_correctonly(:,2);
                    
                    % HH20160214: Patch calculation for LEFT & RIGHT choices of 0 headings (because I failed to do this
                    % in the original CP_HH and I feel hesitant to redo all the batch files from A to Z right now...).
                    % Note, sadly, that this could be NaN because I also failed to pack the spike counts into
                    % spike_counts_allheadings_grouped if the number of left/right choices were fewer than 3... WTF...

                    zero_spikes_left_right = [mean(this_raw.spike_counts_allheadings_grouped{zero_index, 1}); mean(this_raw.spike_counts_allheadings_grouped{zero_index, 2})];
                    this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); zero_spikes_left_right; this_tuning_correctonly(zero_index+1:end)];
 
%                     this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); this_tuning_correctonly(zero_index+1:end)]; % No zero version
                    
                    % Align preferred direcions to 'rightward'. So if PREF == 1 for a certain cell, we
                    % flip all its tuning curve.
                    if flip_tuning
                        this_tuning = flipud(this_tuning);
                        this_tuning_correctonly = flipud(this_tuning_correctonly);
                    end
                    
                    % Undate dynamic range
                    if sum(pp == time_for_dynamic_range) > 0
                        dynamic_min = min(dynamic_min, min(this_tuning_correctonly(:)));
                        dynamic_max = max(dynamic_max, max(this_tuning_correctonly(:)));
                        
                        dynamic_min_own_range(k) = min(dynamic_min_own_range(k), min(this_tuning_correctonly(:)));
                        dynamic_max_own_range(k) = max(dynamic_min_own_range(k), max(this_tuning_correctonly(:)));
                    end
                    
                    % Pack data
                    try
                        tuning_pack{1}{k,pp}(:,i) = this_tuning';
                        tuning_pack{2}{k,pp}(:,i) = this_tuning_correctonly';
                    catch
                        tuning_pack{1}{k,pp}(:,i) = nan;
                        tuning_pack{2}{k,pp}(:,i) = nan;
                    end
                end
            end
            
            % Normalization according to dynamic range of each neuron
            % (lowest in all tuning curve = 0, highest = 1)
            offset = dynamic_min;
            gain = dynamic_max - dynamic_min;
            
            offset_own_range = dynamic_min_own_range;
            gain_own_range = dynamic_max_own_range - dynamic_min_own_range;
            
            for pp = 1:length(CP_ts{j})
                for k = 1:3
                    if modalities_share_same_dynamic_range
                        % Three modalities share the same dynamic range
                        tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) - offset;
                        tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) / gain;
                        
                        tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) - offset;
                        tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) / gain;
                    else
                        % Three modalities have their own dynamic ranges
                        tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) - offset_own_range(k);
                        tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) / gain_own_range(k);
                        
                        tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) - offset_own_range(k);
                        tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) / gain_own_range(k);
                    end
                end
            end
            
        end
        
        % Get means and sems
        tuning_mean_all = nan(length(unique_heading),length(CP_ts{j}),3);
        tuning_sem_all = tuning_mean_all;
        tuning_mean_correctonly = nan(length(unique_heading_two_zeros),length(CP_ts{j}),3);
        tuning_sem_correctonly = tuning_mean_correctonly;
        tuning_sig_fit_correctonly = nan(3,length(CP_ts{j}),3);
        
        % For sigmoid fitting. HH20150528
%         sigfunc = @(sig_para, x)(sig_para(1)+ sig_para(2)./ (1 + exp(-x/sig_para(3))));
        sigfunc = @(sig_para, x)(sig_para(1)+ sig_para(2) * normcdf(x,0,sig_para(3)));
        
        function err = cost_function(q,data_cum)
            x = data_cum(:,1);
            y = data_cum(:,2);
            
            z = sigfunc(q,x);
            err = norm(z-y);
        end

        for pp = 1:length(CP_ts{j})
            
            for k = 1:3
                % Mean and sem
                % Patch calculation for LEFT & RIGHT choices of 0 headings. HH20160214
                this_tuning_all = tuning_pack{1}{k,pp}(:,~isnan(tuning_pack{1}{k,pp}(1,:)));
                tuning_mean_all(:,pp,k) = nanmean(this_tuning_all,2);
                tuning_sem_all(:,pp,k) = nanstd(this_tuning_all,[],2)./sqrt(sum(~isnan(this_tuning_all),2));
                
                this_tuning_correctonly =  tuning_pack{2}{k,pp}(:,~isnan(tuning_pack{2}{k,pp}(1,:)));
                tuning_mean_correctonly(:,pp,k) = nanmean(this_tuning_correctonly,2);
                tuning_sem_correctonly(:,pp,k) = nanstd(this_tuning_correctonly,[],2)./sqrt(sum(~isnan(this_tuning_correctonly),2));;
                
                % Fitting sigmoid function. HH20150528
                yy = nanmean(this_tuning_correctonly,2);
%                 tuning_sig_fit_correctonly(:,pp,k) = nlinfit(unique_heading, yy, sigfunc, [min(yy) range(yy) 1]);

                % Begin optimization
                quick = fminsearch(@(q)cost_function(q,[unique_heading_two_zeros yy]),[min(yy) range(yy) 1]);
                
                % Output
                tuning_sig_fit_correctonly(:,pp,k) = quick;
               
            end
            
        end
        
        %%
        
        % Plotting fitting slopes
        slope = squeeze(tuning_sig_fit_correctonly(3,:,:));
        slope(abs(slope)>200) = nan;
        figure(); plot(CP_ts{j},slope,'Linew',2); SetFigure();
        
        
        % Plotting tuning curves at three time points
        set(figure(2099+figN),'pos',[9 206 1210 754],'name',['Tuning curve, j = ' num2str(j)]); clf; figN = figN+1;
        
        
        for pp = 1:length(tuning_time_phase)
            
            for k = 1:3  % For each stim type
                
                % Plotting
                subplot(2,length(tuning_time_phase),pp); hold on; ylabel('All');
                plot(unique_heading,tuning_mean_all(:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
                h = errorbar(unique_heading,tuning_mean_all(:,tuning_time_phase(pp),k),tuning_sem_all(:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2);
                errorbar_tick(h,10000);
                title([tuning_time_phase_title{pp} ', n = ' num2str(size(this_tuning_all,2))]);
                axis tight; xlim(xlim*1.1);
                
                subplot(2,length(tuning_time_phase),pp + length(tuning_time_phase));  hold on; ylabel('Correct only');
                plot(unique_heading_two_zeros,tuning_mean_correctonly(:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
                h = errorbar(unique_heading_two_zeros,tuning_mean_correctonly(:,tuning_time_phase(pp),k),tuning_sem_correctonly(:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2);
                try errorbar_tick(h,10000); catch end;
                title([tuning_time_phase_title{pp} ', t = ' num2str(CP_ts{j}(tuning_time_phase(pp))),...
                    ' n = ' num2str(size(this_tuning_correctonly,2))]);
                axis tight; xlim(xlim*1.1);
                
                SetFigure(15);
                
                % Plot sigmoid fitting
                xx = linspace(min(unique_heading),max(unique_heading),100);
                plot(xx,sigfunc(tuning_sig_fit_correctonly(:,tuning_time_phase(pp),k),xx),'color',colors(k,:),'linew',3);
                
            end
        end
        
        %% Tuning: Animation
        
                % Plotting tuning curves at three time points
                set(figure(2099+figN),'pos',[91 15 687 562],'name',['Tuning curve, j = ' num2str(j)]); clf; figN = figN+1;
                ylabel('Correct only');
                SetFigure(15); drawnow;
        
                % Preallocate movie structure.
                mov_tuning(1:length(CP_ts{j})) = struct('cdata', [], 'colormap', []);
        
                for pp = 1:length(CP_ts{j})
                    hold off;
                    for k = 1:3  % For each stim type
                        % Plotting
                        plot(unique_heading_two_zeros,tuning_mean_correctonly(:,pp,k),'o','markersize',9,'color',colors(k,:),'LineWid',1);
                        hold on;
                        h = errorbar(unique_heading_two_zeros,tuning_mean_correctonly(:,pp,k),tuning_sem_correctonly(:,pp,k),'color',colors(k,:),'LineWid',1);
                        try errorbar_tick(h,10000); catch end
                        
                        % Plot sigmoid fitting
                        xx = linspace(min(unique_heading),max(unique_heading),100);
                        plot(xx,sigfunc(tuning_sig_fit_correctonly(:,pp,k),xx),'color',colors(k,:),'linew',3);

                    end
                    axis([-7 7 0.25 0.8]);
                    title(sprintf('t = %g, \\sigma = %g %g %g',CP_ts{j}(pp),squeeze(tuning_sig_fit_correctonly(3,pp,:))));
                    axis tight;
        
                    mov_tuning(pp) = getframe(gcf);
        
                    drawnow;
                    pause(0.02);
                end
                
                try
                    movie2avi(mov_tuning,'LIP_tuning_evolve_j_1_3.avi','compression','none');
                catch
                    disp('==== Animation not saved ====');
                end
        
        %% Tuning: Hotgram
        to_evolve = {tuning_mean_all tuning_mean_correctonly};
        to_evolve_title = {'All' 'Correct only'};
        unique_heading_in_use = {unique_heading unique_heading_two_zeros};
        
        for co = 1:2
            set(figure(2099+figN),'pos',[82 99 927 855],'name','Evolve of tuning curves'); clf; figN = figN+1;
            
            [X,Y] = meshgrid(CP_ts{j},unique_heading_in_use{co});
            [Xq,Yq] = meshgrid(linspace(min(CP_ts{j}),max(CP_ts{j}),1000),linspace(min(unique_heading_in_use{co}),max(unique_heading_in_use{co}),100));
            
            
            for k = 1:3  % For each stim type
                subplot_tight(1,3,k,0.02,0.1);
%                 surf(Xq,Yq,interp2(X,Y,to_evolve{co}(:,:,k),Xq,Yq),'edgecolor','none'); view(90,-90); hold on;
                contourf(Xq,Yq,interp2(X,Y,to_evolve{co}(:,:,k),Xq,Yq),40,'edgecolor','k'); view(90,-90); hold on;
                
                if k>1 ;set(gca,'xtick',[]); end
                if k==2 ; title(to_evolve_title{co}); end
                
                caxis([min(to_evolve{co}(:)) 1.05*max(to_evolve{co}(:))]); % Use the same color range
                xlim([min(CP_ts{j}) max(CP_ts{j})]); ylim([min(unique_heading_in_use{co}) max(unique_heading_in_use{co})]);
                
                for tt = 1:3
                    plot3([1 1] * time_markers{j}(1,tt),ylim,-3*[1 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
                end
                
                load('colormap_for_tuning_hotgram','colormap_for_tuning_hotgram');
                colormap(colormap_for_tuning_hotgram);
                
                %                 axis tight;
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
                
            end
            
            SetFigure();
        end
    end

    marker_size = 11;
    tcell_cross_size = 15;

    function f3p1(debug)      % Correlations 1. Mem-sac v.s. CD / CP
        if debug
            dbstack;
            keyboard;
        end
        
        % mean(A(:,2:3),2)+0.4,mean(ChoiceDiv_All(selectCells,rate_ts > -500 & rate_ts < 00,1),2),'o'
        
        % Maybe I should calculate a corresponding ChoiceDiv for MemSac at PREF/NULL target locations??  @HH20150419
        
        j = 2;  % Not much time-sensitive, so I use j = 2 here.
        
        nHist = 12;
        
        % ------- Pre-sac CDiv vs. MemSac -------
        
        LinearCorrelation(MemSac_indicator(select_bottom_line),...
            [  mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),...
            mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),...
            mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,3),2)],...
            'Xlabel',MemSac_indicator_txt,'Ylabel','Pre-sac CDiv','FaceColors',{'b','r','g'},'Markers',{'o'},...
            'LineStyles',{'b-','r-','g-'},'MarkerSize',marker_size,'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
        set(gcf,'name',['j = ' num2str(j)]);
        
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(MemSac_indicator(select_tcells),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(MemSac_indicator(select_tcells),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,3),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        SetFigure(20);
        
        % ------- Post-sac CDiv vs. MemSac -------
        
        LinearCorrelation(MemSac_indicator(select_bottom_line),...
            [  mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > 100 & rate_ts{j} < 500,1),2),...
            mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > 100 & rate_ts{j} < 500,2),2),...
            mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > 100 & rate_ts{j} < 500,3),2)],...
            'Xlabel',MemSac_indicator_txt,'Ylabel','Post-sac CDiv','FaceColors',{'b','r','g'},'Markers',{'o'},...
            'LineStyles',{'b-','r-','g-'},'MarkerSize',marker_size,'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
        set(gcf,'name',['j = ' num2str(j)]);
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > 100 & rate_ts{j} < 500,1),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(MemSac_indicator(select_tcells),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > 100 & rate_ts{j} < 500,2),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(MemSac_indicator(select_tcells),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > 100 & rate_ts{j} < 500,3),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        SetFigure(20);
        
        % LinearCorrelation(mean(MemSac_PREFmNULL_PSTH(select_bottom_line,memsac_ts > -300 & memsac_ts < -100),2),...
        %     [  mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),...
        %        mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),...
        %        mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,3),2)],...
        %         'Xlabel','MemSac Rate Diff (pre)','Ylabel','Pre-sac Cdiv','FaceColors',{'b','r','g'},'Markers',{'o'},...
        %         'LineStyles',{'b-','r-','g-'},'MarkerSize',12,'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
        
        % =================== Correlations 2. Mem-sac v.s. CP
        
        j = 2; % Not much time-sensitive, so I use j = 2 here.
        CP_interest = -500;
        
        CP_interval = CP_ts{j} == CP_interest ;
        
        for k = 1:3
            CP_NS{k} = find_bottom_line(CP_p{j}(select_bottom_line,CP_interval,k) > 0.05);
            CP_S{k} = find_bottom_line(CP_p{j}(select_bottom_line,CP_interval,k) <= 0.05);
        end
        
        
        h = LinearCorrelation(...
            {   MemSac_indicator(CP_NS{1}),...
            MemSac_indicator(CP_S{1}),...
            MemSac_indicator(CP_NS{2}),...
            MemSac_indicator(CP_S{2}),...
            MemSac_indicator(CP_NS{3}),...
            MemSac_indicator(CP_S{3}),...
            },...
            {  mean(CP{j}(CP_NS{1},CP_interval,1),2),...
            mean(CP{j}(CP_S{1},CP_interval,1),2),...
            mean(CP{j}(CP_NS{2},CP_interval,2),2),...
            mean(CP{j}(CP_S{2},CP_interval,2),2),...
            mean(CP{j}(CP_NS{3},CP_interval,3),2),...
            mean(CP{j}(CP_S{3},CP_interval,3),2)},...
            'CombinedIndex',[3 12 48],...
            'Xlabel',MemSac_indicator_txt,'Ylabel',['CP at Sac ' num2str(CP_interest)],...
            'FaceColors',{'none','b','none','r','none','g'},'Markers',{'o'},...
            'LineStyles',{'b:','b:','r:','r:','g:','g:','b-','r-','g-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
        set(gcf,'name',['j = ' num2str(j)]);
        delete([h.group(1:6).line]);
        plot(xlim,[0.5 0.5],'k--');         SetFigure(20);
        
        
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),mean(CP{j}(select_tcells,CP_interval,1),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(MemSac_indicator(select_tcells),mean(CP{j}(select_tcells,CP_interval,2),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(MemSac_indicator(select_tcells),mean(CP{j}(select_tcells,CP_interval,3),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        h = LinearCorrelation({  mean(ChoiceDiv_All{j}(CP_NS{1},rate_ts{j} > -400 & rate_ts{j} < -100,1),2)+0.5,...
            mean(ChoiceDiv_All{j}(CP_S{1},rate_ts{j} > -400 & rate_ts{j} < -100,1),2)+0.5,...
            mean(ChoiceDiv_All{j}(CP_NS{2},rate_ts{j} > -400 & rate_ts{j} < -100,2),2)+0.5,...
            mean(ChoiceDiv_All{j}(CP_S{2},rate_ts{j} > -400 & rate_ts{j} < -100,2),2)+0.5,...
            mean(ChoiceDiv_All{j}(CP_NS{3},rate_ts{j} > -400 & rate_ts{j} < -100,3),2)+0.5,...
            mean(ChoiceDiv_All{j}(CP_S{3},rate_ts{j} > -400 & rate_ts{j} < -100,3),2)+0.5},...
            {  mean(CP{j}(CP_NS{1},CP_interval,1),2),...
            mean(CP{j}(CP_S{1},CP_interval,1),2),...
            mean(CP{j}(CP_NS{2},CP_interval,2),2),...
            mean(CP{j}(CP_S{2},CP_interval,2),2),...
            mean(CP{j}(CP_NS{3},CP_interval,3),2),...
            mean(CP{j}(CP_S{3},CP_interval,3),2)},...
            'CombinedIndex',[3 12 48],...
            'Xlabel','Choice Div (pre) + 0.5','Ylabel',['CP at Sac ' num2str(CP_interest)],...
            'FaceColors',{'none','b','none','r','none','g'},'Markers',{'o'},...
            'LineStyles',{'b:','b:','r:','r:','g:','g:','b-','r-','g-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',nHist,'YHist',nHist,'SameScale',1); figN = figN + 1;
        delete([h.group(1:6).line]);
        plot(xlim,[0.5 0.5],'k--'); plot([0.5 0.5],ylim,'k--');         SetFigure(20);
        
        
        set(gcf,'name',['j = ' num2str(j)]);
        % Annotate tcells
        plot(.5+mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),mean(CP{j}(select_tcells,CP_interval,1),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(.5+mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),mean(CP{j}(select_tcells,CP_interval,2),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(.5+mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,3),2),mean(CP{j}(select_tcells,CP_interval,3),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        
        % ------------ Memsac DDI and Choice Preference
        
        k = 1;
        tt = 1;
        cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
        
        MemSac_DDI_selected = MemSac_indicator(select_bottom_line);
        
        h = LinearCorrelation({
            MemSac_DDI_selected((~cpref_sig));
            MemSac_DDI_selected((cpref_sig));
            MemSac_DDI_selected((~cpref_sig));
            MemSac_DDI_selected((cpref_sig));
            MemSac_DDI_selected((~cpref_sig));
            MemSac_DDI_selected((cpref_sig));
            },...
            {
            abs(Choice_pref_all(1,~cpref_sig,tt)) ;
            abs(Choice_pref_all(1,cpref_sig,tt)) ;
            abs(Choice_pref_all(2,~cpref_sig,tt)) ;
            abs(Choice_pref_all(2,cpref_sig,tt)) ;
            abs(Choice_pref_all(3,~cpref_sig,tt)) ;
            abs(Choice_pref_all(3,cpref_sig,tt)) ;
            },...
            'CombinedIndex',[3 12 48],...
            'Xlabel',MemSac_indicator_txt,'Ylabel','abs(Choice preference)',...
            'FaceColors',{'none','b','none','r','none','g'},'Markers',{'o'},...
            'LineStyles',{'b:','b:','r:','r:','g:','g:','b-','r-','g-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','grouped','SameScale',0,'Method','Spearman'); figN = figN + 1;

% h = LinearCorrelation({
%             MemSac_DDI_selected;
%             MemSac_DDI_selected
%             MemSac_DDI_selected;
%             },...
%             {
%             abs(Choice_pref_all(1,:,tt)) ;
%             abs(Choice_pref_all(2,:,tt)) ;
%             abs(Choice_pref_all(3,:,tt)) ;
%             },...
%             'CombinedIndex',[],...
%             'Xlabel',MemSac_indicator_txt,'Ylabel','abs(Choice preference)',...
%             'FaceColors',{'b','r','g'},'Markers',{'o'},...
%             'LineStyles',{'b-','r-','g-'},'MarkerSize',marker_size,...
%             'figN',figN,'XHist',20,'YHist',20,...
%             'XHistStyle','stacked','YHistStyle','grouped','SameScale',0,'Method','Spearman'); figN = figN + 1;
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        delete([h.group(1:6).line]);
        % Annotate tcells
        plot(MemSac_DDI_selected(select_tcells(select_bottom_line)),abs(Choice_pref_all(1,select_tcells(select_bottom_line),tt)),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(MemSac_DDI_selected(select_tcells(select_bottom_line)),abs(Choice_pref_all(2,select_tcells(select_bottom_line),tt)),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(MemSac_DDI_selected(select_tcells(select_bottom_line)),abs(Choice_pref_all(3,select_tcells(select_bottom_line),tt)),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(MemSac_DDI_selected(:),abs(Choice_pref_all(k,:,tt)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_bottom_line});
        
        
    end
    function f3p2(debug)      % Correlations 2. Choice preference v.s. modality preference (Anne Fig.2f-h; Fig.3a)
        if debug
            dbstack;
            keyboard;
        end
        
        %% ---------------- 1 Choice pref v.s. modality pref (Anne Fig.3a)-------------------
        tt = 1;
        
        for k = 1:3
            set(figure(figN),'name','CDiv vs. MDiv','pos',[17 514 1151 449]);
            cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
            mpref_sig = Modality_pref_p_value_all(1,:,tt) < 0.05;
            
            h = LinearCorrelation({
                (Modality_pref_all(1, ~cpref_sig & ~mpref_sig,tt)) ;
                (Modality_pref_all(1, xor(cpref_sig , mpref_sig),tt));
                (Modality_pref_all(1, cpref_sig & mpref_sig,tt))},...
                {
                (Choice_pref_all(k,~cpref_sig & ~mpref_sig,tt)) ;
                (Choice_pref_all(k,xor(cpref_sig , mpref_sig),tt)) ;
                (Choice_pref_all(k,cpref_sig & mpref_sig,tt)) },...
                'CombinedIndex',[7],...
                'Ylabel','Choice preference (pre)','Xlabel','Modality preference (pre)',...
                'FaceColors',{'none',colors(k,:)*0.2 + [0.8 0.8 0.8],colors(k,:)},'Markers',{'o'},...
                'LineStyles',{'k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,'Method','Spearman'); figN = figN + 1;
            delete([h.group(1:3).line h.diag]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       SetFigure(15);
            axis([-1 1 -1 1])
            set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1);
        
           
            % Annotate tcells
            plot(Modality_pref_all(1,select_tcells(select_bottom_line),tt),...
                Choice_pref_all(k,select_tcells(select_bottom_line),tt),'+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(Modality_pref_all(1,:,tt),Choice_pref_all(k,:,tt),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
            
            %     for i = 1:3
            %         set(h.group(i).dots,'color',colors(k,:));
            %     end
        end
        
        %% ---------------- 2 Choice pref (vest and visual) -------------------
        
        tt = 1;
        
        set(figure(figN),'name','CDiv (visual) vs. CDiv (vest)','pos',[17 514 1151 449]);
        
        cpref_sig_1 = Choice_pref_p_value_all(1,:,tt) < 0.05;
        cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
        
        % Abs() or not?
        Choice_pref_all_temp = (Choice_pref_all);
        
        h = LinearCorrelation({
            (Choice_pref_all_temp(2,~cpref_sig_1 & ~cpref_sig_2,tt)) ;
            (Choice_pref_all_temp(2, xor(cpref_sig_1,cpref_sig_2),tt)) ;
            (Choice_pref_all_temp(2,cpref_sig_1 & cpref_sig_2,tt))
            },...
            {
            (Choice_pref_all_temp(1,~cpref_sig_1 & ~cpref_sig_2,tt)) ;
            (Choice_pref_all_temp(1,xor(cpref_sig_1,cpref_sig_2),tt)) ;
            (Choice_pref_all_temp(1,cpref_sig_1 & cpref_sig_2,tt))
            },...
            'CombinedIndex',[7],...
            'Ylabel','Vestibular choice preference','Xlabel','Visual choice preference',...
            'FaceColors',{'none',[0.8 0.8 0.8],'k'},'Markers',{'o'},...
            'LineStyles',{'k:','k:','k:','k-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,'Method','Spearman'); figN = figN + 1;
        delete([h.group(1:3).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        
        % Annotate tcells
        plot(Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt),Choice_pref_all_temp(1,select_tcells(select_bottom_line),tt),...
            '+','markersize',tcell_cross_size,'color','r','linew',2);
        
        axis([-1 1 -1 1])
        set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(Choice_pref_all(2,:,tt),Choice_pref_all(1,:,tt),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
        
        clear Choice_pref_all_temp;
        
        %     for i = 1:3
        %         set(h.group(i).dots,'color',colors(k,:));
        %     end
        
        %% ---------------- 3 Choice pref (single vs comb) -------------------
        
        tt = 1;
        
        set(figure(figN),'name','CDiv (single) vs. CDiv (comb)','pos',[17 514 1151 449]);
        
        cpref_sig_1 = Choice_pref_p_value_all(1,:,tt) < 0.05;
        cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
        cpref_sig_3 = Choice_pref_p_value_all(3,:,tt) < 0.05;
        
        % Abs() or not?
        Choice_pref_all_temp = (Choice_pref_all);
        
        h = LinearCorrelation({
            (Choice_pref_all_temp(1,~cpref_sig_1 & ~cpref_sig_3,tt)) ;
            (Choice_pref_all_temp(1,xor(cpref_sig_1 ,cpref_sig_3),tt)) ;
            (Choice_pref_all_temp(1,cpref_sig_1 & cpref_sig_3,tt)) ;
            (Choice_pref_all_temp(2,~cpref_sig_2 & ~cpref_sig_3,tt)) ;
            (Choice_pref_all_temp(2,xor(cpref_sig_2 ,cpref_sig_3),tt)) ;
            (Choice_pref_all_temp(2,cpref_sig_2 & cpref_sig_3,tt))
            
            },...
            {
            (Choice_pref_all_temp(3,~cpref_sig_1 & ~cpref_sig_3,tt)) ;
            (Choice_pref_all_temp(3, xor(cpref_sig_1 , cpref_sig_3),tt)) ;
            (Choice_pref_all_temp(3,cpref_sig_1 & cpref_sig_3,tt)) ;
            (Choice_pref_all_temp(3,~cpref_sig_2 & ~cpref_sig_3,tt)) ;
            (Choice_pref_all_temp(3, xor(cpref_sig_2 , cpref_sig_3),tt)) ;
            (Choice_pref_all_temp(3,cpref_sig_2 & cpref_sig_3,tt))
            
            },...
            'CombinedIndex',[7 56],...
            'Xlabel','Single modality choice preference','Ylabel','Combined choice preference',...
            'FaceColors',{'none',[0.8 0.8 1],[0.2 0.2 1],'none',[1 0.8 0.8],'r'},'Markers',{'o'},...
            'LineStyles',{'k:','k:','k:','k:','k:','k:','b-','r-'},'MarkerSize',marker_size,...
            'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
            'SameScale',1,'Method','Spearman'); figN = figN + 1;
        
        delete([h.group(1:6).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        axis([-1 1 -1 1]); set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1); axis square;
        
        % Annotate tcells
        plot((Choice_pref_all_temp(1,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot([Choice_pref_all_temp(1,:,tt) Choice_pref_all_temp(2,:,tt)],[Choice_pref_all_temp(3,:,tt) Choice_pref_all_temp(3,:,tt)],'visible','off'); hold on;
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_cpref_mpref; select_cpref_mpref],2});
        
        clear Choice_pref_all_temp;
        
        %             for i = 1:3
        %                 set(h.group(i).dots,'color',colors(k,:));
        %             end
        
        %% ---------------- 4 Choice pref (pre vs post) -------------------
        
        set(figure(figN),'name','CDiv (pre-sac) vs. CDiv (post-sac)','pos',[17 514 1151 449]);
        
        k = 3;
        
        cpref_sig_pre = Choice_pref_p_value_all(k,:,1) < 0.05;
        cpref_sig_post = Choice_pref_p_value_all(k,:,2) < 0.05;
        
        h = LinearCorrelation({
            (Choice_pref_all(k,~cpref_sig_pre & ~cpref_sig_post,1)) ;
            (Choice_pref_all(k, xor(cpref_sig_pre , cpref_sig_post),1)) ;
            (Choice_pref_all(k,cpref_sig_pre & cpref_sig_post,1)) ;
            },...
            {
            (Choice_pref_all(k,~cpref_sig_pre & ~cpref_sig_post,2)) ;
            (Choice_pref_all(k, xor(cpref_sig_pre , cpref_sig_post),2)) ;
            (Choice_pref_all(k,cpref_sig_pre & cpref_sig_post,2)) ;
            
            },...
            'CombinedIndex',[7],...
            'Xlabel','Choice preference: decision formation','Ylabel','Choice preference: movement',...
            'FaceColors',{'none',colors(k,:)*0.2 + [0.8 0.8 0.8],colors(k,:)},'Markers',{'o'},...
            'LineStyles',{'k:','k:','k:','k-'},'MarkerSize',marker_size,...
            'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
            'SameScale',1,'Method','Spearman'); figN = figN + 1;         SetFigure(20);
        
        
        delete([h.group(1:3).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
        axis([-1 1 -1 1]); set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1); axis square;
        
        % Annotate tcells
        plot(Choice_pref_all(k,select_tcells(select_bottom_line),1),Choice_pref_all(k,select_tcells(select_bottom_line),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(Choice_pref_all(k,:,1),Choice_pref_all(k,:,2),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
        
        %     for i = 1:3
        %         set(h.group(i).dots,'color',colors(k,:));
        %     end
        
    end
    function f3p3(debug)      % Correlations 3. Psychophysics v.s. CD
        if debug
            dbstack;
            keyboard;
        end
        %%
        % T-cell is important here because Polo's behavior was getting worse while
        % I was getting better at finding T-cells, so...
        select_psycho = select_bottom_line & select_tcells; %& ~[zeros(80,1);ones(138-80,1)];
        
        j = 1;
        
        % Enhancement of cDiv in combined condition
        % enhance_cdiv = max(ChoiceDiv_ModDiffer{1}(:,0 <= rate_ts{j} & rate_ts{j} <= 1500,2),[],2); % Comb - vis
        
        t_begin = 700; t_end = 800;
        enhance_cdiv = nanmean(ChoiceDiv_ModDiffer{1}(:,700 <= rate_ts{j} & rate_ts{j} <= 800, 3 ),2); % Comb - vis
        
        h = LinearCorrelation({
            Psy_pred_ratio(select_psycho)
            Psy_pred_ratio(select_psycho)
            Psy_pred_ratio(select_psycho)},...
            {
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 1 ),2);
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 2 ),2);
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 3 ),2);
            },...
            'FaceColors',{'b','r','k'},'Markers',{'o'},...
            'LineStyles',{'b-','r-','k-'},'MarkerSize',marker_size,...
            'Ylabel',sprintf('\\Delta CDiv (%g - %g ms, 3-1, 3-2, 1-2)',t_begin,t_end),'Xlabel','Psycho prediction ratio',...
            'MarkerSize',12,...
            'figN',figN,'XHist',15,'YHist',15,'logx',1,...
            'XHistStyle','stacked','YHistStyle','stacked','Method','Pearson'); figN = figN + 1;
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        
    end

    function f4p0(debug)      % Pack PCA_A
        if debug
            dbstack;
            keyboard;
        end
        %%
        % Pack all data into matrix A: [memsacDDI, vest_div, vis_div, comb_div]
        
        % Now it becomes A: [memsacDDI, vest_div, vis_div, comb_div, 2-1 mod_div, 3-1 mod_div, 3-2 mod_div]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sort_time_interval1 = [find(rate_ts{j_PCA_A} > time_markers{j_PCA_A}(1,1), 1)  find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,3),1)-1]; % Stim on to saccade on
        sort_time_interval2 = [find(rate_ts{j_PCA_A} > time_markers{j_PCA_A}(1,3),1)  length(rate_ts{j_PCA_A})]; % Post sac
        sort_time_interval3 = [1  length(CP_ts{j_PCA_A})];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Exclude any cells that have NaNs
        % for k = 1:3
        %     selectCells = selectCells & ~isnan(ChoiceDiv_All(:,1,k));
        % end
        
        
        % --- Sort A according to different methods for visualizing ---
        
        A_CP = reshape(CP{j_PCA_A}(select_bottom_line,:,:),sum(select_bottom_line),[])-0.5;
        A_choicediv= reshape(ChoiceDiv_All{j_PCA_A}(select_bottom_line,:,:),sum(select_bottom_line),[]);
        
        % A_memSac = xls_num{1}(selectCells,header.DDI_LS:header.DDI_post) - 0.4; % All memSac data
        % LS, M, Pre, Co, Post
        A_memSac = MemSac_DDI(select_bottom_line,2:end) - 0.4;  % Read memsac from .mat file. HH20150413
        
        
        % Now I add mod_div. HH20150415
        A_moddiv = reshape(ModDiv_All{j_PCA_A}(select_bottom_line,:,:),sum(select_bottom_line),[]);
        
        PCA_A = [A_memSac A_choicediv A_moddiv]; % A_CP];
        
    end
    function f4p1(debug)      % 1. Hot-gram
        if debug
            dbstack;
            keyboard;
        end
        
        % If data not packed
        if isempty(PCA_A)
            f4p0(0);
        end
        
        sort_method = {
            [1 4];
            %     [3 4];
            
            %     5 + 0*length(rate_ts{j_PCA_A}) + sort_time_interval1;
            %     5 + 1*length(rate_ts{j_PCA_A}) + sort_time_interval1;
            5 + 2*length(rate_ts{j_PCA_A}) + sort_time_interval1;
            
            %     5 + 0*length(rate_ts{j_PCA_A}) + sort_time_interval2;
            %     5 + 1*length(rate_ts{j_PCA_A}) + sort_time_interval2;
            5 + 2*length(rate_ts{j_PCA_A}) + sort_time_interval2;
            
            5 + 3* length(rate_ts{j_PCA_A}) + sort_time_interval1;
            %     5 + 4* length(rate_ts{j_PCA_A}) + sort_time_interval1;
            %     5 + 5* length(rate_ts{j_PCA_A}) + sort_time_interval1;
            
            
            %     5 + 6* length(rate_ts{j_PCA_A}) + 0*length(CP_ts{j_PCA_A}) + sort_time_interval3;
            %     5 + 6* length(rate_ts{j_PCA_A}) + 1*length(CP_ts{j_PCA_A}) + sort_time_interval3;
            %     5 + 6* length(rate_ts{j_PCA_A}) + 2*length(CP_ts{j_PCA_A}) + sort_time_interval3;
            };
        
        for sort_ind = 1:length(sort_method)
            
            sort_begin = sort_method{sort_ind}(1);
            sort_end = sort_method{sort_ind}(end);
            
            mean_for_sort = mean(PCA_A(:,sort_begin:sort_end),2);
            mean_for_sort(isnan(mean_for_sort)) = -inf;  % Nan goes last
            
            [~,sort_order] = sort(mean_for_sort);
            
            % Enlarge mem-sac part for clarity
            A_forplot = [reshape(repmat(A_memSac,enlarge_factor,1),size(A_memSac,1),[]) A_choicediv A_moddiv A_moddiv(:,end)];% A_CP A_CP(:,end)]; % To ensure the last value to be shown
            A_forplot = A_forplot(sort_order,:);
            
            A_forplot = [A_forplot; A_forplot(end,:)]; % To ensure the last cell to be shown
            
            A_forplot(isnan(A_forplot)) = -0.5; % This workaround is to avoid that cells which are next the NaNs cannot be displayed.
            
            [nn, tt] = meshgrid(1:size(A_forplot,1), 1:size(A_forplot,2));
            
            set(figure(2099+figN),'name',['CDiv and MDiv, hot-gram, j_PCA_A = ' num2str(j_PCA_A)],'pos',[409 182 913 775]); clf; figN = figN+1;
            h1 = surf(tt,nn,A_forplot','Edgecolor','none');
            view(0,90);
            
            % Annotate time separations
            hold on;
            CD_begin_at = enlarge_factor * 5 + 1;
            
            
            % Stim on / stim off / sac on
            for ttt = 1:3
                for tt = 0:5
                    plot3(CD_begin_at + tt*length(rate_ts{j_PCA_A}) + find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,ttt),1) * [1 1],...
                        [1 size(A_forplot,1)],[1 1],'k','linesty',marker_for_time_markers{j_PCA_A}{ttt},'linewid',1.5);
                end
            end
            
            % End
            plot3([CD_begin_at CD_begin_at],[1 size(A_forplot,1)],[1 1],'k','linewid',3);
            for tt = 1:6
                plot3(CD_begin_at + tt*[length(rate_ts{j_PCA_A})  length(rate_ts{j_PCA_A})],[1 size(A_forplot,1)],[1 1],'k','linewid',3);
            end
            
            % Annotate tcells
            cell_loc = select_tcells(select_bottom_line);
            h2 = plot3(0,find(cell_loc(sort_order))+.5,1,'+','markersize',5,'color','k','linew',1.5);
            
            % Annotate sort methods
            temp_beg = @(x)((x<=5)*((x-1)*enlarge_factor) + 1 + (x>5)*(CD_begin_at+(x-5)-1));
            temp_end = @(x)((x<=5)*((x)*enlarge_factor) + 1 + (x>5)*(CD_begin_at+(x-5)));
            
            sort_begin_forplot = temp_beg(sort_begin);
            sort_end_forplot = temp_end(sort_end);
            
            plot3([sort_begin_forplot sort_end_forplot],[size(A_forplot,1)+1 size(A_forplot,1)+1],[1 1],'r','linewid',5);
            xlabel('Temporal features');
            ylabel('Cell Number');
            
            set(gca,{'xtick','ytick','ztick'},{[],10:10:size(A_forplot,1),[]},'linewidth',0.00000001);
            axis tight; colorbar; ylim([1 size(PCA_A,1)+2.5]); SetFigure();
            
            % Show individual cell selected from the figure. HH20150424
            [~,where] = sort(sort_order);
            h_line = plot(zeros(1,length(where)),where+.5,'visible','off'); hold on;
            set([gca; h1; h2(:); h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_PCA_B});
            
        end
        
    end
    function f4p2(debug)      % 2. PCA_A (Eigenfeature)
        if debug;  dbstack;  keyboard;  end
        
        % If data not packed
        if isempty(PCA_A)
            f4p0;
        end
        
        % Use the original A instead of A_forplot
        dataPCA = PCA_A;
        
        dataPCA(sum(isnan(dataPCA),2)>0,:) = [];
        
        % [~,sort_order] = sort(mean(dataPCA(:,3:4),2)); % Sort according to mem-sac
        [~,sort_order] = sort(mean(dataPCA(:, 3:4 + 0*length(rate_ts{j_PCA_A}) + sort_time_interval1),2));
        dataPCA = dataPCA(sort_order,:);
        
        % Do PCA
        [PCA_A_PC, ~] = pca(dataPCA);
        
        % The first three eigenvectors that have the first three largest
        % eigenvalues.
        PC1 = PCA_A_PC(:,1);
        PC2 = PCA_A_PC(:,2);
        % PC3 = sortedEigVectors(:,3);
        
        % Projecting the raw data onto the first three eigenvectors
        projPC1 = dataPCA * PC1;
        projPC2 = dataPCA * PC2;
        % projPC3 = dataPCA(:,2:end) * PC3;
        
        set(figure(2099+figN),'name',['Scatter of Eigenfeature, j_PCA_A = ' num2str(j_PCA_A)],'pos',[41 445 543 451]); clf; figN = figN+1;
        
        colorOrder = jet(length(projPC1));
        
        scatter(projPC1, projPC2, 150, colorOrder,'fill');
        grid on; hold on;
        xlabel('PC1'); ylabel('PC2');
        SetFigure();
        
        set(figure(2099+figN),'pos',[602 446 978 450],'name',['Eigenfeature, j_PCA_A = ' num2str(j_PCA_A)]); clf; figN = figN+1;
        
        for ss = 1:2
            subplot(2,1,ss);
            
            eval(['plot(PC' num2str(ss) ',''lineW'',2)']);
            
            hold on;
            
            plot([6 6],ylim,'k','linewid',1);
            plot(xlim,[0 0],'k--');
            
            % Stim on / stim off / sac on
            for ttt = 3
                for tt = 0:5
                    plot(6 + tt*length(rate_ts{j_PCA_A}) + find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,ttt),1) * [1 1],...
                        ylim,'k','linesty',marker_for_time_markers{j_PCA_A}{ttt},'linewid',1.5);
                end
            end
            
            % End
            plot([6 6],ylim,'k','linewid',3);
            for tt = 1:6
                plot(6 + tt*[length(rate_ts{j_PCA_A})  length(rate_ts{j_PCA_A})],ylim,'k','linewid',3);
            end
            
            axis tight; axis off
        end
        
        
        SetFigure();
        
    end

    function f5p0(debug)      % Do PCA_B analysis
        if debug;  dbstack;  keyboard;  end
        
        
        %% Do PCA_B analysis
        % Pack all data (PSTH) into matrix B: [vest_pref, vest_null, vis_pref, vis_null, comb_pref, comb_null, -250_pre, -250_null, 250_pre, 250_null]
        % Only decicion period is included
        
        find_PCA_B = find(select_for_PCA_B);
        
%         PCA_B = nan(sum(select_for_PCA_B),6 * sum(PCA_B_time_range));
        PCA_B = nan(sum(select_for_PCA_B), stim_type_num*2 * sum(PCA_B_time_range));
        
        for i = 1:sum(select_for_PCA_B)
            raw_PSTH = group_result(find_PCA_B(i)).mat_raw_PSTH.PSTH{j_PCA_B,ALL_CHOICE,1}.ys;
%             if size(raw_PSTH,1)==6
            if size(raw_PSTH,1)==stim_type_num*2
                PCA_B(i,:) = reshape(raw_PSTH(:,PCA_B_time_range)',[],1)';
            end
        end
        
        PCA_B(sum(isnan(PCA_B),2) > 0,:) = [];
        
        % Do PCA
        [weights_PCA_B_PC, score, latent, ~, PCA_B_explained] = pca(PCA_B');
        
        % Manually rotate PC1 and PC2 (they're also orthogonal, so there's no
        % change in PCA results). @HH20150424
%         theta = 55 / 180 * pi; % Rotate counterclockwise the two PC bases
%         THETA = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
%         weights_PCA_B_PC(:,1:2) = weights_PCA_B_PC(:,1:2) * THETA;
%         weights_PCA_B_PC(:,2) = -weights_PCA_B_PC(:,2);
        
        % Projecting the raw data onto the first several eigenvectors
        for dim = 1:denoised_dim
%             PCA_B_projPC{dim} = (reshape(weights_PCA_B_PC(:,dim)' * PCA_B,[],6))';
            PCA_B_projPC{dim} = (reshape(weights_PCA_B_PC(:,dim)' * PCA_B,[],stim_type_num*2))';
            %     projPC{dim} = reshape(score(:,dim),[],6)';
        end
        
    end
    function f5p1(debug)      % 1. Weights and correlations  % HH20150413
        if debug;  dbstack;  keyboard;  end
        
        
        % If PCA_B has neven been done, do it first
        if isempty(PCA_B_projPC)
            f5p0(0);
        end
        
        % ------- 3-d PCA weights ---------
        set(figure(2099+figN),'pos',[13 53 836 702],'name','PCA weights in 3d'); clf; figN = figN+1; hold on;
        
        tcell_in_bottom_line = select_tcells(select_for_PCA_B);
        plot3(weights_PCA_B_PC(tcell_in_bottom_line,1),weights_PCA_B_PC(tcell_in_bottom_line,2),weights_PCA_B_PC(tcell_in_bottom_line,3),'r+','markersize',15,'linew',1.5);
        
        h_line = plot3(weights_PCA_B_PC(:,1),weights_PCA_B_PC(:,2),weights_PCA_B_PC(:,3),'ko','markersize',13,'linew',1.5);
        grid off;
        SetFigure(); xlabel('Contribution to Eigen-neuron 1 (choice)'); ylabel('Contribution to Eigen-neuron 2 (modality)'); zlabel('PC3'); grid on; axis tight;
        
        % Show individual cell selected from the figure. HH20150424
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell,h_line,select_for_PCA_B});
        
        % -------- Weights vs cell properties ----------
        Weights_Property_Correlation(weights_PCA_B_PC(:,[1 2]),...
            {'Contribution to Eigen-neuron 1','Contribution to Eigen-neuron 2'},select_for_PCA_B);
        
        %
        %         % ===========   1. PCA weights for naive 2-D space ============
        %         set(figure(2099+figN),'pos',[13 53 836 702],'name',['Variance explained, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        %
        %         tcell_in_bottom_line = select_tcells(select_for_PCA_B);
        %         plot3(PCA_B_PC(tcell_in_bottom_line,1),PCA_B_PC(tcell_in_bottom_line,2),PCA_B_PC(tcell_in_bottom_line,3),'r+','markersize',15,'linew',1.5);
        %
        %         h_line = plot3(PCA_B_PC(:,1),PCA_B_PC(:,2),PCA_B_PC(:,3),'ko','markersize',13,'linew',1.5);
        %         grid off;
        %         SetFigure(); xlabel('Contribution to Eigen-neuron 1 (choice)'); ylabel('Contribution to Eigen-neuron 2 (modality)'); zlabel('PC3'); grid on; axis tight;
        %
        %         % Show individual cell selected from the figure. HH20150424
        %         set([gca h_line],'ButtonDownFcn',{@Show_individual_cell,h_line,select_for_PCA_B});
        %
        %
        %         h = LinearCorrelation({
        %             (PCA_B_PC(:,1));
        %             },...
        %             {
        %             PCA_B_PC(:,2) ;
        %             },...
        %             'CombinedIndex',[],...
        %             'Ylabel','Abs(choice preference)','Xlabel','Contribution to Eigen-neuron 1',...
        %             'FaceColors',{'k'},'Markers',{'o'},...
        %             'LineStyles',{'k-'},'MarkerSize',12,...
        %             'figN',figN,'XHist',20,'YHist',20,...
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %
        %         % Annotate tcells
        %         plot(PCA_B_PC(select_tcells(select_for_PCA_B),1),PCA_B_PC(select_tcells(select_for_PCA_B),2),'+','markersize',16,'color','r','linew',2);
        %
        %         % Show individual cell selected from the figure. HH20150424
        %         set([gca h.group.dots],'ButtonDownFcn',{@Show_individual_cell,h.group.dots,select_for_PCA_B});
        %
        %
        %         % ===========   2. Correlations =============
        %         % ---------  PCA weight 1 and Choice preference
        %         tt = 1;
        %         k = 1;
        %
        %         set(figure(figN),'name','PCA weight 1 vs. Choice div','pos',[17 514 1151 449]);
        %
        %         cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
        %
        %         h = LinearCorrelation({
        %             (PCA_B_PC(~cpref_sig,1));
        %             (PCA_B_PC(cpref_sig,1));
        %             },...
        %             {
        %             abs(Choice_pref_all(k,~cpref_sig ,tt)) ;
        %             abs(Choice_pref_all(k,cpref_sig,tt)) ;
        %             },...
        %             'CombinedIndex',[3],...
        %             'Ylabel','Contribution to Eigen-neuron 2 (modality)','Xlabel','Contribution to Eigen-neuron 1 (choice)',...
        %             'FaceColors',{'none','k'},'Markers',{'o'},...
        %             'LineStyles',{'k-'},'MarkerSize',12,...
        %             'figN',figN,'XHist',20,'YHist',20,...
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
        %         delete([h.group(1:2).line]);
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %
        %         % Annotate tcells
        %         plot(PCA_B_PC(select_tcells(select_for_PCA_B),1),abs(Choice_pref_all(k,select_tcells(select_for_PCA_B),tt)),'+','markersize',16,'color','r','linew',2);
        %
        %         % Show individual cell selected from the figure. HH20150424
        %         h_line = plot(PCA_B_PC(:,1),abs(Choice_pref_all(k,:,tt)),'visible','off'); hold on;
        %         set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
        %
        %         % ---------  Memsac DDI and PCA weight 1
        %
        %         MemSac_DDI_phase = [2:4];
        %
        %         h = LinearCorrelation({
        %             mean(MemSac_DDI(select_for_PCA_B,MemSac_DDI_phase),2);
        %             },...
        %             {
        %             (PCA_B_PC(:,1));
        %             },...
        %             'CombinedIndex',[],...
        %             'Xlabel',['Memsac DDI (' num2str(MemSac_DDI_phase) ')'],'Ylabel','Contribution to Eigen-neuron 1',...
        %             'FaceColors',{'k'},'Markers',{'o'},'Markers',{'o'},...
        %             'LineStyles',{'k-'},'MarkerSize',12,...
        %             'figN',figN,'XHist',20,'YHist',20,...
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
        %
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %
        %         % Annotate tcells
        %         plot(mean(MemSac_DDI(select_tcells,MemSac_DDI_phase),2),PCA_B_PC(select_tcells(select_for_PCA_B),1),...
        %             '+','markersize',16,'color','r','linew',2);
        %
        %         % Show individual cell selected from the figure. HH20150424
        %         h_line = plot(mean(MemSac_DDI(select_for_PCA_B,MemSac_DDI_phase),2),(PCA_B_PC(:,1)),'visible','off'); hold on;
        %         set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_PCA_B});
        %
        %         % ------------ Memsac DDI and Choice Preference
        %
        %         k = 1;
        %         tt = 1;
        %         MemSac_DDI_phase = [2:4];
        %         cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
        %
        %         MemSac_DDI_selected = MemSac_DDI(select_bottom_line,:);
        %         h = LinearCorrelation({
        %             mean(MemSac_DDI_selected((~cpref_sig),MemSac_DDI_phase),2);
        %             mean(MemSac_DDI_selected((cpref_sig),MemSac_DDI_phase),2);
        %             },...
        %             {
        %             abs(Choice_pref_all(k,~cpref_sig,tt)) ;
        %             abs(Choice_pref_all(k,cpref_sig,tt)) ;
        %             },...
        %             'CombinedIndex',[3],...
        %             'Xlabel',['Memsac DDI (' num2str(MemSac_DDI_phase) ')'],'Ylabel','Abs(choice preference)',...
        %             'FaceColors',{'none','k'},'Markers',{'o'},...
        %             'LineStyles',{'k-'},'MarkerSize',12,...
        %             'figN',figN,'XHist',20,'YHist',20,...
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
        %
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %         delete([h.group(1:2).line]);
        %         % Annotate tcells
        %         plot(mean(MemSac_DDI_selected(select_tcells(select_for_PCA_B),MemSac_DDI_phase),2),abs(Choice_pref_all(k,select_tcells(select_for_PCA_B),tt)),...
        %             '+','markersize',16,'color','r','linew',2);
        %
        %         % Show individual cell selected from the figure. HH20150424
        %         h_line = plot(mean(MemSac_DDI_selected(:,MemSac_DDI_phase),2),abs(Choice_pref_all(k,:,tt)),'visible','off'); hold on;
        %         set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_PCA_B});
        %
        %         % ------------ PCA weight 2 and Modality preference
        %         tt = 1;
        %         k = 1; % Vis-vest
        %
        %         set(figure(figN),'name','PCA weight 1 vs. Choice div','pos',[17 514 1151 449]);
        %
        %         mpref_sig = Modality_pref_p_value_all(1,:,tt) < 0.05;
        %
        %         h = LinearCorrelation({
        %             (Modality_pref_all(k,~mpref_sig,tt)) ;
        %             (Modality_pref_all(k,mpref_sig,tt)) ;
        %             },...
        %             {
        %             (PCA_B_PC(~mpref_sig,2));
        %             (PCA_B_PC(mpref_sig,2));
        %             },...
        %             'CombinedIndex',[3],...
        %             'Xlabel','Modality preference (vis - vest)','Ylabel','Contribution to Eigen-neuron 2',...
        %             'FaceColors',{'none','k'},'Markers',{'o'},...
        %             'LineStyles',{'k-'},'MarkerSize',12,...
        %             'figN',figN,'XHist',20,'YHist',20,...
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
        %         delete([h.group(1:2).line]);
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %
        %         % Annotate tcells
        %         plot((Modality_pref_all(k,select_tcells(select_bottom_line),tt)), PCA_B_PC(select_tcells(select_bottom_line),2),...
        %             '+','markersize',16,'color','r','linew',2);
        %
        %         % Show individual cell selected from the figure. HH20150424
        %         h_line = plot((Modality_pref_all(k,:,tt)),PCA_B_PC(:,2),'visible','off'); hold on;
        %         set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_PCA_B});
        
        
    end
    function f5p2(debug)      % 2. 1-D Trajectory
        if debug;  dbstack;  keyboard;  end
        
        % If PCA_B has neven been done, do it first
        if isempty(PCA_B_projPC)
            f5p0(0);
        end
        
        % ============  1. Variance explained =================
        set(figure(2099+figN),'pos',[936 491 733 463],'name',['Variance explained, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        plot((1:length(PCA_B_explained))', cumsum(PCA_B_explained),'o-','markersize',8,'linew',1.5);
        plot((1:denoised_dim)',cumsum(PCA_B_explained(1:denoised_dim)),'ro-','markersize',8,'linew',1.5,'markerfacecol','r');
        plot([0 1],[0 PCA_B_explained(1)],'r-','linew',1.5);
        plot(xlim,[1 1]*sum(PCA_B_explained(1:denoised_dim)),'r--');
        text(denoised_dim,sum(PCA_B_explained(1:denoised_dim))*0.9,[num2str(sum(PCA_B_explained(1:denoised_dim))) '%'],'color','r');
        SetFigure(); xlabel('Num of principal components'); ylabel('Explained variability (%)'); ylim([0 100]);
        
        % ============  2. 1-D Trajectory of Eigen-neurons =================
        
        set(figure(2099+figN),'pos',[462 67 1207 889],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        ds = [1:9];
        
        for d = 1:length(ds)
            ranges(d) = range(PCA_B_projPC{ds(d)}(:));
        end
        
        % Plotting
        for d = 1:length(ds)
            % Normalize to [-1,1]
            gain = max(ranges) / (2 * 0.8) ;
            offset = mean(PCA_B_projPC{ds(d)}(:));
            norm_proj_PC_this = (PCA_B_projPC{ds(d)}-offset)/gain;
            
            h = subplot(fix(sqrt(length(ds))),ceil(length(ds)/fix(sqrt(length(ds)))),d);
            SeriesComparison(shiftdim(norm_proj_PC_this',-1),PCA_B_times,'Colors',{'b','b','r','r',[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'axes',h);
            axis tight; ylim([-1 1]); xlim([-200 1860])
            for tt = 1:3
                plot([1 1] * time_markers{j_PCA_B}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_PCA_B}{tt},'linew',1.5);
            end
            set(gca,'ytick',[-1 0 1]);
            xlabel([]); ylabel('Amplitude (a.u.)'); title(['Eigen-neuron ' num2str(ds(d))]); legend off;
        end
        % set(get(gcf,'children'),'ylim',[min_ylim max_ylim]);
        SetFigure(15);
        
        
        %         figure(904);  clf; hold on;
        %         for k = 1:3
        %             plot(PCA_B_times ,PCA_B_projPC{1}(k*2-1,:)-PCA_B_projPC{1}(k*2,:),'color',colors(k,:),'linew',2);
        %         end
        %         SetFigure;
        
        
    end
    function f5p3(debug)      % 3. 2-D Trajectory  % HH20150413
        if debug;  dbstack;  keyboard;  end
        
        % If PCA_B has neven been done, do it first
        if isempty(PCA_B_projPC)
            f5p0(0);
        end
        
        % ======== 2D ========= %
        set(figure(2099+figN),'pos',[18 170 898 786],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        
        which_two_dimension = [1,2];
        
        for k = 1:stim_type_num
            
            % Time markers
            start_time = 1; % Start point
            % Pref
            h_pref(k) = plot(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,start_time),PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,start_time),'o','color',colors(k,:),'markersize',20,'markerfacecol',colors(k,:));
            % Null
            h_null(k) = plot(PCA_B_projPC{which_two_dimension(1)}(k*2,start_time),PCA_B_projPC{which_two_dimension(2)}(k*2,start_time),'o','color',colors(k,:),'markersize',20,'linew',3);
            
            % Pref
            plot(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,:),PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,:),'-','color',colors(k,:),'linew',3);
            % Null
            plot(PCA_B_projPC{which_two_dimension(1)}(k*2,:),PCA_B_projPC{which_two_dimension(2)}(k*2,:),'--','color',colors(k,:),'linew',3);
            
        end
        
        axis tight;  grid off;
        axis off
        % xlabel('PC1'); ylabel('PC2');
        
        % xlims = xlim; ylims = ylim;
        % h_text = text(xlims(1),ylims(2),'');
        
        
        % ======= Euclidean distance =========
        h_timeaxis = axes('pos', [0.026 0.616 0.344 0.363] ,'color','none');  hold on
        
        prefs = 1:2:2*stim_type_num;
        nulls = 2:2:2*stim_type_num;
        
        % distance = sqrt((PCA_B_projPC{1}(prefs,:) - PCA_B_projPC{1}(nulls,:)).^2 + (PCA_B_projPC{2}(prefs,:) - PCA_B_projPC{2}(nulls,:)).^2 + (PCA_B_projPC{3}(prefs,:) - PCA_B_projPC{3}(nulls,:)).^2);
        % distance = distance';
        
        distance = sum(reshape(cell2mat((cellfun(@(x)(x(prefs,:)-x(nulls,:)).^2',PCA_B_projPC,'uni',0))),[],stim_type_num,denoised_dim),3);
        
        % Normalize distance
        % distance = (distance - repmat(min(distance,[],1),size(distance,1),1))./repmat(range(distance),size(distance,1),1);
        % set(figure(2099+figN),'pos',[18 355 766 601],'name',['Euclidean distance, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        
        for k = 1:stim_type_num
            plot(repmat(PCA_B_times',1,stim_type_num),distance(:,k),'color',colors(k,:),'linew',3);
        end
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_PCA_B}(1),Gauss_vel(:,2)*max(ylim)/4,'--','linew',3,'color',[0.6 0.6 0.6]);
        axis tight
        
        plot(repmat(PCA_B_times',1,stim_type_num),sum(distance(:,[1 2]),2),'k--','linew',2);
        ylim([0 max(distance(:))]);
        
        for tt = 1:3
            plot([1 1] * time_markers{j_PCA_B}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_PCA_B}{tt},'linew',2);
        end
        
        text(mean(xlim)/3,max(ylim)*.8, sprintf('Euclidean distance\nin %g-d space \n(N = %g)', denoised_dim, sum(select_for_PCA_B)));
        
        SetFigure(15);
        
        % =========== Animation ===========
        
        h_timeline = plot([PCA_B_times(1) PCA_B_times(1)],ylim,'m','linew',3);
        
        set(h_timeaxis,'xlim',[PCA_B_times(1)  PCA_B_times(end)],'ytick',[],'ycolor','w');
        
        uicontrol('Style','pushbutton','Unit','norm','Position',[0.032 0.022 0.085 0.054],...
            'String','Play','FontSize',13,...
            'Callback',{@Population_dynamics,0}); % Play
        uicontrol('Style','pushbutton','Unit','norm','Position',[0.132 0.022 0.085 0.054],...
            'String','Record','FontSize',13,...
            'Callback',{@Population_dynamics,10}); % Record
        %     uicontrol('Style','pushbutton','Unit','norm','Position',[0.229 0.022 0.052 0.054],...
        %         'String','<<','FontSize',13,...
        %         'Callback',{@Population_dynamics,-5}); % Backward
        %     uicontrol('Style','pushbutton','Unit','norm','Position',[0.292 0.022 0.052 0.054],...
        %         'String','<','FontSize',13,...
        %         'Callback',{@Population_dynamics,-1}); % Forward
        %     uicontrol('Style','pushbutton','Unit','norm','Position',[0.352 0.022 0.052 0.054],...
        %         'String','>','FontSize',13,...
        %         'Callback',{@Population_dynamics,1}); % Backward
        %     uicontrol('Style','pushbutton','Unit','norm','Position',[0.412 0.022 0.052 0.054],...
        %         'String','>>','FontSize',13,...
        %         'Callback',{@Population_dynamics,5}); % Forward
        uicontrol('Style','slider','Unit','norm','Position',[0.026 0.53 0.344 0.03],...
            'Min',1,'Max',length(PCA_B_times),'Value',1,...
            'Callback',{@Population_dynamics,999});
        
        % ======= Momentum ========
        figure(532);  clf;
        
        plot(repmat(PCA_B_times(1:end-1)',1,stim_type_num), diff(distance),'linew',3); hold on;
        axis([-50 1800 -100 250]);
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_PCA_B}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',2,'color',[0.6 0.6 0.6]);
        for tt = 1:3
            plot([1 1] * time_markers{j_PCA_B}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_PCA_B}{tt},'linew',2);
        end
        SetFigure(15);
        
        % Note here I use a double nested function. @HH20150425
        % Play population dynamics. HH20150424
        function Population_dynamics(~,~,flag)
            switch flag
                case {0,10}   % Play & record
                    % Preallocate movie structure.
                    mov(1:size(PCA_B_projPC{1},2)) = struct('cdata', [],...
                        'colormap', []);
                    
                    for t_ind = 1 : size(PCA_B_projPC{1},2)
                        for kk = 1:3
                            % Update
                            set(h_pref(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}((kk-1)*2+1,t_ind),'ydata',PCA_B_projPC{which_two_dimension(2)}((kk-1)*2+1,t_ind));
                            set(h_null(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}(kk*2,t_ind),'ydata',PCA_B_projPC{which_two_dimension(2)}(kk*2,t_ind));
                        end
                        
                        %     set(h_text,'string',num2str(fix(rate_ts{j_PCA_B}(tt)/10)*10));
                        set(h_timeline,'xdata',PCA_B_times(t_ind)*[1 1]);
                        
                        if flag
                            mov(t_ind) = getframe(gcf);
                        else
                            pause(0.01);
                        end
                        
                        drawnow;
                    end
                    
                    if flag ; movie2avi(mov,'LIP_dynamics_j_1.avi','compression','none'); end
                    
                case 999  % From slider bar
                    t_ind = round(get(gcbo,'value'));
                    
                    for kk = 1:3
                        % Update
                        set(h_pref(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}((kk-1)*2+1,t_ind),'ydata',PCA_B_projPC{which_two_dimension(2)}((kk-1)*2+1,t_ind));
                        set(h_null(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}(kk*2,t_ind),'ydata',PCA_B_projPC{which_two_dimension(2)}(kk*2,t_ind));
                    end
                    
                    %     set(h_text,'string',num2str(fix(rate_ts{j_PCA_B}(tt)/10)*10));
                    set(h_timeline,'xdata',PCA_B_times(t_ind)*[1 1]);
                    
                    
                    %             otherwise % Step back or forward
                    %                 % Get current time
                    %                 curr_t = get(h_timeline,'xdata');
                    %                 curr_t = curr_t(1);
                    %                 [~,curr_t_ind] = min(abs(PCA_B_times-curr_t));
                    %
                    %                 % Update to new time
                    %                 curr_t_ind = curr_t_ind + flag
                    %                 if 0 < curr_t_ind && curr_t_ind <= sum(PCA_B_time_range)
                    %                     for kk = 1:3
                    %                         % Update
                    %                         set(h_pref(kk),'xdata',projPC{which_two_dimension(1)}((kk-1)*2+1,curr_t_ind),'ydata',projPC{which_two_dimension(2)}((kk-1)*2+1,curr_t_ind));
                    %                         set(h_null(kk),'xdata',projPC{which_two_dimension(1)}(kk*2,curr_t_ind),'ydata',projPC{which_two_dimension(2)}(kk*2,curr_t_ind));
                    %                     end
                    %
                    %                     %     set(h_text,'string',num2str(fix(rate_ts{j_PCA_B}(tt)/10)*10));
                    %                     set(h_timeline,'xdata',PCA_B_times(curr_t_ind)*[1 1]);
                    %                 end
            end
        end
        
    end


%%%%%%%%%%%%%%%%%%%  6. SVM Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ SVM Training ------
select_for_SVM = select_bottom_line;
j_for_SVM = 1;
min_reps = 15; % For screening cells. (Anne used 20)
training_reps = 25; % I artificially added some trials by bootstrapping if reps < training_reps.
% This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
SVM_training_epoch = [700 1100]; % Anne used 500~700 ms
n_bootstrap_weights = 1000; % Anne used 1000
n_bootstrap_threshold = 50; % Anne used 25

% ------ SVM Testing ------
SVM_testing_span = 100; % Anne used 100 ms
SVM_testing_step = 50;
n_bootstrap_test = 300; % Anne used 1000

% ------ Initialization -----
pseudo_trial_pool_training = [];  pseudo_trial_pool_testing = [];
SVM_testing_tcenters = min(rate_ts{j_for_SVM}) + SVM_testing_span/2 : SVM_testing_step : max(rate_ts{j_for_SVM}) - SVM_testing_span/2;
answers_choice = []; answers_modality = []; reps_actual = [];
weights_svm_choice_allbootstrap = []; weights_svm_modality_allbootstrap = []; angles = [];
weights_svm_choice_mean = []; weights_svm_modality_mean = [];
thres_choice = []; thres_modality = []; select_for_SVM_actual = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function f6p0(debug)      % SVM 1: Training @HH20150429
        if debug;  dbstack;  keyboard;  end
        
        
        % Override
        %         min_reps = 15; % Anne used 20
        %         training_reps = 25; % I artificially added some trials by bootstrapping if reps < training_reps.
        %                             % This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
        
        SVM_training_epoch = [0 1700];    % This make it comparable with PCA results
        
        find_for_SVM = find(select_for_SVM);
        SVM_training_epoch_ind = min(SVM_training_epoch) <= rate_ts{j_for_SVM} & rate_ts{j_for_SVM} <= max(SVM_training_epoch);
        
        % Regroup pseudo-trials pool
        pseudo_trial_pool_training = cell(length(find_for_SVM),4);
        pseudo_trial_pool_testing = cell(length(find_for_SVM),4,length(SVM_testing_tcenters));
        
        for ii = 1:length(find_for_SVM)
            %             raw_this = group_result(find_for_SVM(ii)).mat_raw_PSTH.PSTH{j_for_SVM,ALL_CHOICE,1}.raw([1 2 3 4]); % Vest & vis
            raw_this = group_result(find_for_SVM(ii)).mat_raw_PSTH.PSTH{j_for_SVM,ALL_CHOICE,1}.raw([1 2 5 6]); % Vest & Comb
            %             raw_this = group_result(find_for_SVM(ii)).mat_raw_PSTH.PSTH{j_for_SVM,ALL_CHOICE,1}.raw([3 4 5 6]); % Vis & Comb
            
            % By default, SVM for choice is related to pref and null
            % This make it comparable with PCA results and reasonable for TDR
            
            % Flip to contralateral vs ipsilateral
            %             if ~group_result(find_for_SVM(ii)).if_contralateral
            %                 raw_this = raw_this([2 1 4 3]);  % +/-,+/-
            %             end
            
            % Flip to right vs left
            %             if group_result(find_for_SVM(ii)).PREF_PSTH == 1
            %                 raw_this = raw_this([2 1 4 3]);
            %             end
            
            % Training pool (one certain time window)
            means_this = cellfun(@(x)mean(x(:,SVM_training_epoch_ind),2),raw_this,'uniformOutput',false);
            [pseudo_trial_pool_training{ii,:}] = means_this{:};
            
            % Testing pool (shifting time windows)
            for ttt = 1:length(SVM_testing_tcenters)
                test_epoch_ind = SVM_testing_tcenters(ttt) - SVM_testing_span/2 <= rate_ts{j_for_SVM} & rate_ts{j_for_SVM} <= SVM_testing_tcenters(ttt) + SVM_testing_span/2;
                means_this = cellfun(@(x)mean(x(:,test_epoch_ind),2),raw_this,'uniformOutput',false);
                
                [pseudo_trial_pool_testing{ii,:,ttt}] = means_this{:};
            end
            
        end
        
        % Find cells who cross the minimal repetitions for each condition (Anne)
        lengths = cellfun(@(x)length(x),pseudo_trial_pool_training);
        select_for_SVM_actual = all(lengths >= min_reps,2);
        pseudo_trial_pool_training = pseudo_trial_pool_training(select_for_SVM_actual,:);
        pseudo_trial_pool_testing = pseudo_trial_pool_testing(select_for_SVM_actual,:,:);
        
        lengths = lengths(select_for_SVM_actual,:);
        
        % I artificially added some trials by bootstrapping if reps < training_reps.
        % This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
        
        who_needs_add_some_trials = find(lengths < training_reps);
        
        for aa = 1:length(who_needs_add_some_trials)
            this_who = who_needs_add_some_trials(aa);
            
            % Select trials
            real_reps = length(pseudo_trial_pool_training{this_who});
            n_to_add = training_reps - real_reps;
            trials_to_add = randperm(real_reps,n_to_add);
            
            % Add trials
            pseudo_trial_pool_training{this_who} = [pseudo_trial_pool_training{this_who}; pseudo_trial_pool_training{this_who}(trials_to_add)];
            [sub1,sub2] = ind2sub(size(pseudo_trial_pool_training),this_who);
            
            for ttt = 1:length(SVM_testing_tcenters)
                pseudo_trial_pool_testing{sub1,sub2,ttt} = [pseudo_trial_pool_testing{sub1,sub2,ttt}; pseudo_trial_pool_testing{sub1,sub2,ttt}(trials_to_add)];
            end
        end
        
        % Recalculate lengths
        lengths = cellfun(@(x)length(x),pseudo_trial_pool_training);
        reps_actual = min(lengths);
        
        % -------- Teacher signals --------
        % For choices, 1 = Contralateral, 0 = Ipsilateral
        answers_choice = [ones(1,reps_actual(1)) zeros(1,reps_actual(2)) ones(1,reps_actual(3)) zeros(1,reps_actual(4))]';
        % For modality, 0 = Vestibular, 1 = Visual
        answers_modality = [zeros(1,reps_actual(1)) zeros(1,reps_actual(2)) ones(1,reps_actual(3)) ones(1,reps_actual(4))]';
        
        %% =============  SVM Bootstrapping ============
        % Parallel computing
        if matlabpool('size') == 0 ;  matlabpool;  end
        training_set = cell(n_bootstrap_weights);
        tic; disp('Starting SVM training...');
        
        parfor_progress(n_bootstrap_weights);
        
        parfor nn = 1:n_bootstrap_weights
            % Randomly select trials without replacement from the pseudo-trials pool
            pseudo_trial_pool_perm = cellfun(@(x)x(randperm(size(x,1)),1),pseudo_trial_pool_training,'uniform',false);
            
            for conds = 1:4
                training_set{nn}{conds} = cell2mat(cellfun(@(x)x(1:reps_actual(conds),:),pseudo_trial_pool_perm(:,conds),'uniform',false)');
            end
            
            firing_for_this_training = cell2mat(training_set{nn}');
            
            % Train two SVM classifiers (choice and modality)
            svm_training_choice(nn) = svmtrain(firing_for_this_training,answers_choice);
            svm_training_modality(nn) = svmtrain(firing_for_this_training,answers_modality);
            
            parfor_progress;
        end
        parfor_progress(0);
        toc; disp('Done!');
        
        %% Averaged weights (bagging)
        weights_svm_choice_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_choice,'uniform',0));
        weights_svm_modality_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_modality,'uniform',0));
        
        % Orthogonal test
        angles = acos(sum(weights_svm_choice_allbootstrap.*weights_svm_modality_allbootstrap,1)./...
            sqrt(sum(weights_svm_choice_allbootstrap.^2,1)./ sum(weights_svm_modality_allbootstrap.^2,1)))/pi*180;
        
        weights_svm_choice_mean = mean(weights_svm_choice_allbootstrap,2);
        weights_svm_choice_mean = weights_svm_choice_mean / max(abs(weights_svm_choice_mean));
        % svm_weights_choice_sem = std(svm_weights_choice,[],2)/sqrt(n_bootstrap);
        
        weights_svm_modality_mean = mean(weights_svm_modality_allbootstrap,2);
        weights_svm_modality_mean = weights_svm_modality_mean / max(abs(weights_svm_modality_mean));
        % svm_weights_modality_sem = std(svm_weights_modality,[],2)/sqrt(n_bootstrap);
        
        %% Find proper threshold for SVM
        training_set = cell(n_bootstrap_threshold);
        bestZ_choice = nan(n_bootstrap_threshold,1);
        bestZ_modality = nan(n_bootstrap_threshold,1);
        tic; disp('Find proper threshold for SVM decoders...');
        
        parfor nn = 1:n_bootstrap_threshold
            % Randomly select trials without replacement from the pseudo-trials pool
            pseudo_trial_pool_perm = cellfun(@(x)x(randperm(size(x,1)),1),pseudo_trial_pool_training,'uniform',false);
            
            for conds = 1:4
                training_set{nn}{conds} = cell2mat(cellfun(@(x)x(1:reps_actual(conds),:),pseudo_trial_pool_perm(:,conds),'uniform',false)');
            end
            
            proj_choice_on_choice = cell2mat(training_set{nn}') * weights_svm_choice_mean;
            proj_modality_on_modality = cell2mat(training_set{nn}') * weights_svm_modality_mean;
            
            % Find the best threshold (the most northwestern point on ROC curve)
            [~,bestZ_choice(nn)] = rocN(proj_choice_on_choice(logical(answers_choice)),proj_choice_on_choice(~logical(answers_choice)));
            [~,bestZ_modality(nn)] = rocN(proj_modality_on_modality(logical(answers_modality)),proj_modality_on_modality(~logical(answers_modality)));
            
        end
        toc; disp('Done!');
        
        thres_choice = mean(bestZ_choice);
        thres_modality = mean(bestZ_modality);
        
        % ================ SVM Training End ===============
        
    end
    function f6p1(debug)      % SVM 2: Plotting SVM weights
        if debug;  dbstack;  keyboard;  end
        
        
        if isempty(thres_choice)
            f6p0(0);
        end
        
        % ------ Plotting -----
        set(figure(611),'pos',[89 26 761 606]); clf
        
        find_for_SVM = find(select_for_SVM);
        find_for_SVM_actual_of_all = find_for_SVM(select_for_SVM_actual);
        who_are_tcells = select_tcells(find_for_SVM_actual_of_all);
        
        % Weights for choice decoder
        subplot(2,2,1);
        [~,sortt] = sort(abs(weights_svm_choice_mean),'descend');
        set(bar(weights_svm_choice_mean(sortt),0.5),'facecolor','k','edgecolor','none'); hold on;
        set(bar(find(who_are_tcells(sortt)),...
            weights_svm_choice_mean(sortt(who_are_tcells(sortt))),0.5),'facecolor','r','edgecolor','none');
        
        xlim([-2 sum(select_for_SVM_actual)+1]); ylim([-1 1]); title('Choice');
        set(gca,'xtick',[1 sum(select_for_SVM_actual)]);
        % errorbar(svm_weights_choice_mean(sortt),svm_weights_choice_sem(sortt),'.');
        
        % Weights for modality decoder
        subplot(2,2,2);
        [~,sortt] = sort(abs(weights_svm_modality_mean),'descend');
        set(bar(weights_svm_modality_mean(sortt),0.5),'facecolor','k','edgecolor','none'); hold on;
        set(bar(find(who_are_tcells(sortt)),...
            weights_svm_modality_mean(sortt(who_are_tcells(sortt))),0.5),'facecolor','r','edgecolor','none');
        
        xlim([-2 sum(select_for_SVM_actual)+1]);  ylim([-1 1]); title('Modality');
        set(gca,'xtick',[1 sum(select_for_SVM_actual)]);
        % errorbar(svm_weights_modality_mean(sortt),svm_weights_modality_sem(sortt),'.');
        
        % Correlation
        monkeys = xls_num{1}(:,header.Monkey);  % Temporarily here. HH20150914
        monkey1 = monkeys == 5;
        monkey2 = monkeys == 10;
        
        ax = subplot(2,2,3);
        h = LinearCorrelation({
            (weights_svm_modality_mean(monkey1(select_for_SVM)));
            (weights_svm_modality_mean(monkey2(select_for_SVM)))
            },...
            {
            (weights_svm_choice_mean(monkey1(select_for_SVM))) ;
            (weights_svm_choice_mean(monkey2(select_for_SVM)))
            },...
            'CombinedIndex',[3],...
            'Xlabel','Weights in modality decoder','Ylabel','Weights in choice decoder',...
            'FaceColors',{'k'},'Markers',{'o','^'},...
            'LineStyles',{'k-'},'MarkerSize',6,...
            'axes',ax,'XHist',0,'YHist',0,...
            'SameScale',1,'Method','Spearman');
        delete([h.diag h.group(1:2).line]); legend off;
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
        text(min(xlim)+0.2,min(ylim)+0.2,sprintf('r^2 = %g, p = %g',h.group(3).r_square,h.group(3).p),'fonts',11);
        
        % Annotate tcells
        h_line = plot(weights_svm_modality_mean(who_are_tcells),...
            weights_svm_choice_mean(who_are_tcells),...
            '+','markersize',8,'color','r','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        select_for_SVM_actual_of_all = zeros(N,1);
        select_for_SVM_actual_of_all(find_for_SVM_actual_of_all) = 1;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h.group(1).dots, select_for_SVM_actual_of_all});
        
        % Angles between weights
        subplot(2,2,4);
        [n,x] = hist(angles,20); bar(x,n,'facecolor','k','edgecolor','k'); hold on;
        plot([median(angles) median(angles)],ylim,'k-');
        plot([90 90],ylim,'k--'); axis tight; xlim([0 180]); set(gca,'xtick',0:90:180);
        title(sprintf('Median angle = %g',median(angles)));
        
        SetFigure(15);
        
        % ---------  Memsac DDI and weights for choice decoder
        set(figure(612),'pos',[89 26 761 606]); clf
        
        MemSac_DDI_phase = [2:4];
        
        h = LinearCorrelation({
            MemSac_indicator(select_for_SVM);
            },...
            {
            (weights_svm_choice_mean);
            },...
            'CombinedIndex',[],...
            'Xlabel',MemSac_indicator_txt,'Ylabel','Weight in choice decoder',...
            'FaceColors',{'k'},'Markers',{'o'},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',612,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman');
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),weights_svm_choice_mean(select_tcells(select_for_SVM),1),...
            '+','markersize',16,'color','r','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(MemSac_indicator(select_for_SVM),(weights_svm_choice_mean(:,1)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_SVM});
        
        
        % ---------  Weights for eigen-neuron 1 and weights for SVM choice decoder
        if isempty(weights_PCA_B_PC)
            f5p0(0);
        end
        
        set(figure(613),'pos',[89 26 761 606]); clf
                
        h = LinearCorrelation({
            weights_PCA_B_PC(:,1);
            },...
            {
            (weights_svm_choice_mean);
            },...
            'CombinedIndex',[],...
            'Xlabel','Weight for eigen-neuron 1','Ylabel','Weight in choice decoder',...
            'FaceColors',{'k'},'Markers',{'o'},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',613,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman');
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(weights_PCA_B_PC(select_tcells(select_for_SVM),1),weights_svm_choice_mean(select_tcells(select_for_SVM),1),...
            '+','markersize',16,'color','r','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        set([gca h.group.dots],'ButtonDownFcn',{@Show_individual_cell, h.group.dots, select_for_SVM});
        
    end
    function f6p2(debug)      % SVM 3: Testing SVM
        if debug;  dbstack;  keyboard;  end
        
        
        if isempty(thres_choice)
            f6p0(0);
        end
        
        corr_rate_choice_by_choice = nan(n_bootstrap_test,length(SVM_testing_tcenters));
        corr_rate_modality_by_modality = corr_rate_choice_by_choice;
        corr_rate_choice_by_modality = corr_rate_choice_by_choice;
        corr_rate_modality_by_choice = corr_rate_choice_by_choice;
        
        tic; disp('Testing SVM decoders...');
        
        progressbar('Testing SVM');
        
        for nn = 1:n_bootstrap_test  % For each bootstrapping
            % Generate rand perm indices
            pseudo_trial_pool_perm_ind = cellfun(@(x)randperm(size(x,1)),pseudo_trial_pool_training,'uniform',false);
            
            for ttt = 1:length(SVM_testing_tcenters)
                
                for conds = 1:4
                    for ii = 1:size(pseudo_trial_pool_perm_ind,1)
                        testing_set{nn}{conds}(:,ii) = pseudo_trial_pool_testing{ii,conds,ttt}(pseudo_trial_pool_perm_ind{ii,conds}(1:reps_actual(conds)));
                    end
                end
                
                % Projection on SVM weights
                proj_choice_on_choice = cell2mat(testing_set{nn}') * weights_svm_choice_mean;
                proj_modality_on_modality = cell2mat(testing_set{nn}') * weights_svm_modality_mean;
                
                % Classification by thresholds
                class_test_choice = proj_choice_on_choice >= thres_choice;
                class_test_modality = proj_modality_on_modality >= thres_modality;
                
                % Correct rates
                corr_rate_choice_by_choice(nn,ttt) = sum(answers_choice == class_test_choice)/length(answers_choice);
                corr_rate_modality_by_modality(nn,ttt) = sum(answers_modality == class_test_modality)/length(answers_modality);
                corr_rate_choice_by_modality(nn,ttt) = sum(answers_choice == class_test_modality)/length(answers_choice);
                corr_rate_modality_by_choice(nn,ttt) = sum(answers_modality == class_test_choice)/length(answers_modality);
                
            end
            
            progressbar(nn/n_bootstrap_test);
        end
        
        toc; disp('Done!');
        
        %% --------- Plotting ---------
        
        set(figure(621),'pos',[114 223 930 353]); clf
        subplot(1,2,1);
        errorbar(SVM_testing_tcenters,nanmean(corr_rate_choice_by_choice),...
            nanstd(corr_rate_choice_by_choice),'color','m','linew',2); hold on;
        errorbar(SVM_testing_tcenters,nanmean(corr_rate_modality_by_choice),...
            nanstd(corr_rate_modality_by_choice),'color',[.87 .49 0],'linew',2);
        legend('Decodes choice','Decodes modality','location','best');
        
        axis tight; ylim([0 1]); plot(xlim,[0.5 0.5],'k--'); title('Choice decoder');
        
        % Time markers
        for tt = 1:3
            plot([1 1] * time_markers{j_for_SVM}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_for_SVM}{tt},'linew',1);
        end
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_for_SVM}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
        subplot(1,2,2);
        errorbar(SVM_testing_tcenters,nanmean(corr_rate_modality_by_modality),...
            nanstd(corr_rate_modality_by_modality),'color',[.87 .49 0],'linew',2); hold on;
        errorbar(SVM_testing_tcenters,nanmean(corr_rate_choice_by_modality),...
            nanstd(corr_rate_choice_by_modality),'color','m','linew',2);
        
        axis tight; ylim([0 1]); plot(xlim,[0.5 0.5],'k--'); title('Modality decoder');
        
        % Time markers
        for tt = 1:3
            plot([1 1] * time_markers{j_for_SVM}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_for_SVM}{tt},'linew',1);
        end
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_for_SVM}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
        SetFigure(15);
        
        %%
        %         keyboard
    end
    function f6p3(debug)      % SVM 4: Mean PSTH projectd on weight vector (SVM-weighted sum)
        if debug;  dbstack;  keyboard;  end
        
        
        if isempty(thres_choice)
            f6p0(0);
        end
        
        find_for_SVM = find(select_for_SVM);
        find_for_SVM_actual_of_all = find_for_SVM(select_for_SVM_actual);
        select_for_SVM_actual_of_all = false(N,1);
        select_for_SVM_actual_of_all(find_for_SVM_actual_of_all) = true;
        
        projs = {weights_svm_choice_mean,weights_svm_modality_mean};
        selects = {select_for_SVM_actual_of_all,select_for_SVM_actual_of_all};
        
        h = Weighted_sum_PSTH(projs,{'Weighted by svm\_choice\_mean','Weighted by svm\_modality\_mean'},selects);
        set(h.fig,'name','Projected on SVM decoders (Correct only, all choices)');
        
    end

    function f7p1(debug)      % Regression 1. Comb = w1 Vest + w2 Vis. @HH20150427
        if debug;  dbstack;  keyboard;  end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j_for_reg = 1;
        ts = rate_ts{j_for_reg};
%         select_for_regression = select_bottom_line;
        select_for_regression = select_tcells;

        t_span = 1000;    t_step = 50; % For time evolution
        toi = [950 1500] +  time_markers{j_for_reg}(1,1);   % Time of interests for weights illustration
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        find_for_regression = find(select_for_regression);
        t_centers = ts(1) + t_span / 2 : t_step : ts(end)-t_span / 2;
        
        weight_vest_vis = nan(length(find_for_regression),2 * 2 + 2,length(t_centers)); % Weights (3), p_values (3), r^2, p-value of r^2
        measure_pred_all = nan(length(find_for_regression) * 10,2,length(t_centers));
        
        progressbar('Number of cells');
        
        for i = 1:length(find_for_regression)  % For each cell
            
            for tcc = 1:length(t_centers) % For each time bin
                
                t_range = t_centers(tcc) - t_span / 2 <= ts & ts <= t_centers(tcc) + t_span / 2;
                
                r = [];
                for k = 1:3
                    r(:,k) = mean(group_result(find_for_regression(i)).mat_raw_PSTH.PSTH{j_for_reg,2,k}.ys(:,t_range),2);
                end
                
                r = r - repmat(nanmean(r,1),size(r,1),1); % Mean removed
                
                % GLM fit
                [b,~,stat] = glmfit(r(:,1:2) ,r(:,3),'normal','link','identity','constant','off');
                [~,~,~,~,stat_reg] = regress(r(:,3),[ones(size(r,1),1) r(:,1:2)]);
                
                weight_vest_vis(i,:,tcc) = [b' stat.p' stat_reg([1 3])]; % Weights, r^2 and p-value
                
                yfit = r(:,1:2) * b;
                measure_pred_all((i-1)*10 + 1:(i-1)*10+size(r,1),:,tcc) = [r(:,3) yfit];
                
                % r = r([end:-2:2 1:2:end-1],:);
                % yfit = glmval(b,r(:,1:2),'identity');
                % figure(348); hold on; plot(1:10,r(:,3),'o',1:10,yfit,'r-');
                
            end
            progressbar(i/length(find_for_regression));
        end
        
        %% ========= Plotting ========
        
        % ----Time-dependent regression parameters
        set(figure(711),'pos',[20 390 674 572]);  clf; hold on;
        mean_paras = squeeze(mean(weight_vest_vis(:,[1 2 end-1],:),1))';
        sem_paras = squeeze(std(weight_vest_vis(:,[1 2 end-1],:)))'/sqrt(size(weight_vest_vis,1));
        
        temp_col = {'b','r','k'};
        temp_marker = {'s','o','o'}; % {'','',''};
        for pp = 1 : 3
            h = shadedErrorBar(repmat(t_centers',1,1),mean_paras(:,pp),sem_paras(:,pp),{[temp_marker{pp} temp_col{pp} '-'],'linew',2,'markersize',10});
            hold on;
        end
        
        set(legend('w_{vest}','w_{vis}','R^2'),'location','best'); axis tight; %ylim([0 max(ylim)])
        xlabel(sprintf('Center of %g ms window',t_span)); ylabel('Median');  SetFigure();
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_for_SVM}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);

        
        for tt = 1:3
            plot([1 1] * time_markers{j_for_reg}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_for_reg}{tt},'linew',2);
        end
        
        
        % ---- R^2 and weight. Figure 5 of Gu 2008
        
        set(figure(612),'pos',[14 96 1319 864]);  clf;
        
        for toii = 1:2
            [~,toi_ind] = min(abs(toi(toii)-t_centers));
            
            sig_ind = weight_vest_vis(:,end,toi_ind)<0.05;
            nsig_ind = weight_vest_vis(:,end,toi_ind)>=0.05;
            
            % Predicted v.s. measured
            ax0 = subplot(2,3,1 + (toii-1)*3);
            h = LinearCorrelation( measure_pred_all(:,1,toi_ind), measure_pred_all(:,2,toi_ind),...
                'CombinedIndex',[],...
                'Xlabel','Measured response','Ylabel','Predicted response',...
                'FaceColors',{[0.9 0.9 0.9]},'Markers',{'.'},...
                'LineStyles',{'k-'},'MarkerSize',10,...
                'XHist',0,'YHist',0,...
                'XHistStyle','grouped','YHistStyle','grouped','SameScale',1,'Method','Spearman','axes',ax0);
            text(min(xlim)+2,min(ylim)+2,sprintf('r^2 = %g, p = %g',h.group(1).r_square,h.group(1).p),'fonts',13);
            legend off;
            
            % Distribution of R^2
            ax1 = subplot(2,3,2 + (toii-1)*3);
            h = HistComparison({weight_vest_vis(sig_ind,end-1,toi_ind),weight_vest_vis(nsig_ind,end-1,toi_ind)},...
                'EdgeColors',{'k','k'},'FaceColors',{[.3 .3 .3 ],[0.9 0.9 0.9]},'XCenters',0.05:0.1:1,'MeanType','Median','Axes',ax1);
            xlabel('R^2'); ylabel('Number of neurons'); xlim([0 1]); ylim([0 23]);
            title(sprintf('t_{center} = %g',t_centers(toi_ind)));
            
            % Vestibular / visual weights
            ax2 = subplot(2,3,3 + (toii-1)*3);
            plot(weight_vest_vis(:,1,toi_ind),weight_vest_vis(:,2,toi_ind),'o');
            
            h = LinearCorrelation({
                weight_vest_vis(nsig_ind,1,toi_ind);
                weight_vest_vis(sig_ind,1,toi_ind);
                },...
                {
                weight_vest_vis(nsig_ind,2,toi_ind) ;
                weight_vest_vis(sig_ind,2,toi_ind) ;
                },...
                'CombinedIndex',[],...
                'Xlabel','Vestibular weight','Ylabel','Visual weight',...
                'FaceColors',{[0.9 0.9 0.9],[.3 .3 .3 ]},'Markers',{'o'},...
                'LineStyles',{'k--','k-'},'MarkerSize',9,...
                'XHist',20,'YHist',10,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,'Method','Spearman','axes',ax2);
            
            if toii == 1 ; xlabel(''); end
            delete([h.group(1).line h.diag]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--'); plot(xlim,[1 1],'k--'); plot([1 1],ylim,'k--');
            text(min(xlim)+0.2,min(ylim)+0.2,sprintf('r^2 = %g, p = %g',h.group(2).r_square,h.group(2).p),'fonts',13);
            legend off;
            set(gca,'xtick',-10:1:10,'ytick',-10:1:10);
            
            % Annotate tcells
            h_t = plot(weight_vest_vis(select_tcells(select_for_regression),1,toi_ind),...
                weight_vest_vis(select_tcells(select_for_regression),2,toi_ind),'+','markersize',10,'color','r','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(weight_vest_vis(:,1,toi_ind), weight_vest_vis(:,2,toi_ind),'visible','off'); hold on;
            set([gca [h.group.dots] h_line h_t],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_regression});
            
        end
        
        SetFigure(18);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%  8. TDR   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights_TDR_PCA_SVM_mean = [];
weights_TDR_PCA_SVM_allbootstrap = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function f8p1(debug)      % TDR: PCA + SVM
        if debug;  dbstack;  keyboard;  end
        
        % --------  Targeted dimensionality reduction ----------
        if isempty(weights_PCA_B_PC)
            f5p0(0);
        end
        
        if isempty(weights_svm_choice_mean)
            f6p0(0);
        end
        
        % ------ Mean weight --------
        % Denoised matrix from PCA
        D = weights_PCA_B_PC(:,1:denoised_dim) * weights_PCA_B_PC(:,1:denoised_dim)';
        
        % Project SVM weights in subspace spanned by the first denoised_dim principal components
        beta = [weights_svm_choice_mean weights_svm_modality_mean];
        beta_pca = D * beta;
        
        % Orthogonalization by QR-decomposition
        [Q,~] = qr(beta_pca);
        weights_TDR_PCA_SVM_mean = Q(:,1:2);   weights_TDR_PCA_SVM_mean(:,1) = - weights_TDR_PCA_SVM_mean(:,1);
        
        % Contribution of naive PCs
        pc_contribution =  weights_PCA_B_PC(:,1:denoised_dim)' * weights_TDR_PCA_SVM_mean;
        figure(64); bar(pc_contribution,1) ; SetFigure(); set(gca,'xtick',1:denoised_dim);

        % ------ Bootstrap weights --------
        weights_TDR_PCA_SVM_allbootstrap = nan(size(weights_svm_choice_allbootstrap,1),2,size(weights_svm_choice_allbootstrap,2));
        
        for bb = 1:size(weights_svm_choice_allbootstrap,2)
            
            % Project SVM weights in subspace spanned by the first denoised_dim principal components
            beta = [weights_svm_choice_allbootstrap(:,bb) weights_svm_modality_allbootstrap(:,bb)];
            beta_pca = D * beta;
            
            % Orthogonalization by QR-decomposition
            [Q,~] = qr(beta_pca);
            weights_TDR_PCA_SVM_allbootstrap(:,:,bb) = Q(:,1:2);   
            
            % Normalized by sum, not 2-norm
            weights_TDR_PCA_SVM_allbootstrap(:,1,bb) = (weights_TDR_PCA_SVM_allbootstrap(:,1,bb))/sum((weights_TDR_PCA_SVM_allbootstrap(:,1,bb)),1); 
            weights_TDR_PCA_SVM_allbootstrap(:,2,bb) = (weights_TDR_PCA_SVM_allbootstrap(:,2,bb))/sum((weights_TDR_PCA_SVM_allbootstrap(:,2,bb)),1); 
            
        end
        
        % -------- Weights vs cell properties ----------
        Weights_Property_Correlation(weights_TDR_PCA_SVM_mean,...
            {'Weights for TDR axis 1 (choice)','Weights for TDR axis 2 (modality)'}, select_for_PCA_B);
        
    end
    function f8p2(debug)      % TDR: PCA + SVM, correct trials / angles
        if debug;  dbstack;  keyboard;  end
        
        if isempty(weights_TDR_PCA_SVM_mean)
            f8p1(0);
        end
        
        % ------ Projection of raw data on task-related axes -------
        %         h = Weighted_sum_PSTH({weights_TDR_PCA_SVM_mean(:,1),weights_TDR_PCA_SVM_mean(:,2)},{'Weighted by TDR\_PCA\_SVM\_choice','Weighted by TDR\_PCA\_SVM\_modality'},...
        %             {select_for_SVM select_for_SVM});
        
        h = Weighted_sum_PSTH({squeeze(weights_TDR_PCA_SVM_allbootstrap(:,1,:)),squeeze(weights_TDR_PCA_SVM_allbootstrap(:,2,:))},...
            {'Weighted by TDR\_PCA\_SVM\_choice','Weighted by TDR\_PCA\_SVM\_modality'},...
            {select_for_SVM select_for_SVM});
        
        %
        figure(904);  clf; hold on;
        for k = 1:3
            plot(rate_ts{1}(1:end-1) ,diff(h.projected_all_mean{1}(k*2-1,:)-h.projected_all_mean{1}(k*2,:)),'color',colors(k,:),'linew',2);
        end
        SetFigure;
        
        %% ------ Tuning curves ------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
        tuning_centers = [ mean(time_markers{j}(1,1:2)) + 150  % Stim center + sensory delay
            mean(time_markers{j}(1,3)) - group_result(end).mat_raw_PSTH.binSize_CP/2   % Pre-sac epoch
            mean(time_markers{j}(1,3)) + group_result(end).mat_raw_PSTH.binSize_CP/2];  % Post-sac epoch
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,center_t_ind] = min(abs(rate_ts{j} - tuning_centers(1)));
        [~,pre_t_ind] = min(abs(rate_ts{j} - tuning_centers(2)));
        [~,post_t_ind] = min(abs(rate_ts{j} - tuning_centers(3)));
        
        tuning_time_phase = [center_t_ind pre_t_ind post_t_ind];
        tuning_time_phase_title = {'Stimulus', 'Pre-sac', 'Post-sac'};
        
        unique_heading = group_result(end).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
        unique_heading = [unique_heading(unique_heading<0);-0.00001; 0.00001; unique_heading(unique_heading>0)]; % Add another 0 heading
        
        % Get means and sems
%         tuning_mean_correctonly = permute(h.projected_angles_mean{j},[2 1 3]);
        tuning_mean_correctonly = h.projected_angles_mean{j};
        tuning_sem_correctonly = h.projected_angles_sem{j};
        
        smooth_time_win = 500;
        for hh = 1:size(tuning_mean_correctonly,1)
            for kk = 1:size(tuning_mean_correctonly,3)
                tuning_mean_correctonly(hh,:,kk) = smooth(tuning_mean_correctonly(hh,:,kk), round(smooth_time_win/mean(diff(rate_ts{1}))));
            end
        end
        
        % Sort to match 'unique_heading'
        tuning_mean_correctonly =  [tuning_mean_correctonly(end:-2:1,:,:); tuning_mean_correctonly(1:2:end,:,:)];
        
        %         for pp = 1:length(rate_ts{j})
        %
        %             for k = 1:3
        %
        %                 % Mean and sem
        %                 this_tuning_correctonly = tuning_pack{2}{k,pp}(:,~isnan(tuning_pack{2}{k,pp}(1,:)));
        %                 tuning_mean_correctonly(:,pp,k) = mean(this_tuning_correctonly,2);
        %
        %                 % Fitting sigmoid function
        %
        %             end
        %
        %         end
        
        % Plotting tuning curves at three time points
        
        set(figure(2099+figN),'pos',[9 546 1378 414],'name',['Tuning curve, j = ' num2str(j)]); clf; figN = figN+1;
        
        for pp = 1:length(tuning_time_phase)
            
            for k = 1:3  % For each stim type
                
                % Plotting
                subplot(1,length(tuning_time_phase),pp );  hold on; ylabel('Correct only');
                plot(unique_heading(end/2+1:end),tuning_mean_correctonly(end/2+1:end,tuning_time_phase(pp),k),'-o','markersize',9,'color',colors(k,:),'markerfacecol',colors(k,:),'LineWid',2);
                plot(unique_heading(1:end/2),tuning_mean_correctonly(1:end/2,tuning_time_phase(pp),k),'-o','markersize',9,'color',colors(k,:),'markerfacecol','none','LineWid',2);
                
                h_tmp = errorbar(unique_heading,tuning_mean_correctonly(:,tuning_time_phase(pp),k),tuning_sem_correctonly(:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2);
                errorbar_tick(h_tmp,10000);
                
                title(num2str(tuning_centers(pp)));
                axis tight;
                
                SetFigure(20);
            end
        end
        
        %% Tuning: Animation
        
%         set(figure(2099+figN),'pos',[9 309 788 650],'name',['Tuning curve, j = ' num2str(j)]); clf; figN = figN+1;
%         ylabel('Correct only');
%         SetFigure(15); drawnow;
%         
%         % Preallocate movie structure.
%         mov_tuning(1:length(rate_ts{j})) = struct('cdata', [], 'colormap', []);
%         
%         for pp = 1:length(rate_ts{j})
%             
%             hold off;
%             
%             for k = 1:3  % For each stim type
%                 
%                 % Plotting
%                 plot(unique_heading(end/2+1:end),tuning_mean_correctonly(end/2+1:end,pp,k),'-o','markersize',9,'color',colors(k,:),'markerfacecol',colors(k,:),'LineWid',2);
%                 hold on;
%                 plot(unique_heading(1:end/2),tuning_mean_correctonly(1:end/2,pp,k),'-o','markersize',9,'color',colors(k,:),'markerfacecol','none','LineWid',2);
%                 
%                 h_tmp = errorbar(unique_heading,tuning_mean_correctonly(:,pp,k),tuning_sem_correctonly(:,pp,k),'.','color',colors(k,:),'LineWid',2);
%                 errorbar_tick(h_tmp,10000);
%             end
%             
%             title(num2str(rate_ts{j}(pp)));
%             SetFigure(15);
%             
%             %             axis tight;
%             
%             mov_tuning(pp) = getframe(gcf);
%             
%             drawnow;
%             pause(0.01);
%         end
%         
%         movie2avi(mov_tuning,'LIP_tuning_evolve_j_1_TDR.avi','compression','none');
        
        %% Tuning: Hotgram
        
        to_evolve =  {tuning_mean_correctonly};
        to_evolve_title = {'Correct only'};
        
        for co = 1
            set(figure(2099+figN),'pos',[82 99 927 855],'name','Evolve of tuning curves'); clf; figN = figN+1;
            
            [X,Y] = meshgrid(rate_ts{j},unique_heading);
            [Xq,Yq] = meshgrid(linspace(min(rate_ts{j}),max(rate_ts{j}),1000),linspace(min(unique_heading),max(unique_heading),100));
            
            
            for k = 1:3  % For each stim type
                subplot_tight(1,3,k,0.02,0.1);
%                 surf(Xq,Yq,interp2(X,Y,to_evolve{co}(:,:,k),Xq,Yq),'edgecolor','none'); view(90,-90); hold on;
                contourf(Xq,Yq,interp2(X,Y,to_evolve{co}(:,:,k),Xq,Yq),40,'edgecolor','k'); view(90,-90); hold on;
                
                if k>1 ;set(gca,'xtick',[]); end
                if k==2 ; title(to_evolve_title{co}); end
                
                caxis([min(to_evolve{co}(:)) 1.05*max(to_evolve{co}(:))]); % Use the same color range
                xlim([min(CP_ts{j}) max(CP_ts{j})]); ylim([min(unique_heading) max(unique_heading)]);
                
                for tt = 1:3
                    plot3([1 1] * time_markers{j}(1,tt),ylim,-3*[1 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
                end
                
                load('colormap_for_tuning_hotgram_TDR','colormap_for_tuning_hotgram_TDR');
                colormap(colormap_for_tuning_hotgram_TDR);
                
                %                 axis tight;
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
                
            end
            
            SetFigure();
        end
        
    end


    function f9p1(debug)      % Cell counter
        if debug;  dbstack;  keyboard;  end
        
        % Plot cell_num v.s. experimental days to see if I've been too lazy (which is beyond doubt)...
        dates = xls_num{1}(:,header.Date);
        day_counter = datenum(num2str(dates),'yyyymmdd');
        
        set(figure(9999),'name','My Progress','position',[17 258 649 693]); clf;
        plot([day_counter; today],[cumsum(select_all); sum(select_all)],'k-','linew',3); hold on;
        plot([day_counter; today],[cumsum(select_sus); sum(select_sus)],'r-','linew',3);
        plot([day_counter; today],[cumsum(select_tcells); sum(select_tcells)],'g-','linew',3);
        legend({'SUs + MUs','SUs (bottom line)','T SUs'},'location','best');
        
        plot(xlim,[100 100],'r--');
        plot([today today],ylim,'k-');
        
        xlims = xlim;
        set(gca,'XTick',linspace(xlims(1),xlims(2),12));
        datetick('x',26,'KeepTicks');
        
        rotateXLabels(gca,30);
        axis tight;
        SetFigure();
        
        % After running this part of code, I am literally crying...
        % Come on man!!!!!!!
        
        % Plot prediction ratio
        set(figure(9998),'name','Prediction ratio','position',[686 539 900 412]); clf;
        h_line = plot(Psy_pred_ratio,'ko','markerfacecol','k','markersize',10); hold on;
        plot(find(select_tcells),Psy_pred_ratio(select_tcells),'+r','linew',2,'markersize',10);
        plot(xlim,[1 1],'k--'); ylabel('Prediction ratio'); set(gca,'yscale','log','ytick',[0.6 1 2]);
        SetFigure(15);
        
        % Plot threshold for single cues
%          Psy_pred_ratio_vestibular_visual = sort(Psy_pred_ratio_vestibular_visual,2);
        plot(Psy_pred_ratio_vestibular_visual(:,1) ,'b-','linew',2);
        plot(Psy_pred_ratio_vestibular_visual(:,2) ,'r-','linew',2);

        
        % Show individual cell selected from the figure. HH20150424
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell,h_line,true(size(Psy_pred_ratio))});
    end
    function f9p2(debug)      % Others 1. Targ first vs Target Last (All SU + MU)
        if debug;  dbstack;  keyboard;  end
        
        %%
        
        j = 1;
        selectCells = (xls_num{1}(:,header.Units_RealSU) >=0 ) & (xls_num{1}(:,header.HD_TargFirst)~=0) & (~isnan(ChoiceDiv_All{j}(:,1,k)));
        
        figure(2099+figN); clf; figN = figN+1;
        
        for k = 1:3
            ys = nanmean(ChoiceDiv_All{j}(selectCells,:,k));
            errors = nanstd(ChoiceDiv_All{j}(selectCells,:,k))/sqrt(sum(selectCells));
            %             h1 = shadedErrorBar(rate_ts{j},ys,errors * 0,{'Color',colors(k,:)},transparent);
            
            h1.mainLine  = plot(rate_ts{j},ys,'Color',colors(k,:));
            set(h1.mainLine,'LineWidth',3);
            
            hold on;
        end
        
        str1 = sprintf('Choice divergence (All MU & SU, N = %g)',sum(selectCells));
        
        % for tt = 1:3
        %     plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
        % end
        % plot(xlim, [0 0],'k--');
        % SetFigure(); axis tight;
        
        
        %%% Choice divergence: Targ last
        
        selectCells = (xls_num{1}(:,header.HD_TargFirst)==0) & (~isnan(ChoiceDiv_All{j}(:,1,k)));
        
        % figure(2099+figN); clf; figN = figN+1;
        
        for k = 1:3
            ys = mean(ChoiceDiv_All{j}(selectCells,:,k));
            errors = std(ChoiceDiv_All{j}(selectCells,:,k))/sqrt(sum(selectCells));
            
            %             h2 = shadedErrorBar(rate_ts{j},ys,errors * 0,{'Linestyle','-.','Color',colors(k,:) * 0.5 + [.5 .5 .5] },transparent);
            
            h2.mainLine = plot(rate_ts{j},ys,'Linestyle',':','Color',colors(k,:) );
            set(h2.mainLine,'LineWidth',3)
            hold on;
        end
        
        str2 = sprintf('Choice divergence (Target Last, N = %g)',sum(selectCells));
        axis tight;
        
        for tt = 1:3
            plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
        end
        
        SetFigure();
        
        l = legend([h1.mainLine h2.mainLine],str1,str2);
        set(l,'FontSize',10,'box','off','location','best');
        
        plot(xlim, [0 0],'k--');
        
        %}
        
        %{
    %% Mean divergence curve and its derivation
    % figure(2099+figN); clf; figN = figN+1;
    %
    % mean_div = smooth(mean(ys_CD(1,:),1),5,'rloess');
    % plot(mean(ys_CD(1,:),1),'o'); hold
    % plot(mean_div,'o');
    % figure(2099+figN); clf; figN = figN+1;
    % plot(rate_ts{j}(1:end-1),diff(mean_div));
    
        %}
    end

    function f9p9(debug)      % Test
        if debug;  dbstack;  keyboard;  end
        
        
        weights_normalized_PSTH
        keyboard;
    end

%     keyboard;
%}
% catch err
%     err
%     opentoline(err.stack(1).file,err.stack(1).line);
%     keyboard
% end

    function Weights_Property_Correlation(weights,txt,select)  % Compare different population weights with cell properties
        % weights: [choice weight, modality weight]
        
        % ===========   1. Choice vs modality weights ============
        
        if size(weights,2) == 2
            h = LinearCorrelation({
                (weights(:,1));
                },...
                {
                weights(:,2) ;
                },...
                'CombinedIndex',[],...
                'Ylabel',txt{2},'Xlabel',txt{1},...
                'FaceColors',{'k'},'Markers',{'o'},...
                'LineStyles',{'k-'},'MarkerSize',12,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
            
            % Annotate tcells
            plot(weights(select_tcells(select),1),weights(select_tcells(select),2),'+','markersize',16,'color','r','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            set([gca h.group.dots],'ButtonDownFcn',{@Show_individual_cell,h.group.dots,select});
            
        end
        
        % ===========   2. Correlations =============
        % ---------  Choice weights and Choice preference
        tt = 1;
        k = 1;
        
        set(figure(figN),'name','Choice weights and Choice preference','pos',[17 514 1151 449]);
        
        cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
        
        h = LinearCorrelation({
            (weights(~cpref_sig,1));
            (weights(cpref_sig,1));
            },...
            {
            abs(Choice_pref_all(k,~cpref_sig ,tt)) ;
            abs(Choice_pref_all(k,cpref_sig,tt)) ;
            },...
            'CombinedIndex',[3],...
            'Ylabel','abs(Choice preference)','Xlabel',txt{1},...
            'FaceColors',{'none','k'},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',figN,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(weights(select_tcells(select),1),abs(Choice_pref_all(k,select_tcells(select),tt)),'+','markersize',16,'color','r','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(weights(:,1),abs(Choice_pref_all(k,:,tt)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
        
        % ---------  Memsac DDI and Choice Weight
        monkeys = xls_num{1}(:,header.Monkey);  % Temporarily here. HH20150914
        monkey1 = monkeys == 5;
        monkey2 = monkeys == 10;
        
%         h = LinearCorrelation({
%             MemSac_indicator(select);
%             },...
%             {
%             (weights(:,1));
%             },...
%             'CombinedIndex',[],...
%             'Xlabel',MemSac_indicator_txt,'Ylabel',txt{1},...
%             'FaceColors',{'k'},'Markers',{'o'},'Markers',{'o'},...
%             'LineStyles',{'k-'},'MarkerSize',12,...
%             'figN',figN,'XHist',20,'YHist',20,...
%             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
%         
%         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
%         
%         % Annotate tcells
%         plot(MemSac_indicator(select_tcells),weights(select_tcells(select),1),...
%             '+','markersize',16,'color','r','linew',2);
%         
%         % Show individual cell selected from the figure. HH20150424
%         h_line = plot(MemSac_indicator(select),(weights(:,1)),'visible','off'); hold on;
%         set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select});
        
        h = LinearCorrelation({
            MemSac_indicator(select & monkey1);
            MemSac_indicator(select & monkey2);
            },...
            {
            (weights(monkey1(select),1));
            (weights(monkey2(select),1));
            },...
            'CombinedIndex',[3],...
            'Xlabel',MemSac_indicator_txt,'Ylabel',txt{1},...
            'FaceColors',{'k'},'Markers',{'o','^'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',figN,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
        delete([h.group(1:2).line]);
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),weights(select_tcells(select),1),...
            '+','markersize',16,'color','r','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(MemSac_indicator(select),(weights(:,1)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select});        
        
        % ------------ Modality preference and Modality weight
        
        if size(weights,2) == 2
            
            tt = 1;
            k = 1; % Vis-vest
            
            set(figure(figN),'name','Modality weight vs Modality preference','pos',[17 514 1151 449]);
            
            mpref_sig = Modality_pref_p_value_all(1,:,tt)' < 0.05;
            
            h = LinearCorrelation({
                (Modality_pref_all(k,~mpref_sig & monkey1(select),tt)) ;
                (Modality_pref_all(k,~mpref_sig & monkey2(select),tt)) ;
                (Modality_pref_all(k,mpref_sig & monkey1(select),tt)) ;
                (Modality_pref_all(k,mpref_sig & monkey2(select),tt)) ;
                },...
                {
                (weights(~mpref_sig & monkey1(select),2));
                (weights(~mpref_sig & monkey2(select),2));
                (weights(mpref_sig & monkey1(select),2));
                (weights(mpref_sig & monkey2(select),2));
                },...
                'CombinedIndex',15,...
                'Xlabel','Modality preference (vis - vest)','Ylabel',txt{2},...
                'FaceColors',{'none','none','k','k'},'Markers',{'o','^'},...
                'LineStyles',{'k-'},'MarkerSize',12,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
            delete([h.group(1:4).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
            
            % Annotate tcells
            plot((Modality_pref_all(k,select_tcells(select),tt)), weights(select_tcells(select),2),...
                '+','markersize',16,'color','r','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot((Modality_pref_all(k,:,tt)),weights(:,2),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select});
        end
        
    end
    function h = Weighted_sum_PSTH(weights,txt,selects)     % Plot weighted sum of firing rates
        
        % ------  PSTH correct only, all choices -----------
        set(figure(figN),'pos',[95 95 646 544]); clf; figN = figN + 1;
        
        if ~iscell(weights)  ;  weights = {weights};     end
        if ~iscell(selects) ;  selects = {selects};   end
        
        for pp = 1:length(weights)
            
            % Now this function can receive all bootstrap weights for computing error bars. @HH20150525
            nbootstrap = size(weights{1},2);
            
            for j = 1:2
                PSTH_projected{j} = nan(nbootstrap,size(PSTH_all_raw{j},2),size(PSTH_all_raw{j},3));

                for bb = 1:nbootstrap 
                    proj_this = weights{pp}(:,bb) / sum(weights{pp}(:,bb)); % Normalization to let the firing rate make sense

                    % Projections
                    for kk = 1:size(PSTH_all_raw{j},3)
                        PSTH_projected{j}(bb,:,kk) = proj_this'*(PSTH_all_raw{j}(logical(selects{pp}),:,kk));
                    end
                end
            end
            
            h_subplot = subplot(1,length(weights),pp);
            
            % Note here the errorbars should be STD instead of SEM. (Not independent sampling, but bootstrapping)
            h_series = SeriesComparison({PSTH_projected{1} PSTH_projected{2}},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',mat2cell(reshape(repmat(colors,1,2)',3,[])',ones(stim_type_num*2,1)),'LineStyles',{'-','--'},...
                'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Weighted sum of firing','axes',h_subplot);
            hold on;    legend off;     axis tight;
            xlim([min(xlim) + 200 max(xlim)-300]);   % ylim([0.1 .7]);
                      
%             % Time markers
%             for tt = 1:3
%                 plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
%             end
%             
%             xlim([rate_ts{j}(10) rate_ts{j}(end-10)]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{1}(1), min(ylim) + Gauss_vel(:,2) * range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
            h.fig = gcf;
            h.projected_all_allbootstrap{pp} = PSTH_projected;
            h.projected_all_mean{pp} = h_series.means;
            
            title([txt{1} ', n = ' num2str(sum(selects{1}))]);
            
        end
        
        SetFigure(15);
        
        % ---------- Different angles ----------
        
        colors_angles = colormap(gray);
        colors_angles = colors_angles(round(linspace(length(colors_angles)-15,1,5)),:);
        colors_angles = reshape(repmat(colors_angles,1,2)',3,[])';
        colors_angles = mat2cell(colors_angles,ones(10,1));
        
        set(figure(figN),'name',[txt{1} ', n = ' num2str(sum(selects{1}))],'pos',[795 78 858 879]); clf; figN = figN + 1;
        h_subplot = tight_subplot(stim_type_num,2,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
        
        for k = 1:stim_type_num
            
            for j = 1:2
                
                h.projected_angles_allbootstrap{j} = nan(nbootstrap,size(PSTH_correct_angles_raw{j},2),size(PSTH_correct_angles_raw{j},3),size(PSTH_correct_angles_raw{j},4));
                h.projected_angles_diff{j} = nan(nbootstrap,size(PSTH_correct_angles_raw{j},2),size(PSTH_correct_angles_raw{j},3)/2,size(PSTH_correct_angles_raw{j},4));
                
                for bb = 1:nbootstrap
                    
                    % ------- Different angles -------
                    % Deal with nans of 0 heading for some cells
                    PSTH_correct_angles_raw_this = PSTH_correct_angles_raw{j}(selects{1},:,:,k);
                    PSTH_correct_angles_raw_this(isnan(PSTH_correct_angles_raw_this)) = 0;
                    
                    % Weighted sum with choice weight
                    yyy = reshape(reshape(PSTH_correct_angles_raw_this,size(PSTH_correct_angles_raw_this,1),[])' ...
                        * weights{1}(:,bb), size(PSTH_correct_angles_raw_this,2),[]);
                    h.projected_angles_allbootstrap{j}(bb,:,:,k) = yyy;
                    
                    % -------- Difference ---------
                    yyy_diff = yyy(:,1:2:end) - yyy(:,2:2:end);
                    h.projected_angles_diff{j}(bb,:,:,k) = yyy_diff;
                    
                end
            end
            
            % ------ Plotting ------
%             SeriesComparison({shiftdim(h.projected_angles{1}(:,:,k),-1) shiftdim(h.projected_angles{2}(:,:,k),-1)},...
%                 {rate_ts{1} rate_ts{2} time_markers},...
%                 'Colors',colors_angles,'LineStyles',{'-','--'},...
%                 'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k));
            
            h_series = SeriesComparison({h.projected_angles_allbootstrap{1}(:,:,:,k) h.projected_angles_allbootstrap{2}(:,:,:,k)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',colors_angles,'LineStyles',{'-','--'},...
                'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k));
            
            h.projected_angles_mean{1}(:,:,k) = h_series.means{1};
            h.projected_angles_mean{2}(:,:,k) = h_series.means{2};
            h.projected_angles_sem{1}(:,:,k) = h_series.errors{1};
            h.projected_angles_sem{2}(:,:,k) = h_series.errors{2};
            
            if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            if k == 1
                title([txt{1} ', n = ' num2str(sum(selects{1}))]);
            end
            
            % Gaussian vel
            axis tight;
            plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            xlim([min(xlim) + 200 max(xlim)-300]);   % ylim([0.1 .7]);
            legend off;
            
            % ------ Plotting diff ------
            
            SeriesComparison({h.projected_angles_diff{1}(:,:,:,k) h.projected_angles_diff{2}(:,:,:,k)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',colors_angles(1:2:end),'LineStyles',{'-'},...
                'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Diff','axes',h_subplot(k+stim_type_num));
            
            if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            legend off;
            
            axis tight;
            
            % Gaussian vel
            axis tight;
            plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            xlim([min(xlim) + 200 max(xlim)-300]);   % ylim([0.1 .7]);
            
        end
        SetFigure(15);
        
    end
    function Show_individual_cell(~,~,h_line, select_for_this, couples)    % Show individual cell selected from the figure. @HH20150424
        
        if nargin < 5
            couples = 1; % HH20150515 Show coupled cells in the figure
        end
        
        persistent h_marker;
        
        allX = get(h_line,'xData');
        allY = get(h_line,'yData');
        n_single_group = length(allX)/couples;
        
        % ------- Recover the cell number -------
        
        if ismember('control',get(gcf,'currentModifier'))  % Select cell from file name and show in the figure. @HH20150527
            fileN = input('Which cell do you want from the figure?    ','s');
            available_cells = find(select_for_this);
            
            nn = 1; iffind = [];
            while nn <= length(available_cells) && isempty(iffind) % Find the first match
                iffind = strfind(group_result(available_cells(nn)).cellID{1},fileN);
                nn = nn + 1;
            end
            
            if ~isempty(iffind)
                ind = nn - 1;
            else
                fprintf('Are you kidding...?\n');
                return
            end
            
        else   % Select cell from figure
            pos = get(gca,'currentPoint'); posX = pos(1,1); posY = pos(1,2);
            [min_dis,ind] = min(abs(((posX-allX)/range(xlim)).^2+((posY-allY)/range(ylim)).^2));
            if min_dis > (range(xlim)^2+range(ylim)^2)/100 +inf ; return; end
            
            ind = mod(ind-1,n_single_group) + 1; % Deal with coupled cells
        end
        
        real_cell_no = find(cumsum(select_for_this)==ind,1);
        
        % Plotting
        if ~isempty(h_marker) ; try delete(h_marker); catch ; end ;end
        
        all_inds = mod(1:length(allX),n_single_group) == mod(ind,n_single_group); % Deal with coupled cells
        h_marker = plot(allX(all_inds),allY(all_inds),'x','color','m','markersize',15,'linew',3);
        
        % PSTH in HD
        j_this = 1;
        
        ys_this{1} = group_result(real_cell_no).mat_raw_PSTH.PSTH{1,1,1}.ys';
        ys_this{2} = group_result(real_cell_no).mat_raw_PSTH.PSTH{2,1,1}.ys';
        
        SeriesComparison({shiftdim(ys_this{1},-1) shiftdim(ys_this{2},-1)},...
            {rate_ts{1} rate_ts{2} time_markers},...
            'Colors',mat2cell(reshape(repmat(colors,1,2)',3,[])',ones(stim_type_num*2,1)),'LineStyles',{'-','--'},'figN',1463);
        
        axis tight; legend off;
%         set(gcf,'pos',[-517 453 511 393]);
                set(gcf,'pos',[6 568 511 393]);
        
        for t = 1:3
            plot([1 1] * time_markers{j_this}(1,t),ylim,'k','linestyle',marker_for_time_markers{j_this}{t},'linew',1.5);
        end
        title(sprintf('#%g, %s',...
            real_cell_no,group_result(real_cell_no).cellID{1}(9:end)));
        xlabel('Time to saccade onset (ms)');      ylabel('Firing rate');
        xlim(xlim * 0.9);
        SetFigure(10);
        
        % PSTH in Memsac
        Plot_memsac(real_cell_no);
    end
    function Plot_memsac(cell_no)    % Plot associated mem-sac traces @HH20150426
        
%         set(figure(1462),'pos',[-1243 452 712 394]); clf;
                    set(figure(1462),'pos',[8 81 712 394]); clf;
        
        try
            resp_mean = group_result(cell_no).mat_raw_MemSac.resp_mean(2:end);
            resp_err = group_result(cell_no).mat_raw_MemSac.resp_err(2:end);
            p = group_result(cell_no).mat_raw_MemSac.p(2:end);
            vectSum = group_result(cell_no).mat_raw_MemSac.vectSum(2:end);
            temporal_Slice = group_result(cell_no).mat_raw_MemSac.temporal_Slice(2:end,:);
            align_offsets = group_result(cell_no).mat_raw_MemSac.align_offsets;
            align_markers = group_result(cell_no).mat_raw_MemSac.align_markers;  % Desired markers: target onset & target offset & sac onset
            unique_heading = group_result(cell_no).mat_raw_MemSac.unique_heading;
            PREF_target_location = group_result(cell_no).mat_raw_PSTH.PREF_target_location;
            MemSac_interp_PSTH = group_result(cell_no).MemSac_interp_PSTH;
            MemSac_interp_locations = group_result(cell_no).MemSac_interp_locations;
             
            axis_typing = axes('Position',[0.032 0.672 0.459 0.296]);
            axis_polar = axes('Position',[0.039 0.096 0.345 0.586]);
            axis_psth = axes('Position',[0.521 0.17 0.426 0.728]);
            
            % Polar plot
            [~, polarOrder] = sort(max(cell2mat(resp_mean'),[],2),1,'descend'); % The largest goes first in the polar plot
            SetFigure(13);
            
            for sliceN = polarOrder(:)'
                axes(axis_polar);
                
                if strcmp(temporal_Slice{sliceN,5},''); continue; end
                
                if p(sliceN) < 0.05
                    sigMarker = '-'; wid = 2;
                else
                    sigMarker = '-'; wid = 1;
                end;
                
%                 polarwitherrorbar([unique_heading/180*pi; unique_heading(1)/180*pi], ...
%                     [resp_mean{sliceN}'; resp_mean{sliceN}(1)], (p(sliceN) < 0.05) * [resp_err{sliceN}' ; resp_err{sliceN}(1)],[temporal_Slice{sliceN,5} sigMarker],wid);

                set(polar([unique_heading/180*pi; unique_heading(1)/180*pi], ...
                    [resp_mean{sliceN}'; resp_mean{sliceN}(1)],[temporal_Slice{sliceN,5} sigMarker]),'linewidth',2);

                hold on;
                
                if p(sliceN) < 0.05
                    %         h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,[0 vectAmp(sliceN)],[temporal_Slice{sliceN,5} '-']);
                    h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,[max(cell2mat(resp_mean))*0.9 max(cell2mat(resp_mean))*1.3],[temporal_Slice{sliceN,5} '-']);
                    set(h_p,'LineWidth',3);
                else
                    h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,[max(cell2mat(resp_mean))*0.9 max(cell2mat(resp_mean))*1.3],[temporal_Slice{sliceN,5} '-']);
                    set(h_p,'LineWidth',1);
                end
                %
                axes(axis_typing); axis off;
                text(0,- sliceN * 0.15+1.1,sprintf('%s:  \\itp_ = %4.2g',temporal_Slice{sliceN,4},p(sliceN)),'color',temporal_Slice{sliceN,5},...
                    'FontSize',14 * (p(sliceN) < 0.05) + 7 * (p(sliceN) >= 0.05 || isnan(p(sliceN))));
                drawnow;
            end
            
            % ----- PSTH plot -----
            
            % The one which is nearest to the vect_sum of pre-S
            [~,pref] = min(abs(vectSum(3) - MemSac_interp_locations));
            pref = mod(pref-1,length(MemSac_interp_locations)-1)+1;
            null = mod(pref + (length(MemSac_interp_locations)-1)/2 -1, length(MemSac_interp_locations)-1)+1; % The opposite position
            
            ts = group_result(cell_no).mat_raw_MemSac.t_centers{3};
            ys = MemSac_interp_PSTH(:,[pref null]);

            axes(axis_psth); hold on;
            SeriesComparison(shiftdim(ys,-1),ts,'Colors',{'r','r'},'LineStyles',{'-','--'},'axes',axis_psth);
            
            % The two which are nearest to actual PREF location. @HH20150524
            [~,pref] = min(abs(PREF_target_location - MemSac_interp_locations));
            pref = mod(pref-1,length(MemSac_interp_locations)-1)+1;
            null = mod(pref + (length(MemSac_interp_locations)-1)/2 -1, length(MemSac_interp_locations)-1)+1; % The opposite position

            axes(axis_polar); 
            set(polar([MemSac_interp_locations(pref) MemSac_interp_locations(pref)]/180*pi,[0 max(cell2mat(resp_mean))*1.3],'k-'),'linew',2);
            set(polar([MemSac_interp_locations(null) MemSac_interp_locations(null)]/180*pi,[0 max(cell2mat(resp_mean))*1.3],'k--'),'linew',2);
            
            ys = MemSac_interp_PSTH(:,[pref null]);

            axes(axis_psth);
            SeriesComparison(shiftdim(ys,-1),ts,'Colors',{'k','k'},'LineStyles',{'-','--'},'axes',axis_psth);
            
            
            xlabel('Time to saccade onset (ms)');   ylabel('Firing rate');
            plot([0 0],ylim,'k-','linew',1.5);
            axis tight;
            legend off;
            title(sprintf('PREF: %g, %s',PREF_target_location,group_result(cell_no).cellID{2}(30:end)));

            % Annotating temporal slice windows

            for sliceN = 1:size(temporal_Slice,1)
                if temporal_Slice{sliceN,3} == 7 % Precise windows
                    plot([temporal_Slice{sliceN,1} temporal_Slice{sliceN,2}],[max(ylim) max(ylim)],temporal_Slice{sliceN,5},'linew',15);
                else  % Windows that is not so precise
                    % Mean shift
                    meanShift = mean(align_offsets(:,align_markers == temporal_Slice{sliceN,3})-align_offsets(:,align_markers==7),1);
                    plot( [meanShift meanShift],ylim,'k--','LineWidth',1);
                    plot(meanShift+[temporal_Slice{sliceN,1} temporal_Slice{sliceN,2}],[max(ylim) max(ylim)],temporal_Slice{sliceN,5},'linew',15);
                end
            end
            
        catch
        end
    end


end

