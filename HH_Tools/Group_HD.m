%%  LIP HD Pooled Data
% Begin @ HH 201406

function function_handles = Group_HD(XlsData)
% Matlab version
matlab_ver = version;
ver_num = version('-release');
ver_num = str2num(ver_num(1:end-1));

%% Constants
LEFT = 1;
RIGHT = 2;
tmp1 = []; tmp2 = []; tmp3 = []; tmp4 = []; tmp5  = [];

%% Get data
num = XlsData.num;
txt = XlsData.txt;
raw = XlsData.raw;
header = XlsData.header;
stim_type_num = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch address
mat_address = {
    % Major protocol goes here (Address, Suffix)
    
    %     'Z:\Data\Tempo\Batch\20180608_HD_all_IONCluster_LightWeight\','PSTH';  % Git cfb3e647e93
%     'Z:\Data\Tempo\Batch\20180607_HD_all_IONCluster_CDPerm_CP10ms\','PSTH';  % Git 7cf56f23a8   (j=2, 270 time bins) 
    % 'Z:\Data\Tempo\Batch\20180619_HD_all_IONCluster_LightWeight+fisherSimple\','PSTH'; % Git ec5859b24933a                                                            
    'Z:\Data\Tempo\Batch\20180711_HD_all_IONCluster_LightWeight+fisherSimple+PSTHAll','PSTH'; % Git e86e4f76 
    
    %     'Z:\Data\Tempo\Batch\20160918_HD_allAreas_m5_m10_Smooth50ms_NewChoicePref','PSTH';
    %     'Z:\Data\Tempo\Batch\20160908_HD_allAreas_m5_m10_Smooth50ms','PSTH';
    %     'Z:\Data\Tempo\Batch\20160907_HD_allAreas_m5_m10_Smooth20ms','PSTH';
    %     'Z:\Data\Tempo\Batch\20150725_HD_allAreas_m5_m10','PSTH';
    %     'Z:\Data\Tempo\Batch\20150418_LIP_HD_m5_m10_modifiedBatch','PSTH';
    %     'Z:\Data\Tempo\Batch\20150411_LIP_HD_m5','PSTH';
    %     'Z:\Data\Tempo\Batch\20141118_LIP_Decision_WithAUC_ResultPSTHOrderChanged','PSTH';
    
    % Associative protocols
    'Z:\Data\Tempo\Batch\20180608_Memsac_all_IONCluster','MemSac';
    
    %     'Z:\Data\Tempo\Batch\20150725_BP_allAreas_m5_m10','MemSac';
    %     'Z:\Data\Tempo\Batch\20150411_LIP_memSac_m5_m10','MemSac';
    };

% Global Mask for each protocol (for XLS)
mask_all = {
    (strcmp(txt(:,header.Protocol),'HD') | strcmp(txt(:,header.Protocol),'HD_dt')) & (strcmpi(txt(:,header.Area),'LIP') | strcmpi(txt(:,header.Area),'LIP-V')) ...
    & (num(:,header.HD_rep) >= 8) & (num(:,header.Units_RealSU) == 1) & (num(:,header.Monkey) == 5 | num(:,header.Monkey) == 10)...
    & (num(:,header.HD_TargFirst) == 1);  % To speed up loading. To include HD_dt. HH20160918
    
    (strcmp(txt(:,header.Protocol),'MemSac') | strcmp(txt(:,header.Protocol),'2MemSac')) & (strcmpi(txt(:,header.Area),'LIP') | strcmpi(txt(:,header.Area),'LIP-V'));
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
        tolerance_enable = 1;
        depthTol = 0; % Depth tolerance
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
                        
                        % Deal with 2point memsac. @HH20160906
                        if_2pt_memsac = strcmp('2MemSac',xls_txt{pp}(match_i,header.Protocol));
                        
                        if sum(if_2pt_memsac)  % We have 2pt memsac for this cell (almost all Messi cells but none for Polo)
                            % Note now the file_i has a length of 2.
                            % fprintf('Memsac and 2p-Memsac detected\n');
                            file_i(1) = match_i(find(if_2pt_memsac==0,1,'first')); % The first 8pt Memsac
                            file_i(2) = match_i(find(if_2pt_memsac==1,1,'first')); % The first 2pt Memsac
                        else  % No 2pt memsac, real duplication (maybe different eccentricities)
                            fprintf('More than one match have been found for %s,\n Major ID = %s\n', xls_txt{1}{major_i,header.FileNo},cell_info{1}{major_i});
                            disp(cell_info{pp}(match_i));
                            %                 file_choose = input('Which one do you want?');
                            %                 file_i = match_i(file_choose);
                            file_i = match_i(1);  % The first appearance by default.
                        end
                        
                        %                 keyboard;
                    else % Special cases for inexact match (MU-SU)
                        if ~tolerance_enable || ~strcmp(xls_txt{1}(major_i,header.Note),'MU=SU')
                            fprintf('No exact matched files for %s, ID = %s\n', xls_txt{1}{major_i,header.FileNo},cell_info{1}{major_i});
                            not_match = not_match + 1;
                            file_i = NaN;
                        else  % Unit tolerance enabled and "MU=SU"   HH20160920
                            match_i = find(strncmp(cell_info{1}{major_i},cell_info{pp},29));
                            
                            if length(match_i) == 1
                                fprintf('Unit tolerance found:\n   %s (major)\n   %s\n\n', cell_info{1}{major_i}, cell_info{pp}{match_i});
                                file_i = match_i;
                            elseif length(match_i) > 1
                                % Deal with 2point memsac. @HH20160906
                                if_2pt_memsac = strcmp('2MemSac',xls_txt{pp}(match_i,header.Protocol));
                                
                                if sum(if_2pt_memsac)  % We have 2pt memsac for this cell (almost all Messi cells but none for Polo)
                                    % Note now the file_i has a length of 2.
                                    % fprintf('Memsac and 2p-Memsac detected\n');
                                    file_i(1) = match_i(find(if_2pt_memsac==0,1,'first')); % The first 8pt Memsac
                                    file_i(2) = match_i(find(if_2pt_memsac==1,1,'first')); % The first 2pt Memsac
                                    fprintf('Unit tolerance found (and 2pt memsac):\n   %s (major)\n   %s\n\n   %s\n\n', cell_info{1}{major_i}, cell_info{pp}{match_i(1)},cell_info{pp}{match_i(2)});
                                else  % No 2pt memsac, real duplication (maybe different eccentricities)
                                    fprintf('More than one match have been found for %s,\n Major ID = %s\n', xls_txt{1}{major_i,header.FileNo},cell_info{1}{major_i});
                                    disp(cell_info{pp}(match_i));
                                    %                 file_choose = input('Which one do you want?');
                                    %                 file_i = match_i(file_choose);
                                    file_i = match_i(1);  % The first appearance by default.
                                end

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
                        for ii = 1:length(file_i) % For 2pt_memsac. @HH20160906
                            mat_file_name = sprintf('%s_%g',xls_txt{pp}{file_i(ii),header.FileNo},xls_num{pp}(file_i(ii),header.Chan1));
                            mat_file_fullname = [mat_address{pp,1} '\' mat_file_name '_' mat_address{pp,2}];
                            
                            raw = load(mat_file_fullname);
                            
                            group_result(major_i).cellID{pp}{ii} = cell_info{pp}{file_i(ii)};
                            group_result(major_i).fileID = xls_txt{pp}{file_i(ii),header.FileNo}; % For behavior analysis. HH20190702

                            
                            if strcmp('2MemSac',xls_txt{pp}(file_i(ii),header.Protocol)) % 2Mem-sac
                                group_result(major_i).(['mat_raw_2pt' mat_address{pp,2}]) = raw.result;  % Dynamic structure
                            else
                                group_result(major_i).(['mat_raw_' mat_address{pp,2}]) = raw.result;  % Dynamic structure
                            end
                        end
                    catch err
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
ALL_CHOICE = 1; CORRECT_ANGLE = 2; CHOICE_DIFFICULT = 3; OUTCOME = 4; WRONG_ANGLE = 5; CORRECTNWRONG_ANGLE = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
representative_cell = 146; % Last representative cell (HD, not HD_dt) for Messi. HH20160919

for i = 1:length(group_result)
    
    group_result(i).repN = group_result(i).mat_raw_PSTH.repetitionN;
    unique_stim_type = group_result(i).mat_raw_PSTH.unique_stim_type;  % In case we don't have all three conditions
    group_result(i).length_unique_stim_type = length(unique_stim_type);
    group_result(i).unique_heading = group_result(i).mat_raw_PSTH.CP{1,1}.raw_CP_result{1}.Neu_tuning(:,1);
    
    %     if isempty(group_result(i).mat_raw_PSTH)
    %         continue
    %     end
    
    % ------------ Patch in HD_dt cells. HH20160919 -----------
    % To include 18 HD_dt typical cells into my original HD tasks, I align PSTH using the center of the Gaussian movement instead of the
    % actual stim-on for those cells. So the time_alignment{1} should be moved 125 ms forward (in HD_dt task, stim duration = 1750 ms)
    % So now, VSTIM_ON_CD   means --> 1500 ms Gaussian begins (125 ms after the real VSTIM_ON_CD)
    %         VSTIM_OFF_CD  means --> 1500 ms Gaussian ends (125 ms before the real VSTIM_OFF_CD)
    % Here I just mark these cells out.
    if isfield(group_result(i).mat_raw_PSTH,'HD_dt_patched')
        group_result(i).HD_dt_patched = 1;
    else
        group_result(i).HD_dt_patched = 0;
    end
    
    % Time markers (stim on, stim off, sac on)
    group_result(i).time_marker{1} = [0 mean(group_result(i).mat_raw_PSTH.align_offsets_others{1})];
    group_result(i).time_marker{2} = [mean(group_result(i).mat_raw_PSTH.align_offsets_others{2}) 0];
    
    % Psychophysics
    group_result(i).Psy_para = nan(2,3);
    for stim_type = 1:3
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
        
        for stim_type = 1:3  % Always output three conditions
            
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
        
        for stim_type = 1:3  % Always output three conditions
            
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
        
        % Only align to saccade here
        group_result(i).MemSac_ts = group_result(i).mat_raw_MemSac.t_centers{3};
        group_result(i).MemSac_PSTH = reshape([group_result(i).mat_raw_MemSac.result_PSTH_anne_mean{3,:}],length(group_result(i).MemSac_ts),[]);

    end

    if ~isempty(group_result(i).mat_raw_2ptMemSac)
        group_result(i).TwoPtMemSac_p = group_result(i).mat_raw_2ptMemSac.p;
        group_result(i).TwoPtMemSac_DDI = group_result(i).mat_raw_2ptMemSac.DDI;
        group_result(i).TwoPtMemSac_vectSum = group_result(i).mat_raw_2ptMemSac.vectSum;
        group_result(i).TwoPtMemSac_AI = group_result(i).mat_raw_2ptMemSac.activityIndex;
        
        % Only align to saccade here
        group_result(i).TwoPtMemSac_ts = group_result(i).mat_raw_2ptMemSac.t_centers{3};
        group_result(i).TwoPtMemSac_PSTH = reshape([group_result(i).mat_raw_2ptMemSac.result_PSTH_anne_mean{3,:}],...
                                                    length(group_result(i).TwoPtMemSac_ts),[]);
  
    
    end

    %     catch
    %         fprintf('Error in preprocessing %g\n',i);
    %     end
    
    % 6) Spike waveform (finally!). HH20160920
    waveform = str2num(cell2mat(xls_txt{1}(i,header.meanWav)));
    % waveform = smooth(waveform,5)'; waveform = waveform - min(waveform); waveform = waveform./max(waveform);
    dt = 1/25; % 1k ms / 25 kHz
    interp_dt = 1/1000; % Spline interp to 1us. (Keven Johnston 2009 JNS)
    
    if ~isempty(waveform)
        
        waveform = spline(0:dt:dt*(length(waveform)-1),waveform,0:interp_dt:dt*(length(waveform)-1));
        
        waveform_t1 = find(waveform == 1); % Spike peak
        waveform_t_left = find(waveform == min(waveform(1:waveform_t1)),1); % Left trough
        waveform_t_right = find(waveform == min(waveform(waveform_t1:end)),1); % Right trough
        
        group_result(i).Waveform_peakToTrough = (waveform_t_right - waveform_t1)*interp_dt;
        group_result(i).Waveform_peak = waveform_t1;
        group_result(i).Waveform_trough1 = waveform_t_left;
        group_result(i).Waveform_trough2 = waveform_t_right;
        
        % Abnormal waveform (need to be flipped)
        if group_result(i).Waveform_peakToTrough > 0.8  % 0.8 ms was set manually by plotting the distribution
            waveform = 1-waveform;
            waveform_t1 = find(waveform == 1); % Spike peak
            waveform_t_left = find(waveform == min(waveform(1:waveform_t1)),1); % Left trough
            waveform_t_right = find(waveform == min(waveform(waveform_t1:end)),1); % Right trough
            
            group_result(i).Waveform_peakToTrough = (waveform_t_right - waveform_t1)*interp_dt;
            group_result(i).Waveform_peak = waveform_t1;
            group_result(i).Waveform_trough1 = waveform_t_left;
            group_result(i).Waveform_trough2 = waveform_t_right;
            fprintf('Waveform flipped: %s\n',group_result(i).cellID{1}{1});
        end
        
        group_result(i).Waveform_broad = group_result(i).Waveform_peakToTrough > 0.35; % Set manually by plotting the distribution
        
        % Align to waveform_1
        %{
        figure(919); hold on;
        plot(((1:length(waveform))-waveform_t_left)*dt,waveform);
        plot((waveform_t1-waveform_t_left)*dt,waveform(waveform_t1),'or');
        plot((waveform_t_right-waveform_t_left)*dt,waveform(waveform_t_right),'og');
        xlabel('ms');
        %}
    else
        fprintf('No waveform: %s\n', group_result(i).cellID{1}{1});
        group_result(i).Waveform_peakToTrough = nan;
        group_result(i).Waveform_broad = nan;
    end
    
    group_result(i).Waveform = waveform;

    % 7) Cell location.
    
end

% --- Waveform width distribution ---
% This is to determine the threshold of narrow/wide spikes. HH20160920
% figure(); hist([group_result.Waveform_peakToTrough],10);

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
rate_ts = group_result(representative_cell).rate_ts;
CP_ts = group_result(representative_cell).CP_ts;
memsac_ts = group_result(representative_cell).MemSac_ts;

% Initialization

% j-sensitive
for j = 1:2
    % For normalized PSTH plotting (weighted by dynamic ranges for each neuron)
    PSTH_all_Norm{j} = NaN(N,length(rate_ts{j}),6);
    PSTH_correct_angles_Norm{j} = NaN(N,length(rate_ts{j}),1+length(group_result(representative_cell).unique_heading),3);  % Two zero headings
    PSTH_correctNwrong_angles_Norm{j} = NaN(N,length(rate_ts{j}),length(group_result(representative_cell).unique_heading),3);  % One zero heading
    PSTH_outcomes_Norm{j} = NaN(N,length(rate_ts{j}),4,3);
    PSTH_wrong_angles_Norm{j} = NaN(N,length(rate_ts{j}),length(group_result(representative_cell).unique_heading)-1,3); % @HH20150523
    
    % Also pack raw data for different ways of weighted sum PSTH plotting (weighted by SVM / targeted dimensionality reduction, etc.)
    PSTH_all_raw{j} = NaN(N,length(rate_ts{j}),6);
    PSTH_correct_angles_raw{j} = NaN(N,length(rate_ts{j}),1+length(group_result(representative_cell).unique_heading),3);
    PSTH_correctNwrong_angles_raw{j} = NaN(N,length(rate_ts{j}),length(group_result(representative_cell).unique_heading),3);  % One zero heading
    PSTH_outcomes_raw{j} = NaN(N,length(rate_ts{j}),4,3);
    PSTH_wrong_angles_raw{j} = NaN(N,length(rate_ts{j}),length(group_result(representative_cell).unique_heading)-1,3); % @HH20150523
    PSTH_hard_easy_raw_cellBtB4Bk{j} = NaN(N,length(rate_ts{j}),4,3); % HH20160905 For EI calculation. Cell by Time by 4 (diff, easy) by stimtype
    
    CP{j} = NaN(N,length(CP_ts{j}),3);
    
    ChoiceDiv_All{j} = NaN(N,length(rate_ts{j}),3);
    
    ChoiceDiv_All_perm{j}.std = NaN(N,length(rate_ts{j}),3);  % HH20180608
    ChoiceDiv_All_perm{j}.p = NaN(N,length(rate_ts{j}),3);

    ChoiceDiv_Easy{j} = NaN(N,length(rate_ts{j}),3);
    ChoiceDiv_Difficult{j} = NaN(N,length(rate_ts{j}),3);
    ChoiceDiv_EasyMinusDifficult{j} = NaN(N,length(rate_ts{j}),3);
    
    ChoiceDiv_ModDiffer{j} = NaN(N,length(rate_ts{j}),3);  % 3-1, 3-2, 1-2
    ModDiv_All{j} = NaN(N,length(rate_ts{j}),3); % HH20140415   % 2-1, 3-1, 3-2
end

% j-insensitive
group_PREF_target_location = NaN(N,1);
group_PREF_target_location_notebook = NaN(N,1); % From notebook @HH20160906
group_MemSac_DDI = NaN(N,6);
group_MemSac_ps = NaN(N,6);
group_MemSac_AI = NaN(N,6);
group_TwoPtMemSac_DDI = NaN(N,6);
group_TwoPtMemSac_ps = NaN(N,6);
group_TwoPtMemSac_AI = NaN(N,6);
group_MemSac_PREF_Null_DI = NaN(N,6);
group_MemSac_actual_DI = NaN(N,1); % DI of actuall target locations (from 8pt Mem). @HH20150524
group_MemSac_actual_DI_2pt = NaN(N,1); % DI (from 2pt Mem)
group_MemSac_PREFmNULL_LeftRight_PSTH = NaN(N,length(group_result(representative_cell).MemSac_ts));
group_MemSac_PSTH_AngDiff = NaN(N,6);

group_ChoicePreference_pvalue = reshape([group_result(:).ChoicePreference_pvalue]',3,[],3);  % I updated ChoicePref time window. HH20160918
group_ChoicePreference = reshape([group_result(:).ChoicePreference]',3,[],3);  % Stim, Cell No, Pre/Post. Updated time win HH20160918
% group_ChoicePrefernce_pvalue = squeeze(group_ChoicePrefernce_pvalue(3,:,:))'; % HH20160918. 3: stim-on to stim-off

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
        
        for stim_type = 1:3 % Stim_type check
            
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
                
                try
                    PSTH_correctNwrong_angles_raw_this =  group_result(i).mat_raw_PSTH.PSTH{j,CORRECTNWRONG_ANGLE,k}.ys;
                catch
                    PSTH_correctNwrong_angles_raw_this = nan;
                    if ~exist('ww','var')
                        warning('No "Correct + Wrong" in PSTH, you must have used the old data with the new Group_HD.m ...');
                        ww = 1;
                    end
                end
                PSTH_correctNwrong_angles_norm_this =  PSTH_correctNwrong_angles_raw_this - offset;
                PSTH_correctNwrong_angles_norm_this = PSTH_correctNwrong_angles_norm_this / gain;
                
                if size(PSTH_correct_angles_norm_this,1) == size(PSTH_correct_angles_Norm{j},3)
                    PSTH_correct_angles_Norm{j}(i,:,:,k) = PSTH_correct_angles_norm_this';
                    PSTH_correct_angles_raw{j}(i,:,:,k) = PSTH_correct_angles_raw_this';
                    
                    PSTH_correctNwrong_angles_Norm{j}(i,:,:,k) = PSTH_correctNwrong_angles_norm_this';
                    PSTH_correctNwrong_angles_raw{j}(i,:,:,k) = PSTH_correctNwrong_angles_raw_this';
                    
                elseif size(PSTH_correct_angles_norm_this,1) == size(PSTH_correct_angles_Norm{j},3) - 2 % Without zero heading
                    PSTH_correct_angles_Norm{j}(i,:,3:end,k) = PSTH_correct_angles_norm_this';
                    PSTH_correct_angles_raw{j}(i,:,3:end,k) = PSTH_correct_angles_raw_this';
                    
                    PSTH_correctNwrong_angles_Norm{j}(i,:,[1:fix(end/2) fix(end/2)+2:end],k) = PSTH_correctNwrong_angles_norm_this';
                    PSTH_correctNwrong_angles_raw{j}(i,:,[1:fix(end/2) fix(end/2)+2:end],k) = PSTH_correctNwrong_angles_raw_this';
                else
                    disp('No match PSTH_angles_norm...');
                end
                
                if group_result(i).mat_raw_PSTH.PREF == LEFT  % Note that PSTH{6} is grouped by raw heading, not flip to "RIGHT = PREF"
                    PSTH_correctNwrong_angles_Norm{j}(i,:,:,k) = fliplr(squeeze(PSTH_correctNwrong_angles_Norm{j}(i,:,:,k)));
                    PSTH_correctNwrong_angles_raw{j}(i,:,:,k) = fliplr(squeeze(PSTH_correctNwrong_angles_raw{j}(i,:,:,k)));
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
                
                % HH20160905 EI
                PSTH_hard_easy_raw_cellBtB4Bk{j}(i,:,:,k) = group_result(i).mat_raw_PSTH.PSTH{j,CHOICE_DIFFICULT,k}.ys';
            end
            
            % Stim type check already done
            ChoiceDiv_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_ALL{j}(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_All_perm{j}.std(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_ALL_perm{j}.std(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_All_perm{j}.p(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_ALL_perm{j}.p(stim_type,:),smooth_factor_for_divergence);

            
            ChoiceDiv_Easy{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Easy{j}(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_Difficult{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Difficult{j}(stim_type,:),smooth_factor_for_divergence);
            
            ChoiceDiv_EasyMinusDifficult{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Easy{j}(stim_type,:) ...
                - group_result(i).mat_raw_PSTH.ChoiceDivergence_Difficult{j}(stim_type,:),smooth_factor_for_divergence);
            
            
            CP{j}(i,:,stim_type) = group_result(i).CP{j}(stim_type,:);
            CP_p{j}(i,:,stim_type) = group_result(i).CP_p_value{j}(stim_type,:);
            
            ModDiv_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ModalityDivergence{j}(stim_type,:),smooth_factor_for_divergence);
            
        end
        
        ChoiceDiv_ModDiffer{j}(i,:,1) = smooth(ChoiceDiv_All{j}(i,:,3) - ChoiceDiv_All{j}(i,:,1),smooth_factor_for_divergence);
        ChoiceDiv_ModDiffer{j}(i,:,2) = smooth(ChoiceDiv_All{j}(i,:,3) - ChoiceDiv_All{j}(i,:,2),smooth_factor_for_divergence);
        ChoiceDiv_ModDiffer{j}(i,:,3) = smooth(ChoiceDiv_All{j}(i,:,1) - ChoiceDiv_All{j}(i,:,2),smooth_factor_for_divergence);
    end
    
    % j-insensitive
    
    % --- Find the PREF target location ---
    % (1) From the eye trace. See HeadingDis_cum_PSTH_HH.m
    tmp = group_result(i).mat_raw_PSTH.PREF_target_location; % [-90,270]. There was a bug:)
    group_PREF_target_location(i) = mod(tmp,360);  % Changed to [0,360] @HH20160906
    
    % (2) Use the notebook. Note I only record the right target (not necessarily the PREF target). HH20160906
    tmp = xls_num{1}(i,header.HD_TargAng); % [-90,90]
    if group_result(i).PREF_PSTH == LEFT  % Flip 180 degree
        tmp = tmp + 180;   % [-90,270]
    end
    group_PREF_target_location_notebook(i) = mod(tmp,360); % [0,360]

    
    % Mem-sac stuff (Note some are NaNs)
    if ~isempty(group_result(i).MemSac_DDI)
        % -------- Global memsac indicator -------
        group_MemSac_DDI(i,:) = group_result(i).MemSac_DDI;
        group_MemSac_ps(i,:) = group_result(i).MemSac_p;
        group_MemSac_AI(i,:) = group_result(i).MemSac_AI;
        
        % -------- Local memsac indicator --------
        % 1. Left and Right (Pref/null of HD task)
        group_MemSac_PREFmNULL_LeftRight_PSTH(i,:) = (group_result(i).MemSac_PSTH(:,1) - group_result(i).MemSac_PSTH(:,5)) *  sign((group_result(i).PREF_PSTH == 2)-0.5);  % (Right - Left)* Right is pref
        group_MemSac_PSTH_AngDiff(i,:) = mod(group_result(i).MemSac_vectSum - ((group_result(i).PREF_PSTH == 2)*0 + (group_result(i).PREF_PSTH == 1)*180),360);
        
        % 2. Mem's Pref_null DI. @HH20150524
        group_MemSac_PREF_Null_DI(i,:) = group_result(i).mat_raw_MemSac.PREF_NULL_DI;
        
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
        group_result(i).MemSac_interp_PSTH_alignVisON = group_result(i).mat_raw_MemSac.MemSac_interp_PSTH{1}; % Added align to Vis Off
        
        
        % --- Calculate Memsac DI using user defined period ----
        MemSac_actual_DI_period = [3];  % [3,4] % I use memory and pre period to calculate the actual DI
        
        MemSac_temporal_Slice = group_result(i).mat_raw_MemSac.temporal_Slice;
        MemSac_align_offsets = group_result(i).mat_raw_MemSac.align_offsets;
        MemSac_align_markers = group_result(i).mat_raw_MemSac.align_markers;  
        MemSac_ts = group_result(i).MemSac_ts;
        MemSac_actual_DI_time_ind = false(size(MemSac_ts));
        
        for p_ind = 1:length(MemSac_actual_DI_period)
            
            ppp = MemSac_actual_DI_period(p_ind);
            
            if MemSac_temporal_Slice{ppp,3} == 7  % Precise windows (because MemSac_ts has been aligned to 7 (saccade onset))
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
        
        % Use the interpolated data from 8p Memsac
        % Find actual DI for this cell
        [~,pref] = min(abs(group_PREF_target_location(i) - group_result(i).MemSac_interp_locations));
        pref = mod(pref-1,length(group_result(i).MemSac_interp_locations)-1)+1;
        null = mod(pref + (length(group_result(i).MemSac_interp_locations)-1)/2 -1, length(group_result(i).MemSac_interp_locations)-1)+1; % The opposite position
        
        MemSac_PREF_mean = mean(group_result(i).MemSac_interp_PSTH(MemSac_actual_DI_time_ind,pref),1);
        MemSac_NULL_mean = mean(group_result(i).MemSac_interp_PSTH(MemSac_actual_DI_time_ind,null),1);
        
        group_MemSac_actual_DI(i) = (MemSac_PREF_mean - MemSac_NULL_mean) / (MemSac_PREF_mean + MemSac_NULL_mean);
        
    end
    
    group_MemSac_PSTH_AngDiff(group_MemSac_PSTH_AngDiff > 180) = 360 - group_MemSac_PSTH_AngDiff(group_MemSac_PSTH_AngDiff > 180); % Ang diffs are less than 180
    
    % 4. 2-point Memsac  @HH20160906
    if ~isempty(group_result(i).mat_raw_2ptMemSac)
        group_TwoPtMemSac_DDI(i,:) = group_result(i).TwoPtMemSac_DDI;
        group_TwoPtMemSac_ps(i,:) = group_result(i).TwoPtMemSac_p;
        group_TwoPtMemSac_AI(i,:) = group_result(i).TwoPtMemSac_AI;

        MemSac_PREF_mean = mean(group_result(i).TwoPtMemSac_PSTH(MemSac_actual_DI_time_ind,3-group_result(i).PREF_PSTH),1);
        MemSac_NULL_mean = mean(group_result(i).TwoPtMemSac_PSTH(MemSac_actual_DI_time_ind,group_result(i).PREF_PSTH),1);
        group_MemSac_actual_DI_2pt(i) = (MemSac_PREF_mean - MemSac_NULL_mean) / (MemSac_PREF_mean + MemSac_NULL_mean);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- Choose Memsac indicator (affect all Memsac measures) -------------
% 1 = Background, 2 = LS, 3 = Mem, 4 = Pre, 5 = Co, 6 = Post

MemSac_indicator = mean(group_MemSac_DDI(:,[3]),2); % Global DDI
MemSac_indicator_p = group_MemSac_ps(:,3);
MemSac_indicator_txt = 'MemSac\_DDI ([3])';

% MemSac_indicator = group_MemSac_actual_DI; % Actual target locations DI
% MemSac_indicator (~isnan(group_MemSac_actual_DI_2pt)) = group_MemSac_actual_DI_2pt(~isnan(group_MemSac_actual_DI_2pt)); % 2pt Mem is prior
% MemSac_indicator = abs(MemSac_indicator);
% MemSac_indicator_txt = 'MemSac\_actual\_DI ([3])';

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


select_all = []; select_sus = []; select_bottom_line = []; select_tcells = []; select_no_tcells = []; select_bottom_line_all_monkey = [];
select_psy_good = []; select_psy_bad = []; find_bottom_line = [];

Choice_pref_all = []; Choice_pref_p_value_all = []; Modality_pref_all = []; Modality_pref_p_value_all = [];
select_cpref_mpref = []; 

t_criterion_txt = [];

group_position = [];  % All cells
Position_all = []; % Could be one monkey
cell_position();
cell_selection();

%% Cell Selection and Cell Counter
    function cell_selection(t_cell_selection_num)  % Cell Selection and Cell Counter

        if nargin < 1
            t_cell_selection_num = 6; % Default
        end
            
        % --------  @ HH20150413 -------- 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Limitations on petition number & Target first
        select_all_all_monkey = ([group_result.repN]' >= 8) & ([group_result.length_unique_stim_type]' == 3 ) & (xls_num{1}(:,header.HD_TargFirst)~=0);
        
        % + SUs
        select_sus_all_monkey = select_all_all_monkey & (xls_num{1}(:,header.Units_RealSU) > 0) & (xls_num{1}(:,header.Chan1) < 20);
        
        % Bottom line for most figures
        select_bottom_line_all_monkey = select_sus_all_monkey; 
        
        % + T(ypical) CellsZ
        t_cell_selection_criteria = ...
        {  %                   Logic                                   Notes
           % Bottom-line
            'Bottom-line (all)', ones(size(select_sus_all_monkey)); % Just all bottom-line cells
           % Mem-sac based
            'Manually assigned mem-sac (original)',    xls_num{1}(:,header.HD_MemSac) >= 0.8;
            'Memory p value',   group_MemSac_ps(:,3) < 0.05 | group_TwoPtMemSac_ps(:,3) < 0.05;
            'Memory p value & actual DDI', (group_MemSac_ps(:,3)<0.05 & abs(group_MemSac_actual_DI) >= 0.3) | (group_TwoPtMemSac_ps(:,3)<0.05 & abs(group_TwoPtMemSac_DDI(:,3)) >=0.3); 
            'Global DDI of mem-sac',   mean(group_MemSac_DDI(:,[3 4]),2) >= 0.5;
           % HD based
            'ChoicePreference p value (any)', any(group_ChoicePreference_pvalue(:,:,3) < 0.05,1)';  % 3, stim-on to stim-off
            'ChoicePreference p value (all)', all(group_ChoicePreference_pvalue(:,:,3) < 0.05,1)';  
            'ChoicePreference p value (vest)', (group_ChoicePreference_pvalue(1,:,3) < 0.05)';  % 3, stim-on to stim-off
            'ChoicePreference p value (vis)', (group_ChoicePreference_pvalue(2,:,3) < 0.05)';  % 3, stim-on to stim-off
            'ChoicePreference p value (comb)', (group_ChoicePreference_pvalue(3,:,3) < 0.05)';  % 3, stim-on to stim-off
           % For debugging
            'HD_dt patched', [group_result(:).HD_dt_patched]';
            'non HD_dt patched', ~[group_result(:).HD_dt_patched]';
            'Broad-spike',[group_result(:).Waveform_broad]';
            'Narrow-spike',~[group_result(:).Waveform_broad]';
            'Broad-spike & ChoicePreference p value (any)',[group_result(:).Waveform_broad]' & any(group_ChoicePreference_pvalue(:,:,3) < 0.05,1)';
            'Narrow-spike & ChoicePreference p value (any)',~[group_result(:).Waveform_broad]' & any(group_ChoicePreference_pvalue(:,:,3) < 0.05,1)';
            'Broad-spike & ChoicePreference p value (all)',[group_result(:).Waveform_broad]' & all(group_ChoicePreference_pvalue(:,:,3) < 0.05,1)';
            'Narrow-spike & ChoicePreference p value (all)',~[group_result(:).Waveform_broad]' & all(group_ChoicePreference_pvalue(:,:,3) < 0.05,1)';
            'Shallow layer',group_position(:,5) < 1100;
            'Deep layer',group_position(:,5) > 1100;
            'Shallow & Broad', group_position(:,5) < 1100 & [group_result(:).Waveform_broad]';
            'Shallow & Narrow', group_position(:,5) < 1100 & ~[group_result(:).Waveform_broad]';
            'Deep & Broad',group_position(:,5) >= 1100 & [group_result(:).Waveform_broad]';
            'Deep & Narrow',group_position(:,5) >= 1100 & ~[group_result(:).Waveform_broad]';
            'Combined > Max (Vis, Vest) and Comb p < 0.05',(group_ChoicePreference_pvalue(3,:,3) < 0.05)' & (abs(group_ChoicePreference(3,:,3)) > abs(group_ChoicePreference(1,:,3)))' & (abs(group_ChoicePreference(3,:,3)) > abs(group_ChoicePreference(2,:,3)))';
        };
        
        select_tcells_all_monkey = select_sus_all_monkey & t_cell_selection_criteria{t_cell_selection_num,2};   
        select_no_tcells_all_monkey = select_sus_all_monkey & ~ t_cell_selection_criteria{t_cell_selection_num,2};  
        t_criterion_txt = t_cell_selection_criteria{t_cell_selection_num,1};
        
        % ---------
        
        
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
        set(h_all,'string',sprintf('%7d%7d%7d\n',cell_nums'),'fontsize',13);
        h_t_criterion = findall(gcbf,'tag','t_criterion');
        set(h_t_criterion,'string',{t_cell_selection_criteria{:,1}});
        set(h_t_criterion,'value',t_cell_selection_num);
        
        %         h_su = findall(gcbf,'tag','num_su');
        %         set(h_su,'string',num2str(sum(select_sus)));
        %         h_tcell = findall(gcbf,'tag','num_t_cells');
        %         set(h_tcell,'string',num2str(sum(select_tcells)));
        
        % -------- Update/refresh some related datasets that are influenced by cell_selection ---------
        % For Choice and modality preference
        select_cpref_mpref = select_bottom_line;
        Choice_pref_all = reshape([group_result(select_cpref_mpref).ChoicePreference]',3,[],3);  % Stim, Cell No, Pre/Post. Updated time win HH20160918
        Choice_pref_p_value_all = reshape([group_result(select_cpref_mpref).ChoicePreference_pvalue]',3,[],3); % Updated time win. HH20160918
        Modality_pref_all = reshape([group_result(select_cpref_mpref).ModalityPreference]',3,[],3);
        Modality_pref_p_value_all = reshape([group_result(select_cpref_mpref).ModalityPreference_pvalue]',3,[],3);
        
        Position_all = group_position(select_bottom_line,:);
        
        select_for_SVM = select_bottom_line;
        select_for_PCA_B = select_bottom_line;
        
        PCA_A = [];  % Reset PCA_A
        PCA_B_projPC = []; % Reset PCA_B
        thres_choice = []; % Reset SVM training
        weights_PCA_B_PC = []; weights_svm_choice_mean = []; weights_TDR_PCA_SVM_mean = []; % Reset TDR
        
    end

%% Load DrawMapping data to get the cell location. HH20161024
    function cell_position()
        
        GM = -1;    % Gray matter without visual modulation
        LIP = -3;
        VIP = -4;
        MST = -5;
        MT = -6;
        
        drawmapping_data = { % [Monkey,hemi]    % Grid No. of (AP0,Middle)   % Data
            [5,1],[7,0];  % 'Polo_L'
            [5,2],[9,-5]; % 'Polo_R'
            [10,1],[10,-7]; % 'Messi_L'
            [10,2],[10,-4]}; % 'Messi_R'
        
        %  Messi_right
        drawmapping_data{4,3} = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {8,[11,10],[1.9 0.0],[GM 8 17; GM 58 72; VIP 95 118]}
            {40,[13,18],[2.2 0.0-0.3],[GM 31 60; MST 75 103;  MT 120 147]}
            {40,[13,8],[2.05 0.0],[GM 37 59; VIP 67 86]}
            {41,[13,13],[2.05 0.2],[GM 7 18; LIP 31 43], 43}
            {42,[13,12],[1.9 0],[GM 13 56; LIP 74 84],84}
            {43,[13,14],[2.05 0],[GM 8 28; LIP 38 70; MST 92 107; MST 117 135]}
            {44,[13,15],[2.2 0],[GM 7 13; LIP 23 45; MST 69 90; MT 100 130]}
            {45,[13,11],[1.9 0],[GM 5 17; GM 22 63; LIP 73 96], 150}
            {45,[13,2],[1.9 0],[GM 9 49]}
            {46,[11,14],[2.1 0],[GM 3 37; LIP 51 63],63}
            {47,[9,15],[2.2 0+0.1],[GM 10 20; LIP 32 63; MST 114 130]}
            {48,[7,17],[2.2 0],[GM 12 20; LIP 32 63; GM 82 112; GM 132 146]}
            {49,[9,15],[2.2 0],[GM 4 23; LIP 35 40],40}
            {50,[10,14],[2.2 0],[GM 0 22; LIP 46 60],60}
            {53,[13,12],[1.9 0],[GM 13 58; LIP 73 85],85}
            {54,[12,13],[2.05 0],[GM 4 34; LIP 52 67],67}
            {55,[12,12],[2.05 0],[GM 14 51; LIP 69 75],75}
            {56,[12,11],[2.1 0],[GM 10 57; LIP 65 78],78}
            {57,[13,14],[2.05 0.1],[GM 0 22; LIP 33 60],60}
            {58,[13,15],[2.2 0],[LIP 20 39; MST 63 70],70}
            {72,[14,12],[2.05 0],[GM 2 30; LIP 51 75]}
            {73,[14,13],[1.9 0],[GM 13 34; LIP 50 67],67}
            {74,[14,11],[2.05 0],[GM 5 47; LIP 62 70],70}
            {75,[14,14],[1.9+0.15 0],[LIP 33 66; MST 92 95],95}
            {76,[15,10],[1.9 0],[GM 27 62; LIP 75 93],93}
            {77,[15,11],[2.05-0.2 0],[GM 27 59; LIP 72 83],83}
            {78,[15,12],[2.05+0.1 0],[GM 0 20; LIP 34 57],57}
            {79,[16,10],[2.05 0],[GM 2 32; LIP 43 75]}
            }';
        
        %  Messi_left
        drawmapping_data{3,3} = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {1,[15,15],[2.2 0.0],[GM 4 37; GM 56 74; MST 88 117]}
            {2,[15,13],[2.1 0.0],[GM 5 18; GM 31 52; MST 74 96; MST 106 124]}
            {3,[15,17],[2.4 0.0],[GM 26 56; MT 69 122]}
            {3,[15,11],[2.1 0.0 + 0.1],[GM 0 20; GM 29 47; MST 75 96;]}
            {3,[15,9],[2.1 0.0 + 0.1],[GM 6 36; GM 43 67; GM 116 149]}
            {4,[15,19],[2.25 0.0],[GM 7 79],79}
            {4,[13,20],[2.25 0.0],[GM 9 39; MST 89 135; MT 142 150]}
            {4,[13,13],[2.1 0.0],[GM 7 37; GM 44 72; MST 92 122; MT 134 150]}
            {5,[13,17],[2.25 0.0],[GM 5 50; MST 66 90; MT 109 140]}
            {6,[15,7],[2.05 0.0],[GM 33 59; GM 71 90; LIP 108 140]}
            {6,[11,5],[1.8 0.0],[GM 1 26; GM 39 58; GM 66 84; VIP 104 121; VIP 129 143]}
            {7,[11,9],[1.9 0.0],[GM 2 30; GM 65 97; LIP 105 125],125}
            {8,[11,10],[1.9 0.0],[GM 4 33; GM 55 85; LIP 97 123]}
            {9,[9,12],[2.0 0-0.2],[GM 18 39; GM 69 98; LIP 108 128]}
            {10,[9,11],[2.0 0],[GM 7 35;GM 67 95;LIP 106 121],121}
            {11,[13,9],[1.9 0+0.1],[GM 44 73; LIP 82 104]}
            {12,[13,8],[1.9 0],[GM 58 86; LIP 93 102],102}
            {13,[13,7],[1.8 0],[GM 0 13; GM 29 50; GM 70 96; LIP 104 116],116}
            {14,[13,6],[1.7 0],[GM 10 24; GM 37 49; GM 55 70; GM 89 111; LIP 117 129],129}
            {15,[12,8],[1.8 0.0 + 0.1],[GM 6 19; GM 64 93; LIP 102 116],116}
            {16,[12,7],[1.8 0.0 + 0.1],[GM 2 9; GM 31 46; GM 74 97; LIP 105 114; VIP 114 128]}
            {17,[12,9],[1.9 0.0 + 0.1],[GM 0 7; GM 49 79; LIP 86 90],90}
            {18,[12,10],[2.0 0 + 0.1],[GM 3 12; GM 37 64; LIP 73 80],80}
            {19,[12,11],[2.0 0 + 0.1],[GM 3 10; GM 26 55; LIP 64 88],96}
            {20,[11,11],[2.0 0 - 0.05],[GM 4 19; GM 44 79; LIP 88 105],105}
            {21,[14,7],[1.8 0 + 0.1],[GM 20 34; GM 56 81; LIP 89 101],101}
            {22,[14,8],[1.9 0 + 0.1],[GM 42 69; LIP 80 101],120}
            {23,[14,9],[1.9 0 + 0.1],[GM 36 64; LIP 72 79],79}
            {24,[10,12],[2.0 0-0.1],[GM 9 26; GM 48 82;LIP 94 115]}
            {24,[12,10],[2.0 0+0.1],[GM 3 12; GM 36 62; LIP 72 77],77}
            {25,[12,10],[2.0 0+0.1],[GM 36 62; LIP 72 86],86}
            {26,[12,11],[2.0 0+0.1],[GM 30 58; LIP 67 75],75}
            {27,[18,5],[1.7 0.4+0.05],[GM 11 36; LIP 47 65; GM 77 100; GM 109 130]}
            {28,[18,4],[1.7 0.4],[GM 18 42; LIP 53 61; GM 71 91; GM 103 125]}
            {29,[12,12],[2.0 0],[GM 9 18; GM 29 58; LIP 66 74],74}
            {30,[12,13],[2.1 0],[GM 19 39; LIP 50 60],60}
            {31,[11,13],[2.1 0],[GM 0 11; GM 22 55; LIP 65 80],80}
            {32,[11,12],[2.1 0],[GM 0 7; GM 26 57; LIP 69 82],82}
            {33,[11,11],[2.0 0.1],[GM 0 10; GM 34 61; LIP 76 97]}
            {34,[13,11],[2.0 0.1],[GM 26 52; LIP 62 80]}
            {35,[13,10],[2.0 0],[GM 35 65; LIP 75 96],96}
            {36,[14,10],[2.0 0+0.1],[GM 20 46; LIP 55 80],90}
            {37,[13,12],[2.0 0+0.1],[GM 3 8; GM 20 47; LIP 56 82]}
            {38,[14,11],[2.0 0+0.1],[GM 0 38; LIP 48 65]}
            {39,[15,10],[2.0 0],[GM 22 48; LIP 59 81; GM 120 150]}
            {59,[13,11],[2.05 0.0 + 0.1],[GM 13 39; LIP 56 74],74}
            {60,[14,11],[2.0 0],[GM 3 52;LIP 57 60],60}
            {61,[14,10],[2.0 0.0 + 0.1],[GM 21 45; LIP 57 75],75}
            {62,[13,12],[2.0 0.0 + 0.2],[GM 2 36;LIP 45 65],65}
            {64,[11,13],[2.1 0],[GM 15 33; LIP 65 78],78}
            {65,[14,9],[1.9 0.0 + 0.2],[GM 22 51; LIP 63 84],84}
            {66,[14,7],[1.8 0],[GM 29 41; GM 65 89; LIP 99 120],120}
            {66,[17,7],[1.8 0.0 + 0.2],[GM 2 6; GM 35 59; LIP 70 90; LIP 113 150],150}
            {67,[15,8],[1.9 0 + 0.1],[GM 34 62; LIP 74 84],84}
            {68,[15,9],[1.9 0],[GM 47 58; LIP 75 83],83}
            {69,[12,9],[1.9 0 + 0.15],[GM 45 73; LIP 86 97],97}
            {70,[16,7],[1.9 0],[GM 41 69; LIP 77 101; VIP 138 150],150}
            {71,[16,6],[1.9 0 + 0.1],[GM 7 19; GM 39 67; LIP 76 80],80}
            {71,[13,7],[1.9 0],[GM 20 38; GM 62 89; LIP 96 110],110}
            }';
        
        %  Polo_left
        
        drawmapping_data{1,3} = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {62,[6,18],[2.0 0.0],[GM 30 43; LIP 64 72; VIP 72 84; MST 136 150]}
            {63,[6,19],[2.0 0],[GM 4 19; GM 22 39; LIP 56 63],63}
            {64,[6,15],[2.0 0],[GM 15 27; GM 48 69; VIP 81 84],84}
            {65,[6,16],[1.9 0],[GM 52 75; LIP 82 97],97}
            {66,[6,20],[2.1 0],[GM 2 19;LIP 37 60],60}
            {67,[6,21],[2.1 + 0.1 0],[LIP 20 43],43}
            {68,[6,17],[2.0 + 0.1 0],[GM 0 13; GM 25 48; LIP 55 64; LIP 64 75],75}
            {69,[5,19],[2.0  0],[GM 0 41; LIP 53 59],59}
            {70,[5,18],[2.0 0],[GM 0 3; GM 28 53; LIP 62 72],72}
            {71,[5,17],[2.0 0.1],[GM 0 11; GM 27 48; LIP 57 70],70}
            {72,[6,19],[2.1 0],[LIP 46 65]}
            {73,[6,18],[2.1 0],[GM 0 12; GM 26 38; LIP 47 59],59}
            {74,[7,18],[2.0 0],[GM 6 49; LIP 59 73],73}
            {75,[7,17],[2.0 0],[GM 7 54; LIP 62 75],75}
            {76,[7,16],[2.0 0],[GM 2 27; GM 40 61; LIP 73 80; VIP 80 86],86}
            {77,[7,19],[2.1 0],[GM 0 11; LIP 37 53],53}
            {78,[8,17],[2.0 0],[GM 9 50; LIP 66 78],78}
            {79,[8,18],[2.0 0],[GM 0 44; LIP 56 57],57}
            {80,[8,19],[2.1 0],[GM 0 24; LIP 48 60]}
            {81,[7,17],[2.1 0],[GM 0 14; GM 26 46; LIP 59 66]}
            {82,[6,18],[2.1 0],[GM 0 13; GM 27 41; LIP 56 71]}
            {83,[6,12],[2.0 0],[VIP 63 71],71}
            {83,[6,10],[2.0 - 0.1 0],[GM 17 73; VIP 102 123]}
            {84,[6,14],[2.0 0],[GM 5 20; GM 52 81; VIP 92 113]}
            {85,[7,18],[2.1-0.1 0],[GM 8 48; LIP 58 73],73}
            {86,[8,18],[2.1 0],[GM 0 35; LIP 46 65; MST 99 100],100}
            {87,[4,18],[2.0 0],[GM 0 10; GM 31 51; LIP 62 87]}
            {88,[4,17],[2.0 0],[GM 35 58; LIP 71 94 ]}
            {89,[4,19],[2.1-0.1 0],[GM 20 45; LIP 61 87]}
            {90,[3,17],[2.0 0],[GM 35 70; LIP 83 97; VIP 97 110]}
            {91,[3,18],[2.0 0],[GM 0 10; GM 37 60; LIP 73 92; VIP 92 98]}
            {91,[3,19],[2.0 0],[GM 30 45; LIP 63 80]}
            {92,[10,17],[2.0 0],[GM 3 45; LIP 52 69],69}
            {93,[12,15],[2.1 0],[GM 3 27; LIP 37 61]}
            {94,[12,13],[2.0 0],[GM 21 49; LIP 56 72],72}
            {95,[14,12],[1.9 -0.1],[GM 10 25; GM 35 56; LIP 67 89]}
            {96,[16,10],[1.9 -0.1],[GM 9 31; GM 42 70; LIP 100 120; LIP 133 148]}
            {97,[10,14],[2.0 0],[GM 3 10; GM 25 61; LIP 71 81],81}
            {98,[8,15],[2.0 -0.05],[GM 3 22; GM 42 65; LIP 78 89],89}
            {99,[9,15],[2.0 0],[GM 5 62; LIP 73 86],86}
            {100,[11,14],[2.0 0],[GM 3 56; LIP 66 81],81}
            {101,[13,13],[2.0 0],[GM 11 42; LIP 52 54],54}
            {102,[8,11],[2.0 0],[GM 2 17; GM 35 42; VIP 65 89; VIP 95 110]}
            {103,[8,10],[1.9 0],[GM 11 26; GM 36 52; VIP 73 92; VIP 101 119]}
            {104,[8,9],[1.9 0],[GM 14 30; GM 40 59; VIP 82 96 ; VIP 107 115]}
            {105,[8,12],[1.9 0],[VIP 67 91; VIP 103 126]}
            {106,[9,11],[1.8 0-0.1],[GM 30 58; GM 80 102; VIP 117 132],132}
            {107,[9,10],[1.8 0],[GM 16 52; GM 74 91; VIP 107 125]}
            {108,[10,10],[1.8 0],[GM 15 58; GM 68 83; VIP 98 125]}
            {109,[10,9],[1.8 0],[GM 4 25; GM 34 51; GM 72 85; VIP 101 117],117}
            {111,[10,15],[2.1 0+0.1],[GM 0 35; LIP 44 59],59}
            }';
        
        %  Polo_right
        drawmapping_data{2,3} = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {1,[5,15],[2.0 0.0],[GM 0 20; LIP 44 69; MST 147 150]}
            {2,[5,17],[2.0 0.0],[GM 0 17; LIP 32 64],64}
            {3,[5,13],[1.9 0.0],[GM 21 37; LIP 71 96]}
            {3,[5,11],[1.9 0.0],[GM 31 58; LIP 74 97; VIP 97 103]}
            {4,[6,7],[1.7 0.0],[GM 0 21; VIP 78 100; VIP 111 129],129}
            {5,[5,16],[1.7 0.0],[GM 21 51; LIP 76 82],82}
            {6,[5,14],[1.8 0.0],[GM 0 15; GM 30 43; LIP 79 98]}
            {7,[5,12],[1.8 0.0],[GM 38 50; LIP 79 107],107}
            {8,[9,10],[1.8 0.0],[GM 43 73; LIP 91 111],111}
            {9,[9,9],[1.9 0.0],[GM 0 12; GM 35 59; LIP 75 94],94}
            {10,[9,11],[1.9 0.0],[GM 28 56; LIP  75 90],90}
            {11,[9,12],[1.9 0.0],[GM 24 49; LIP 68 69], 69}
            {12,[9,13],[1.9 0.0],[GM 17 35; LIP 53 62; MST 101 105 ], 105}
            {13,[7,11],[1.9 0.0],[GM 0 9; GM 25 47; LIP 61 86; MST 127 149]}
            {14,[7,12],[1.9 0.0],[GM 0 11; GM 26 47; LIP 58 76],76}
            {15,[7,13],[1.9 0.0],[GM 0 11; GM 20 41; LIP 50 68],68}
            {16,[2,10],[1.9 0.0],[VIP 65 86],86}
            {17,[2,11],[1.9 0.0],[GM 49 59 ; VIP 59 77; VIP 90 103]}
            {18,[2,12],[1.9 0.0],[GM 47 75; VIP 88 102]}
            {19,[8,11],[1.9 0.0],[GM 25 53; LIP 61 84; MST 118 138]}
            {20,[3,10],[1.9 0.0],[VIP 57 73],73}
            {21,[3,9],[1.9 0.0],[VIP 68 92],92}
            {22,[3,11],[1.9 0.0], [GM 47 73; LIP 92 99; VIP 99 110]}
            {23,[2,9],[1.9 0.0],[VIP 77 105]}
            {25,[4,9],[1.9 0.0],[GM 68 80; VIP 80 85],85}
            {26,[4,10],[1.9 0],[GM 49 81; VIP 95 107],107}
            {27,[6,11],[1.9 0],[GM 36 66; LIP 79 100],100}
            {28,[6,12],[1.9 0],[GM 29 55; LIP 71 90; MST 142 150]}
            {29,[6,13],[2.1 0],[GM 13 36; LIP 44 60],60}
            {30,[8,10],[1.9 0],[GM 9 25; GM 34 54; LIP 63 73],73}
            {31,[8,9],[1.9 0],[GM 11 22; GM 39 63; LIP 77 84],84}
            {32,[7,10],[1.9 0],[GM 2 17; GM 39 64; LIP 74   98],98}
            {33,[8,12],[2.0 0],[GM 4 34; LIP 45 65; MST 108 110],110}
            {34,[4,12],[2.0 0],[GM 24 43; LIP 67 85],85}
            {35,[4,13],[2.0 0],[GM 19 38; LIP 60 82],82}
            {36,[4,14],[2.0 0],[GM 0 38; LIP 55 78],78}
            {37,[3,13],[2.0 0],[GM 22 43; LIP 72 86]}
            {39,[8,10],[2.0 0],[GM 10 47; LIP 57 81],81}
            {40,[7,11],[2.0 0],[GM 24 45; LIP 58 70; VIP 70 80; MST 125 129],129}
            {38,[3,14],[2.0 0],[GM 18 46; LIP 66 84]}
            {41,[7,12],[2.0 0],[GM 14 37; LIP 48 66], 66}
            {42,[6,14],[2.0 0],[GM 5 30; LIP 45 58],58}
            {43,[6,15],[2.0 0],[GM 6 27; LIP 43 70; MST 144 150]}
            {44,[6,16],[2.0 0],[LIP 30 67; MST 106 135]}
            {45,[3,15],[2.0 0],[GM 11 35; LIP 61 63],63}
            {46,[3,16],[2.0 0],[GM 13 35; LIP 52 56],56}
            {47,[3,18],[2.0 0],[GM 10 38; LIP 53 57],57}
            {48,[3,17],[2.1 0],[GM 0 25; LIP 42 56],56}
            {49,[4,15],[2.0 0],[GM 2 35; LIP 55 71],71}
            {50,[2,14],[2.0 0],[GM 0 2; GM 30 55; LIP 69 73; VIP 73 86];}
            {51,[2,15],[2 0],[GM 27 53; LIP 68 73], 73}
            {52,[2,16],[2.1 0],[GM 0 43; LIP 54 65],65}
            {53,[10,10],[1.9 0],[GM 0 46; VIP 60 68],68}
            {54,[10,12],[2.0 0],[GM 6 15; LIP 33 59; MST 85 90],90}
            {55,[10,13],[2.0 0],[LIP 30 47; MST 92 93],93}
            {56,[8,13],[2.0 0],[GM 0 22; LIP 41 60],60}
            {57,[7,14],[2.0 0],[GM 5 28; LIP 36 58; MST 120 144]}
            {58,[6,13],[2.0 0],[GM 19 33; LIP 46 48],48}
            {59,[6,12],[1.9 0],[GM 40 50; LIP 71 86; MST 148 150]}
            {60,[5,13],[2.0 0],[GM 12 27; LIP 57 62],62}
            {61,[5,14],[2.0 0],[LIP 58 67],67}
            }';
        
        % === Get cell positions ===
        for i = 1:length(group_result)
            % Decode cell position
            id = sscanf(group_result(i).cellID{1}{1},'%gm%gs%gh%gx%gy%gd%g');
            this_monkey_hemi = id([2 4])';
            this_session = id(3);
            this_pos_raw = id(5:7)';
            
            % Find entry in drawmapping_data
            found = 0;
            for mmhh = 1:size(drawmapping_data,1)
                if all(drawmapping_data{mmhh,1} == this_monkey_hemi)
                    xyoffsets = drawmapping_data{mmhh,2};
                    this_mapping_data = drawmapping_data{mmhh,3};
                    
                    for ee = 1:length(this_mapping_data)
                        if all(this_mapping_data{ee}{1} == this_session) && all(this_mapping_data{ee}{2} == this_pos_raw(1:2))
                            found = found + 1;
                            
                            % Anterior-posterior
                            AP = - (this_pos_raw(1) - xyoffsets(1)) * 0.8; % in mm, P - --> AP0 --> + A
                            % Ventral-dorsal
                            VD = (this_pos_raw(2) - xyoffsets(2)) * 0.8 - AP * tan(30/180*pi); % LIPs are inclined by about 30 degree
                            
                            % Depth to the surface of cortical sheet
                            which_GM_is_LIP = find((this_mapping_data{ee}{4}(:,1) == LIP));
                            which_LIP_this_cell_in =  find((this_mapping_data{ee}{4}(which_GM_is_LIP,2) <= ceil(this_pos_raw(3)/100)) & ...
                                (this_mapping_data{ee}{4}(which_GM_is_LIP,3) >= fix(this_pos_raw(3)/100)));
                            
                            if which_LIP_this_cell_in == 1 % This cell is in the first (upper) layer of LIP
                                this_surface = this_mapping_data{ee}{4}(which_GM_is_LIP(which_LIP_this_cell_in),2); % Surface is the start of this GM
                                depth = this_pos_raw(3) - 100 * this_surface; % in um
                            elseif which_LIP_this_cell_in == 2 % This cell is in the second (lower) layer of LIP, the depth should be inversed.
                                % this_surface = this_mapping_data{ee}{4}(which_GM_is_LIP(which_LIP_this_cell_in),3); % Surface is the end of this GM
                                % depth = - (this_pos_raw(3) - 100 * this_surface); % in um
                                depth = nan; % Because I'm not sure whether the end of penetration has reached the (bottom) surface of the second LIP
                            else
                                disp('Check LIP layers!!!'); beep; keyboard
                            end
                            
                            group_result(i).position = [this_monkey_hemi AP VD depth which_LIP_this_cell_in]; % [monkey hemi AP Ventral-Dorsal depth which_LIP_this_cell_in]
                            
                        end
                    end
                end
            end
            
            if found~=1 % Should be exactly 1
                disp('Cell position not found!!');
                beep; keyboard;
            end
        end
        
        group_position = reshape([group_result.position],6,[])';

    end

%% Final Preparation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===========  Common parameters   ============
% Overall time marker (Stim on, Stim off, Sac on)
tt = reshape(cell2mat([group_result.time_marker]),6,[])';

for j = 1:2
    time_markers{j} = [mean(tt(:,(j-1)*3+1:(j-1)*3+3)); std(tt(:,(j-1)*3+1:(j-1)*3+3))];
end

p_critical = 0.05;

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
% PCA_B_time_range = min(rate_ts{j_PCA_B})+100 <= rate_ts{j_PCA_B} & rate_ts{j_PCA_B} <= time_markers{j_PCA_B}(1,2);  % Before stim 
PCA_B_times = rate_ts{j_PCA_B}(PCA_B_time_range);
denoised_dim = 9;
PCA_B = []; weights_PCA_B_PC = []; PCA_B_projPC = []; PCA_B_explained = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ================ Miscellaneous ===================

% set(0,'defaultAxesColorOrder',[0 0 1; 1 0 0; 0 0.8 0.4;]);
% colors = [0 0 1; 1 0 0; 0 0.8 0.4];

set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);
colors = [41 89 204; 248 28 83; 14 153 46]/255;


modality_diff_colors = colors;    modality_diff_colors(3,:) = [0 0 0];

transparent = get(findall(gcbf,'tag','transparent'),'value');
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
    'Mean Rate Metrics', {
    'Correct only, all choices (conventional)', @f1p1;
    '   All single cell traces',@f1p1p6;
    '      Example cells (SuppFig.2)', @f1p1p4;
    '      Divergence time', @f1p1p6p1;
    '   All single cell fittings',@f1p1p7;
    '   Different weighting methods',@f1p1p5;
    'Different headings' , @f1p2;
    '   Comb/max ratio distribution through time', @f1p2p5;
    '   Tuning curve analysis',@f1p2p7;
    '   Partial correlation (grand)',@f1p2p8;
    'Correct / Wrong Trials', @f1p3;
    };
    
    'Variance Metrics',{
    'VarCE',@f10p1;
    
    };
    
    'ROC Metrics',{
    'CP, CDiv, and MDiv',@f2p1;
    'Multisensory Enhancement of CDiv',@f2p2;
    'Easy and Difficult',@f2p3;
    };
    
    'Correlations', {
    'Mem-sac vs. pre/post CDiv/CP, CDiv vs. CP', @f3p1;
    'Mem-sac vs. abs(CPref)', @f3p1p2;
    'Choice Preference vs. Modality Preference', @f3p2;
    'Choice Preference between modalities', @f3p2p2;
    'Choice Preference: pre and post', @f3p2p3;
    'Cell position and type',@f3p4;
    '   Regression',@f3p4p1;
    };
    
    'PCA_choice & modality (Eigen-feature)',{
    'Hot-gram', @f4p1;
    'Cluster and Trajectory', @f4p2;
    };
    
    'PCA_choice & modality (Eigen-neuron)',{
    'Weights and correlations', @f5p1;
    '1-D Trajectory',@f5p2;
    '3-D Trajectory',@f5p3;
    };

    'PCA_heading & choice (Eigen-neuron)',{
    '2-D Trajectory', @f51p1;
    };
    
    'demixed PCA and significant% (2nd reviewer)',{
    '% Significant cells from t-test',@f55p4;
    '% Significant cells from ANCOVA',@f55p5;    
    't+st+dt+sdt, for each stim_type separately', @f55p1;
    '==> t+mt+dt+mdt, stim_type together, heading marginalized', @f55p2;
    't+st+dt+mt, heading+choice+modality, no interaction', @f55p3;
    };
    

    'Linear SVM decoder (choice and modality)',{
    'Training SVM', @f6p0;
    'Weights', @f6p1;
    'Performance (overall)', @f6p2;
    'SVM weighted sum',@f6p3;
    }
    
    'Fisher information of heading',{
    'Simple Fisher like Gu 2010: Sum(slope/mean)', @f6p5p1
    '   + Correlation between modalities',@f6p5p1p1
    '   + Correlation with CPref',@f6p5p1p2
    'Dora''s partial corr and linear regression k', @f6p5p2
    'Training SVM decoders', @f6p5p9;
    } 
        
    'Linear regression',{
    'Comb = w1 Vest + w2 Vis (Fig.5 in Gu 2008)',@f7p1;
    'Fit model traces with real data (Alex Shanghai)',@cell_selection;
    '     1). All heading mean',@f7p2;
    '     2). Each |heading|',@f7p3;
    };
    
    'Targeted Dimensionality Reduction',{
    'PCA + SVM: weights', @f8p1;
    'PCA + SVM: all correct + all angles',@f8p2;
    }
    
    'Others',{
    'Behavior (recording sessions only)', @f9p3;
    '  Psychophysics vs. Neural activity', @f9p4;
    'Cell Counter',@f9p1;
    'Target first vs. Target last',@f9p2;
    '-----------------------------------', @cell_selection;
    'Export Associated Memsac Files',@f9p9;
    };
    
    'NoShow',{@cell_selection};
    
    };

%% ====================================== Function Definitions =============================================%

    function f1p1(debug)      % Rate 1. Correct only, all choices
        if debug  ; dbstack;   keyboard;      end
        
        %% ------- Averaged norm PSTH --------
        set(figure(999),'name','Average PSTH (Correct only, all choices)','pos',[27 63 1449 892]); clf
        h_subplot = tight_subplot(2,3,[0.11 0.05],[0.05 0.1],0.07);
        
        methods_of_select = {
            select_bottom_line, 'All cells';
            select_tcells, 'Typical cells';
            select_no_tcells, 'Non-typical cells'};
        
        %         for j = 1:2
        j = 1;
        for ms = 1:3
            
%             SeriesComparison({shiftdim(ys_this{1},-1) shiftdim(ys_this{2},-1)},...
%                 {rate_ts{1} rate_ts{2} time_markers},...
%                 'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'figN',1463);
%             
            h = SeriesComparison({PSTH_all_Norm{1}(methods_of_select{ms,1},:,:), PSTH_all_Norm{2}(methods_of_select{ms,1},:,:)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),colors(3,:),colors(3,:)},'LineStyles',{'-','--','-','--','-','--'},...
                'ErrorBar',6,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot((ms-1)*2+1),...
                'CompareIndex',[1,3,5;2,4,6],...
                'CompareColor',[mat2cell(colors,ones(3,1))],...
                'Transparent',transparent);
            
%             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end

            xlabel('Time (ms)');
            legend off;
            title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
            axis tight;
            
%             for tt = 1:3
%                 plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
%             end

            axis tight;
            xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  ylim([0.1 0.7]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
            
            % --- Difference (Pref - Null) ---
            PSTH_all_Norm_PrefminusNull{1} = PSTH_all_Norm{1}(methods_of_select{ms,1},:,1:2:end)...
                - PSTH_all_Norm{1}(methods_of_select{ms,1},:,2:2:end);
            PSTH_all_Norm_PrefminusNull{2} = PSTH_all_Norm{2}(methods_of_select{ms,1},:,1:2:end)...
                - PSTH_all_Norm{2}(methods_of_select{ms,1},:,2:2:end);
            
            SeriesComparison({PSTH_all_Norm_PrefminusNull{1}, PSTH_all_Norm_PrefminusNull{2}},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',mat2cell(colors,ones(3,1)),'LineStyles',{'-'},...
                'ErrorBar',6,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot((ms-1)*2+2),...
                'CompareIndex',[1:3,1,2;1:3,3,3],...
                'CompareColor',[mat2cell(colors,ones(3,1));colors(1,:);colors(2,:);colors(3,:)],...
                'Transparent',transparent);
            
            %             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            xlabel('Time (ms)');
            legend off;
            
            axis tight;
            xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);   ylim([-0.1 0.4]);
            
            % Gaussian vel
%             plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/6,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;

            
        end
        %         end
        SetFigure(15);
        
        %% ------- Averaged raw PSTH --------
        set(figure(1999),'name','Average PSTH (Correct only, all choices)','pos',[27 63 1449 892]); clf
        h_subplot = tight_subplot(2,3,[0.11 0.05],[0.05 0.1],0.07);
        
        methods_of_select = {
            select_bottom_line, 'All cells';
            select_tcells, 'Typical cells';
            select_no_tcells, 'Non-typical cells'};
        
        %         for j = 1:2
        j = 1;
        for ms = 1:3
            
%             SeriesComparison({shiftdim(ys_this{1},-1) shiftdim(ys_this{2},-1)},...
%                 {rate_ts{1} rate_ts{2} time_markers},...
%                 'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'figN',1463);
%             
            h = SeriesComparison({PSTH_all_raw{1}(methods_of_select{ms,1},:,:), PSTH_all_raw{2}(methods_of_select{ms,1},:,:)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),colors(3,:),colors(3,:)},'LineStyles',{'-','--','-','--','-','--'},...
                'ErrorBar',6,'Xlabel',[],'Ylabel','Raw firing','axes',h_subplot((ms-1)*2+1),...
                'CompareIndex',[1,3,5;2,4,6],...
                'CompareColor',[mat2cell(colors,ones(3,1))],...
                'Transparent',transparent);
            
%             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end

            xlabel('Time (ms)');
            legend off;
            title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
            axis tight;
            
%             for tt = 1:3
%                 plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
%             end

            axis tight;
            xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 0.7]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
            
            % --- Difference (Pref - Null) ---
            PSTH_all_raw_PrefminusNull_tmp{1} = PSTH_all_raw{1}(methods_of_select{ms,1},:,1:2:end)...
                - PSTH_all_raw{1}(methods_of_select{ms,1},:,2:2:end);
            PSTH_all_raw_PrefminusNull_tmp{2} = PSTH_all_raw{2}(methods_of_select{ms,1},:,1:2:end)...
                - PSTH_all_raw{2}(methods_of_select{ms,1},:,2:2:end);
            
            SeriesComparison({PSTH_all_raw_PrefminusNull_tmp{1}, PSTH_all_raw_PrefminusNull_tmp{2}},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',mat2cell(colors,ones(3,1)),'LineStyles',{'-'},...
                'ErrorBar',6,'Xlabel',[],'Ylabel','Raw firing','axes',h_subplot((ms-1)*2+2),...
                'CompareIndex',[1:3,1,2;1:3,3,3],...
                'CompareColor',[mat2cell(colors,ones(3,1));colors(1,:);colors(2,:);colors(3,:)],...
                'Transparent',transparent);
            
            linearSum = sum(mean(PSTH_all_raw_PrefminusNull_tmp{1}(:,:,1:2),1),3);
            hold on; plot(rate_ts{1}(rate_ts{1}<=1500),linearSum(rate_ts{1}<=1500),'m-','linew',2)
            
            %             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            xlabel('Time (ms)');
            legend off;
            
            axis tight;
            xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);   % ylim([-0.1 0.4]);
            
            % Gaussian vel
%             plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/6,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
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
        Weights_Property_Correlation( weights_straight_mean /sum(weights_straight_mean) ,...
            {'Weighted by straight mean','Weighted by straight mean'},select_bottom_line);
        
        % Norm by dynamic range
        weights_norm = weights_normalized_PSTH(select_bottom_line)';
        Weighted_sum_PSTH( weights_norm,{'Weighted by 1/dynamic range'},select_bottom_line);
        title('Weighted by 1/dynamic range');
        Weights_Property_Correlation(weights_norm/sum(weights_norm),...
            {'Weighted by 1/dynamic range','Weighted by 1/dynamic range'},select_bottom_line);
        
        % Norm by dynamic range with cell selection
        weights_norm_tcells = weights_norm;
        weights_norm_tcells(~select_tcells(select_bottom_line)) = 0;
        Weighted_sum_PSTH( weights_norm_tcells,{'Weighted by 1/dynamic range with cell selection'},select_bottom_line);
        title('Weighted by 1/dynamic range with cell selection');
        Weights_Property_Correlation(weights_norm_tcells/sum(weights_norm_tcells),...
            {'Weighted by 1/dynamic range with cell selection','Weighted by 1/dynamic range with cell selection'},select_bottom_line);
        
    end

    %% =================================================================================================================
    % For f1p1p6 & f1p1p7. @HH20160906 
    % All angles, Pref - null 
    PSTH_all_raw_PrefminusNull{1} = PSTH_all_raw{1}(select_bottom_line,:,1:2:end)...
        - PSTH_all_raw{1}(select_bottom_line,:,2:2:end);
    PSTH_all_raw_PrefminusNull{2} = PSTH_all_raw{2}(select_bottom_line,:,1:2:end)...
        - PSTH_all_raw{2}(select_bottom_line,:,2:2:end);
    [single_cell_plot_order_index,single_cell_plot_order] = sort(mean(abs(Choice_pref_all(:,select_bottom_line(select_cpref_mpref),3))),'descend'); % Stim-on to stim-off
    
    function f1p1p6(debug)
        if debug  ; dbstack;   keyboard;      end
        %%
        % --- Difference (Pref - Null) ---
        
        % 8-degree      
%         PSTH_all_raw_PrefminusNull{1} = squeeze(PSTH_correct_angles_raw{1}(select_bottom_line,:,7,:)...
%             - PSTH_correct_angles_raw{1}(select_bottom_line,:,8,:));
%         PSTH_all_raw_PrefminusNull{2} = squeeze(PSTH_correct_angles_raw{2}(select_bottom_line,:,7,:)...
%             - PSTH_correct_angles_raw{2}(select_bottom_line,:,8,:));
        
%         [~,plot_order] = sort(mean(abs(Choice_pref_all(:,select_bottom_line(select_cpref_mpref),1))),'descend');
        
        %%{
        % === All traces ===
        
        
        set(figure(700),'name','Single cell (Correct only, all choices)','pos',[27 63 1449 892]); clf
        [~,h_subplot] = tight_subplot(4,5,[0.03 0.02]);
        counter = 0;
        ps = nan(3,length(rate_ts{1}),length(single_cell_plot_order));
        in_larger_than_out = ps;

        for nn = 1:length(single_cell_plot_order)
            counter = counter + 1;
            
            if counter > 20
                counter = 1;
                set(figure(700+fix(nn/20)),'name','Single cell (Correct only, all choices)','pos',[27 63 1449 892]); clf
                [~,h_subplot] = tight_subplot(4,5,[0.03 0.02]);
            end
            
            cellNo_in_find_bottom_line = single_cell_plot_order(nn);
            cellNo_origin = find(cumsum(select_bottom_line)==cellNo_in_find_bottom_line,1);
            
%             axes(h_subplot(nn));
            for j = 1:2
                ys_this{j} = group_result(cellNo_origin).mat_raw_PSTH.PSTH{j,1,1}.ys';
                sem_this{j} = group_result(cellNo_origin).mat_raw_PSTH.PSTH{j,1,1}.sem';
                ps_this{j} = group_result(cellNo_origin).mat_raw_PSTH.PSTH{j,1,1}.ps';
            end

            SeriesComparison({shiftdim(ys_this{1},-1) shiftdim(ys_this{2},-1)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'OverrideError',{sem_this{1}, sem_this{2}},...
                'OverridePs',{ps_this{1}, ps_this{2}},'ErrorBar',6,...
                'CompareIndex',[1 3 5;2 4 6],'CompareColor',{colors(1,:),colors(2,:),[0 0.8 0.4]},...
                'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'axes',h_subplot(counter),...
                'Transparent',transparent);

            
            title(h_subplot(counter), sprintf(' %g (ori. #%g, ps: %g, %g, %g; %g)',nn,cellNo_origin,Choice_pref_p_value_all(:,cellNo_in_find_bottom_line,3)',single_cell_plot_order_index(nn)));  % 3, stim-on to stim-off
            legend off;
            axis(h_subplot(counter),'tight');
            ylabel(h_subplot(counter), '');

            % Cache for calcuating the first significant times. HH20161108
            ps(:,:,nn) = ps_this{1}';
            in_larger_than_out(:,:,nn) = (ys_this{1}(:,1:2:end) > ys_this{1}(:,2:2:end))';
            
            % Annotate typical cell
            if select_tcells(cellNo_origin)
                set(gca,'color',hsv2rgb([1 0.1 1]));
            end
            
            % Plot mem-sac trace
            set(gca,'ButtonDownFcn',{@Plot_HD,cellNo_origin});
            
            if counter ~= 20
                set(h_subplot(counter),'xticklabel','');
            end
            
            plot(h_subplot(counter),xlim,[0 0],'k--');
            
        end
        %} 
    end       

    firstDivergenceTime = [];
    function f1p1p6p1(debug)  % Divergence time
        if debug  ; dbstack;   keyboard;      end

        %% === Divergence time ===
        %%%%%%%%%%%%%%%%%%%%%%%
        min_bins_for_1st_sign_time = 25; % 250 ms
        p_critical_for_1st_sign_time = 0.05;
        t_valid_ind = rate_ts{1}>= 200 & rate_ts{1}<=1500; % Only sensory phase
        %%%%%%%%%%%%%%%%%%%%%%%
        
        firstDivergenceTime = nan(length(select_bottom_line),3);
        t_valid = rate_ts{1}(t_valid_ind);
        
        for nn = 1:length(select_bottom_line)
            ys_this = group_result(nn).mat_raw_PSTH.PSTH{1,1,1}.ys'; % j = 1
            ps_this = group_result(nn).mat_raw_PSTH.PSTH{1,1,1}.ps;  % j = 1
            in_larger_than_out = (ys_this(:,1:2:end) > ys_this(:,2:2:end))';

            ps_valid = double((ps_this(:,t_valid_ind) < p_critical_for_1st_sign_time) & in_larger_than_out(:,t_valid_ind));
            
            conv_tmp = conv2(ps_valid,ones(1,min_bins_for_1st_sign_time));
            for k = 1:3
                % This is really brilliant :)
                first_tmp = find(conv_tmp(k,:) == min_bins_for_1st_sign_time,1) - min_bins_for_1st_sign_time + 1;  

                if ~isempty(first_tmp)
                    firstDivergenceTime(nn,k) = t_valid(first_tmp);
                end
            end
        end
        
        % Paired, nans ignored for each two
        first_sign_result = BarComparison(firstDivergenceTime,'figN',556,'Colors',mat2cell(colors,ones(stim_type_num,1)),'PairTTest',1);
        first_sign_ps =  [first_sign_result.ps_ttest(2,3),first_sign_result.ps_ttest(2,4),first_sign_result.ps_ttest(3,4)];
        first_sign_ns =  [first_sign_result.ps_paired_ns(2,3),first_sign_result.ps_paired_ns(2,4),first_sign_result.ps_paired_ns(3,4)];
%         title(sprintf('%g, ',sum(~isnan(firstDivergenceTime)),first_sign_ps));
        title(sprintf('Paired for each two: %g,%g,%g,%g,%g,%g ',first_sign_ns,first_sign_ps));
        ylim([0 1500]); SetFigure(15);
        view([90 90]);
        
        % Paired, only all ~isnan cells are included
        firstDivergenceTime_allSignif = firstDivergenceTime(all(~isnan(firstDivergenceTime),2),:);
        first_sign_result_allchoice = BarComparison(firstDivergenceTime_allSignif,'figN',557,'Colors',mat2cell(colors,ones(stim_type_num,1)),'PairTTest',1);
        first_sign_ps =  [first_sign_result_allchoice.ps_ttest(2,3),first_sign_result_allchoice.ps_ttest(2,4),first_sign_result_allchoice.ps_ttest(3,4)];
        title(sprintf('Paired for "All choice": %g, %g, %g, %g, %g, %g ',sum(~isnan(firstDivergenceTime_allSignif)),first_sign_ps));
        ylim([0 1500]); SetFigure(15);
        view([90 90]);

        % Paired, nans ignored for each two
        first_sign_result = BarComparison(firstDivergenceTime,'figN',558,'Colors',mat2cell(colors,ones(stim_type_num,1)),'PairTTest',0);
        first_sign_ps =  [first_sign_result.ps_ttest(2,3),first_sign_result.ps_ttest(2,4),first_sign_result.ps_ttest(3,4)];
%         title(sprintf('%g, ',sum(~isnan(firstDivergenceTime)),first_sign_ps));
        title(sprintf('Unpaired for "any chioce": %g,%g,%g,%g,%g,%g ',sum(~isnan(firstDivergenceTime)),first_sign_ps));
        ylim([0 1500]); SetFigure(15);
        view([90 90]);
        
        
    end

    %%
    function f1p1p7(debug)
        if debug  ; dbstack;   keyboard;      end
        
        % === Curve fitting ===
        set(figure(800),'name','Single cell fittings (Correct only, all choices)','pos',[27 63 1449 892]); clf
        [~,h_subplot] = tight_subplot(4,5,[0.03 0.02]);
        counter = 0;
        
        % Define sigmoid function for fitting
        func_sigmoid = fittype('aa*1/(1+exp(-(x-tt)/bb))','coeff',{'aa','bb','tt'});
        options = fitoptions(func_sigmoid);
        options.StartPoint = [10 100 600];
        
        aas = nan(length(single_cell_plot_order),stim_type_num);
        bbs = aas; tts = aas;
        
        for nn = 1:length(single_cell_plot_order)
            counter = counter + 1;
            if counter > 20
                counter = 1;
                set(figure(800+fix(nn/20)),'name','Single cell fittings (Correct only, all choices)','pos',[27 63 1449 892]); clf
                [~,h_subplot] = tight_subplot(4,5,[0.03 0.02]);
            end
            
            nnn = single_cell_plot_order(nn);
            
            axes(h_subplot(counter));

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
        
        %% === Statistics
        select_cells = select_tcells;  exclude_cells = setdiff(1:size(bbs,1),select_cells);
        result_bb = BarComparison(bbs(select_cells,:),'figN',556,'Colors',mat2cell(colors,ones(stim_type_num,1)));
        title('bbs'); ylim([1 200]); plot(1:3,bbs(exclude_cells,:),'o:','color',[0.7 0.7 0.7]);
        
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
            % select_bottom_line, 'All cells';
            select_tcells, 'Typical cells';
            };
                
        for ms = 1:size(methods_of_select,1)
            set(figure(999-ms),'name',['Average PSTH (correct only, all headings), ' methods_of_select{ms,2}],'pos',[27 57 919 898]); clf
            h_subplot = tight_subplot(3,2,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
            linkaxes(h_subplot(1:3),'xy')
            linkaxes(h_subplot(4:6),'xy')
            
            for k = 1:3
                
                colors_angles = colormap(gray); 
                colors_angles = ones(size(colors_angles,1),3) - colors_angles .* repmat([1 1 1]-colors(k,:),size(colors_angles,1),1);
                colors_angles = colors_angles(round(linspace(20,length(colors_angles),5)),:);
                colors_angles = reshape(repmat(colors_angles,1,2)',3,[])';
                colors_angles = mat2cell(colors_angles,ones(10,1));
                
                % --- Ramping with different angles ---
                for j = 1:2
%                      yyy{j} = PSTH_correct_angles_Norm{j}(methods_of_select{ms,1},:,:,k);
                    yyy{j} = PSTH_correct_angles_raw{j}(methods_of_select{ms,1},:,:,k);
                    ttt{j} = rate_ts{j};
                    
                    yyy_diff{k}{j} =  yyy{j}(:,:,1:2:end) - yyy{j}(:,:,2:2:end);
                end
                
                
                h = SeriesComparison(yyy,{ttt{1} ttt{2} time_markers},...
                    'Colors',colors_angles,'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                if k == 1
                    title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 max(ylim)]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
                
                % ----  Difference ---
                
                h = SeriesComparison(yyy_diff{k},{ttt{1} ttt{2} time_markers},...
                    'Colors',colors_angles(1:2:end),'LineStyles',{'-'},...
                    'ErrorBar',2,'Xlabel',[],'Ylabel','Diff','axes',h_subplot(k+(2-1)*3));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]); ylim([-.1 max(ylim)]); plot(xlim,[0 0],'k--');
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
            end
            
            axis(h_subplot(end),'tight');
            SetFigure(15); drawnow;
            
            
            % --- Correct and Wrong trials of different angles --- HH20180711 for GuYong's Hangzhou ppt
            
            set(figure(221816-ms),'name',['Average PSTH (correct + wrong, all headings), ' methods_of_select{ms,2}]); clf
            h_subplot = tight_subplot(3,3,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
            set(gcf,'uni','norm','pos',[0.009       0.059       0.732       0.853]);
            linkaxes(h_subplot,'xy')
            
            for k = 1:3
                
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
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                if k == 1
                    title([methods_of_select{ms,2} ', Correct only, n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 max(ylim)]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
                
                % --- Ramping with different angles, Correct + Wrong ---
               
                colors_angles = colormap(gray);
                colors_angles = ones(size(colors_angles,1),3) - colors_angles .* repmat([1 1 1]-colors(k,:),size(colors_angles,1),1);
                colors_angles = colors_angles(round(linspace(20,length(colors_angles),5)),:);
                colors_angles = [flipud(colors_angles); colors_angles(2:end,:)];
                
                for j = 1:2
                    yyy{j} = PSTH_correctNwrong_angles_Norm{j}(methods_of_select{ms,1},:,:,k);
%                     yyy{j} = PSTH_correctNwrong_angles_raw{j}(methods_of_select{ms,1},:,:,k);
                end
                
                h = SeriesComparison(yyy,{rate_ts{1} rate_ts{2} time_markers},...
                    'Colors',colors_angles,'LineStyles',{'--','--','--','--','-','-','-','-','-'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k+(2-1)*3));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                if k == 1
                    title([methods_of_select{ms,2} ', Correct + Wrong, n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  %ylim([0.1 0.8]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                
                % --- Ramping with different angles, Wrong only ---
                % Has been changed to [small pref, small null, ..., large pref, large null]
                
                for j = 1:2
                    yyy{j} = PSTH_wrong_angles_Norm{j}(methods_of_select{ms,1},:,:,k);
%                     yyy{j} = PSTH_wrong_angles_raw{j}(methods_of_select{ms,1},:,:,k);
                end
                
                % Note that for wrong only, here I align to the "stimulus" not the "choice"
                h = SeriesComparison(yyy,{rate_ts{1} rate_ts{2} time_markers},...
                    'Colors',colors_angles([4 4 3 3 2 2 1 1],:),'LineStyles',{'--','-'},...  
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k+(3-1)*3));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                if k == 1
                    title([methods_of_select{ms,2} ', Wrong only, n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  %ylim([0.1 0.8]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                
                legend off;
            end
            
            SetFigure(15); drawnow;
            axis(h_subplot(3),'tight');
            
            
            % --- Multisensory enhancement of different angles --- HH20160126
            set(figure(995-ms),'name',['Enhancement for different angles (correct only, all headings), ' methods_of_select{ms,2}],'pos',[218 35 1674 928]); clf
            h_subplot = tight_subplot(2,3,[0.1 0.05],[0.05 0.07],[0.1 0.05]);
            h_subplot =  reshape(reshape(h_subplot,2,3)',[],1); delete(h_subplot(end));
            unique_abs_heading = unique(abs(group_result(representative_cell).mat_raw_PSTH.heading_per_trial));
            linkaxes(h_subplot,'xy')
            
            for aa = 1:size(yyy_diff{1}{1},3) % Different angles
                
                % === Plotting ===
                yyy_diff_this_angle_all_stim_type = [];
                for j = 1:2
                    for k = 1:3
                        yyy_diff_this_angle_all_stim_type{j}(:,:,k) = yyy_diff{k}{j}(:,:,aa); % Reorganize data
                    end
                end
                
                h = SeriesComparison(yyy_diff_this_angle_all_stim_type, {ttt{1} ttt{2} time_markers},...
                    'Colors',mat2cell(colors,ones(3,1)),'LineStyles',{'-'},...
                    'ErrorBar',2,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(aa),'Transparent',0);
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]); ylim([-.1 .5]); plot(xlim,[0 0],'k--');
                
               
                title(['|heading| = ',num2str(unique_abs_heading(aa))]);
                if k == 1
                    title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                if aa ~= 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end

                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                legend off;     drawnow
                
                % === Linear fitting combined trace. HH20160510 ===
                ts = rate_ts{1};
                t_select = rate_ts{1}> 0 & rate_ts{1} <=1500;
                
                ramping1(:,aa) = h.means{1}(1,t_select);
                ramping2(:,aa) = h.means{1}(2,t_select);
                ramping3(:,aa) = h.means{1}(3,t_select);
                
                w = fminsearch(@(w) sum((w(1)*ramping1(:,aa) + w(2)*ramping2(:,aa) - ramping3(:,aa)).^2), [.5 .5])
                
                plot(ts(t_select),ramping1(:,aa)*w(1)+ramping2(:,aa)*w(2),'k','linew',2);
                plot(ts(t_select),ramping1(:,aa)*sqrt(2)/2+ramping2(:,aa)*sqrt(2)/2,'m--','linew',2);
                text(700,0,num2str(w));
            end
            axis(h_subplot(end-1),'tight');
            
            % == All angle use the same weight. HH20160926 ==
            w_all = fminsearch(@(w) sum((w(1)*ramping1(:) + w(2)*ramping2(:) - ramping3(:)).^2), [.5 .5])
            for aa = 1:size(yyy_diff{1}{1},3) % Different angles
                axes(h_subplot(aa));
                plot(ts(t_select),ramping1(:,aa)*w_all(1)+ramping2(:,aa)*w_all(2),'c','linew',2);  
                text(700,0.1,num2str(w_all),'color','c');
            end

            SetFigure(15); 
        end
        
    end

    function f1p1p4(debug) % Plot Example cells for Supplementary Figure 2
        if debug  ; dbstack;   keyboard;      end
        
        to_plot_ori = {
            % Four example cells in Figure 1
        %{
            [
            61 % m5c313r2u05
            126 % m10c133r3_4_7
            31 % m5c174r1
            120
            ]; % m10c121r3
        %}
            % Other exmples in SuppleFig 2
%         %{
            [ 
            162
            152
            91
            % 114 % Same as 91
            95
            % 128
            29
            149
            6
            74
            130
            
            ]
          %}
            [  % Bad cells
            70
            134
            58
            103
            103
            77
            ]
       
                       
            }; % Manually selected
        
        toPlotFigN = length(to_plot_ori);
        for ff = 1:toPlotFigN
            
            toPlotN = length(to_plot_ori{ff});
            
            figure(5837+ff); clf; % A large figure
            set(gcf,'uni','norm','pos',[0.001   0.04    0.626    0.78/4*toPlotN]);
            hs = tight_subplot(toPlotN,5,[0.06*4/toPlotN 0.06],[0.1 0.1]*4/toPlotN, [0.05 0.05]);
            
            hs = reshape(hs,[],5);
            
            for ee = 1:toPlotN
                Plot_HD([],[],to_plot_ori{ff}(ee),hs(ee,:));
            end
            
            % set(hs(1:end-1,:),'xtick','')
        end
    end


    function f1p2p5(debug)      % Rate 2.5. Different headings, Gu's idea, Enhance index (comb/max ratio) distribution through time
        if debug
            dbstack;
            keyboard;
        end
        
        methods_of_select = {
            select_tcells, 'Typical cells';
        };
        
        unique_abs_heading = unique(abs(group_result(representative_cell).mat_raw_PSTH.heading_per_trial));
        
        for ms = 1:1
            
            j = 1;
            ttt = rate_ts{j};
            
            % Organize data
            yyy_diff_angle_cellBtBheadingBk{j} =  PSTH_correct_angles_raw{j}(methods_of_select{ms,1},:,1:2:end,:)...
                                                - PSTH_correct_angles_raw{j}(methods_of_select{ms,1},:,2:2:end,:);
            yyy_diff_hard_easy_cellBtB2Bk{j} = PSTH_hard_easy_raw_cellBtB4Bk{j}(methods_of_select{ms,1},:,1:2:end,:)...
                                             - PSTH_hard_easy_raw_cellBtB4Bk{j}(methods_of_select{ms,1},:,2:2:end,:);
            
            % Sliding window
            sliding_width = 300; % ms, overlap allowed
            sliding_step = 200;
            ttt_1 = -100;
            sliding_n = ceil((ttt(end)-ttt_1-sliding_width)/sliding_step);
            
            % Hist range
            hist_range = (-1:0.1:4) + 0.05;
            plot_range = [-1.1 4.15];
                
            % -- Grouped by headings --
            for aa = 1:size(yyy_diff_angle_cellBtBheadingBk{j},3) % Different angles
                                
                set(figure(700+aa),'name',['Enhance index, heading = ' num2str(unique_abs_heading(aa)) ', ' methods_of_select{ms,2}],'pos',[2 379 1674 581]); clf
                h_subplot = tight_subplot(1,sliding_n,[0,0.02],[0.2,0.1]);
                maxmaxy = 0;
                
                for ss = 1:sliding_n % Sliding window
                    
                    EI_t_range = (ss-1)*sliding_step + ttt_1 <= ttt & ttt <= (ss-1)*sliding_step + ttt_1 + sliding_width;
                    EI_ts = mean(ttt(EI_t_range));
                    
                    yyy_diff_this_angle_this_bin_cellBstimtype = ...
                        squeeze(mean(yyy_diff_angle_cellBtBheadingBk{j}(:,EI_t_range,aa,:),2));
                
                    EI_this_bin = yyy_diff_this_angle_this_bin_cellBstimtype(:,3)./...,
                        max(yyy_diff_this_angle_this_bin_cellBstimtype(:,1),...
                        yyy_diff_this_angle_this_bin_cellBstimtype(:,2));
%                     EI_this_bin = yyy_diff_this_angle_this_bin_cell_by_stimtype(:,3)./...,
%                         (yyy_diff_this_angle_this_bin_cell_by_stimtype(:,2)+...
%                         0*yyy_diff_this_angle_this_bin_cell_by_stimtype(:,2));  
                    
                    median_EI_this_bin = nanmedian(EI_this_bin);
                    mean_EI_this_bin = nanmean(EI_this_bin);
                    [~,p] = ttest(EI_this_bin,1);
                    
                    axes(h_subplot(ss));
                    hist(EI_this_bin,hist_range);    hold on; axis tight;
                    view(270,90); xlim(plot_range); 
                    
                    set(gca,'xaxislocation','top');
                    
                    if ss<sliding_n
                        %set(gca,'xticklabel',[],'yticklabel',[]);
                        set(gca,'xticklabel',[]);
                    end
                    
                    ylabel(sprintf('%g ~ %g',ttt(find(EI_t_range,1)),ttt(find(EI_t_range,1,'last')))); 
                    title(sprintf('med %.3g\nmean %.3g\np = %.2g',median_EI_this_bin,mean_EI_this_bin,p));
                    
                    maxmaxy = max(maxmaxy,max(ylim));
                    
                    plot([median_EI_this_bin median_EI_this_bin],[0 20],'r-','linew',2);
                    plot([mean_EI_this_bin mean_EI_this_bin],[0 20],'g-','linew',2);
                    
                    % Save EI
                    EI_all_angles_cellBtimeBheading(:,ss,aa) = EI_this_bin;
                end
                
                for ss = 1:sliding_n
%                     ylim(h_subplot(ss),[0 maxmaxy]);
                end
                
                SetFigure(15);
                
            end
            
            % -- Grouped by hard/easy --
            for dd = 1:2
                maxmaxy = 0;
                set(figure(750+dd),'name',['Enhance index, [hard, easy] = ' num2str(dd) ', ' methods_of_select{ms,2}],'pos',[2 379 1674 581]); clf
                h_subplot = tight_subplot(1,sliding_n,[0,0.02],[0.2,0.1]);
                
                for ss = 1:sliding_n % Sliding window
                    
                    EI_t_range = (ss-1)*sliding_step + ttt_1 <= ttt & ttt <= (ss-1)*sliding_step + ttt_1 + sliding_width;
                    
                    yyy_diff_this_difficulty_this_bin_cellBstimtype = ...
                        squeeze(mean(yyy_diff_hard_easy_cellBtB2Bk{j}(:,EI_t_range,dd,:),2));
                    EI_this_bin = yyy_diff_this_difficulty_this_bin_cellBstimtype(:,3)./...,
                        max(yyy_diff_this_difficulty_this_bin_cellBstimtype(:,1),...
                        yyy_diff_this_difficulty_this_bin_cellBstimtype(:,2));

                    median_EI_this_bin = nanmedian(EI_this_bin);
                    mean_EI_this_bin = nanmean(EI_this_bin);
                    [~,p] = ttest(EI_this_bin,1);
                    
                    axes(h_subplot(ss));
                    hist(EI_this_bin,hist_range);    hold on; axis tight;
                    view(270,90); xlim(plot_range);
                    set(gca,'xaxislocation','top');
                    
                    if ss<sliding_n
                        %set(gca,'xticklabel',[],'yticklabel',[]);
                        set(gca,'xticklabel',[]);
                    end

                    ylabel(sprintf('%g ~ %g',ttt(find(EI_t_range,1)),ttt(find(EI_t_range,1,'last')))); 
                    title(sprintf('med %.3g\nmean %.3g\np = %.2g',median_EI_this_bin,mean_EI_this_bin,p));
                    
                    maxmaxy = max(maxmaxy,max(ylim));
                    
                    plot([median_EI_this_bin median_EI_this_bin],[0 20],'r-','linew',2);
                    plot([mean_EI_this_bin mean_EI_this_bin],[0 20],'g-','linew',2);

                    % Save EI
                    EI_all_hardeasy_cellBtimeB2(:,ss,dd) = EI_this_bin;
                   
                end
                
                for ss = 1:sliding_n
%                     ylim(h_subplot(ss),[0 maxmaxy]);
                end
                
                SetFigure(15);               
            end
            
            % -- Correlation between hard and easy EIs ---
            set(figure(750+dd+1),'unit','pix','name',['Enhance index, [hard, easy] = ' num2str(dd) ', ' methods_of_select{ms,2}],'pos',[2 49 1250 911]); clf
            [~,h_subplot] = tight_subplot(fix(sqrt(sliding_n)),ceil(sliding_n/fix(sqrt(sliding_n))),[0.05,0.02],[0.05,0.05]);
            
            for ss = 1:sliding_n % Sliding window
                
                % Easy vs Hard
                EI_t_range = (ss-1)*sliding_step + ttt_1 <= ttt & ttt < (ss-1)*sliding_step + ttt_1 + sliding_width;
                xx = EI_all_hardeasy_cellBtimeB2(:,ss,1);
                yy = EI_all_hardeasy_cellBtimeB2(:,ss,2);
                
                [~,p] = ttest(xx,yy);
                
                valid_EI = xx >0 & xx < 5 & yy >0 & yy < 5; % Just work-around, should select more carefully. @HH20160908
                h = LinearCorrelation(...
                { xx(valid_EI), xx(~valid_EI)},...
                { yy(valid_EI), yy(~valid_EI)},...
                    'Axes',h_subplot(ss),'Xlabel',sprintf('Hard (%g ~ %g), p= %g, x-y= %g',ttt(find(EI_t_range,1)),ttt(find(EI_t_range,1,'last')),p,mean(xx)-mean(yy)),...
                    'Ylabel','Easy','MarkerSize',5,'FaceColors',{'k','none'});
                delete([h.group(2).line]);
                axis([-1 5 -1 5]); axis square; legend off;
                hold on;
                plot(xlim,ylim,'--k');
                
                % Mean + sem
                mean_xx = mean(xx);
                sem_xx = std(xx)/sqrt(length(xx));
                mean_yy = mean(yy);
                sem_yy = std(yy)/sqrt(length(yy));
                plot([mean_xx-sem_xx mean_xx+sem_xx],[mean_yy mean_yy],colors(2,:));
                plot([mean_xx mean_xx],[mean_yy-sem_yy mean_yy+sem_yy],colors(2,:));
                
                % Show individual cell selected from the figure. HH20150424
                h_line = plot(xx,yy,'visible','off'); hold on;
                set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, methods_of_select{ms,1}});

                % EZZoom
                uicontrol('Style','pushbutton','String','tight','unit','norm','Position',[0.024 0.944 0.067 0.034],...
                    'callback','axis(findall(gcf,''type'',''axes''),''tight'')');
                uicontrol('Style','pushbutton','String','[-1 5]','unit','norm','Position',[0.096 0.944 0.067 0.034],...
                    'callback','axis(findall(gcf,''type'',''axes''),[-1 5 -1 5])');
                
            end
        end                
    end

    tuning_pack = [];

    function f1p2p7(debug)      % Rate 2.7. Tuning curve analysis. (Dora's tuning) 20161019
        if debug
            dbstack;
            keyboard;
        end
                
        %% Tuning: At 3 different phases. HH20150414
        % Now this part is dirty and quick.
        % Here are several things to be done:
        %  1. Align to stimulus onset AND saccade onset  (done)
        %  2. Deal with different heading angles
        %  3. Put this part elsewhere
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
        select_for_tuning = select_tcells; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % find_for_tuning = find(select_for_tuning);
        
        
        any_cell = any(group_ChoicePreference_pvalue(:,:,3) < 0.05,1)';
        
        for k = 1:3
            tuning_cell = (group_ChoicePreference_pvalue(k,:,3) < 0.05)';
            tmp = any_cell + tuning_cell;
            tmp(tmp == 0) = [];
            tuning_cell_in_any_cell{k} = find(tmp == 2);
        end
        
        find_for_tuning = find(any_cell);
   
        
        % Time centers (of CP time windows, width = 500 ms)
%         t_stim_center_earlier = mean(time_markers{j}(1,1:2)) -100; % Stim center + sensory delay
%         t_stim_center = mean(time_markers{j}(1,1:2)) + 150; % Stim center + sensory delay
%         t_pre_center = mean(time_markers{j}(1,3)) - group_result(representative_cell).mat_raw_PSTH.binSize_CP/2;   % Pre-sac epoch
%         t_post_center = mean(time_markers{j}(1,3)) + group_result(representative_cell).mat_raw_PSTH.binSize_CP/2;  % Post-sac epoch
        
%         t_stim_center = 550; % Vestibular center
%         t_pre_center = 1100;   % Visual center
%         t_post_center = 1500;  % Post-sac epoch
%         
% %         [~,center_t_ind_earlier] = min(abs(CP_ts{j} - t_stim_center_earlier));
%         [~,center_t_ind] = min(abs(CP_ts{j} - t_stim_center));
%         [~,pre_t_ind] = min(abs(CP_ts{j} - t_pre_center));
%         [~,post_t_ind] = min(abs(CP_ts{j} - t_post_center));

        %         tuning_time_center = [605 750 1100 1500];
        % tuning_time_center = [605 1000 1200 1500];
        tuning_time_center = [655 950 1200 1500];
        
        for tt = 1:length(tuning_time_center)
            [~,tmp] = min(abs(CP_ts{j} - tuning_time_center(tt)));
            tuning_time_phase(tt) = tmp;
        end
        tuning_time_phase_title = {'Vest largest', 'Comb largest','Vis largest', 'Pre-sac'};
        
        % Suppose all the unique_headings are the same!!
        unique_heading = group_result(representative_cell).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
        % HH20160214: Patch calculation for LEFT & RIGHT choices of 0 headings
        zero_index = find(unique_heading == 0);
        unique_heading_two_zeros = [unique_heading(1:zero_index-1); -eps; eps ;unique_heading(zero_index+1:end)];  
%         unique_heading_two_zeros = [unique_heading(1:zero_index-1) ;unique_heading(zero_index+1:end)];  % No zero version

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for tuning curve normalization
        
%         time_for_dynamic_range = find(CP_ts{j} < time_markers{j}(1,2)); % all pre-sac times
        time_for_dynamic_range = [tuning_time_phase(1) - 3 tuning_time_phase(1) + 3]; % Around stimulus center
        modalities_share_same_dynamic_range = 1;
        linear_or_sigmoid = 1; % 1: linear fitting; 2: sigmoid fitting
        min_trials_dora_tuning = 3;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isempty(tuning_pack) % Only run once if needed
            
            progressbar()
            % Pack tuning data
            for i = 1:sum(select_for_tuning)  % For cells
                
                progressbar(i/sum(select_for_tuning))
                
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
                    
                    
                    for k = 1:3
                        this_raw = group_result(find_for_tuning(i)).mat_raw_PSTH.CP{j,k}.raw_CP_result{pp};
                        this_tuning = this_raw.Neu_tuning(:,2);
                        this_tuning_correctonly = this_raw.Neu_tuning_correctonly(:,2);
                        this_tuning_dora = this_raw.Neu_tuning_Dora_matrix_mean';  % HH20161019
                        this_tuning_dora_n = this_raw.Neu_tuning_Dora_matrix_n';
                        
                        % @HH20160214: Patch calculation for LEFT & RIGHT choices of 0 headings (because I failed to do this
                        % in the original CP_HH and I feel hesitant to redo all the batch files from A to Z right now...).
                        % Note, sadly, that this could be NaN because I also failed to pack the spike counts into
                        % spike_counts_allheadings_grouped if the number of left/right choices were fewer than 3... WTF...
                        
                        % @HH20160907. Today I redo all the batch files finally!
                        
                        % Discard dora tuning which has less than "min_trials_dora_tuning" trials. HH20161109
                        this_tuning_dora(this_tuning_dora_n<min_trials_dora_tuning) = nan;
                        
                        if length(this_tuning) == length(unique_heading)
                            zero_spikes_left_right = [mean(this_raw.spike_counts_allheadings_grouped{zero_index, 1}); mean(this_raw.spike_counts_allheadings_grouped{zero_index, 2})];
                            this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); zero_spikes_left_right; this_tuning_correctonly(zero_index+1:end)];
                        else % No zero heading at all in the file
                            zero_spikes_left_right = [nan;nan];
                            this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); zero_spikes_left_right; this_tuning_correctonly(zero_index:end)];
                            this_tuning = [this_tuning(1:zero_index-1); nan ;this_tuning(zero_index:end)];
                            this_tuning_dora = [this_tuning_dora(1:zero_index-1,:); nan nan; this_tuning_dora(zero_index:end,:)];
                        end
                        
                        %                     this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); this_tuning_correctonly(zero_index+1:end)]; % No zero version
                        
                        % Align preferred direcions to 'rightward'. So if PREF == 1 for a certain cell, we
                        % flip all its tuning curve.
                        if flip_tuning
                            this_tuning = flipud(this_tuning);
                            this_tuning_correctonly = flipud(this_tuning_correctonly);
                            
                            % For the dora tuning matrix, if PREF == left, should flip both lr and ud.HH20161019
                            %               -    0    +                             null   0    pref
                            %   L choice    xx   xx   nan     ==>   null choice      xx     xx    nan
                            %   R choice    nan  xx  xx             pref choice      nan    xx    xx
                            this_tuning_dora = flipud(this_tuning_dora);
                            this_tuning_dora = fliplr(this_tuning_dora);
                        end
                        
                        % Update dynamic range
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
                            tuning_pack{3}{k,pp}(:,:,i) = this_tuning_dora';  % HH20161019
                        catch
                            tuning_pack{1}{k,pp}(:,i) = nan;
                            tuning_pack{2}{k,pp}(:,i) = nan;
                            tuning_pack{3}{k,pp}(:,:,i) = nan;
                            keyboard;
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
                            
                            tuning_pack{3}{k,pp}(:,:,i) = tuning_pack{3}{k,pp}(:,:,i) - offset;
                            tuning_pack{3}{k,pp}(:,:,i) = tuning_pack{3}{k,pp}(:,:,i) / gain;  % HH20161019
                        else
                            % Three modalities have their own dynamic ranges
                            tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) - offset_own_range(k);
                            tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) / gain_own_range(k);
                            
                            tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) - offset_own_range(k);
                            tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) / gain_own_range(k);
                            
                            tuning_pack{3}{k,pp}(:,:,i) = tuning_pack{3}{k,pp}(:,:,i) - offset_own_range(k);  % HH20161019
                            tuning_pack{3}{k,pp}(:,:,i) = tuning_pack{3}{k,pp}(:,:,i) / gain_own_range(k);
                        end
                    end
                end
                
            end
        end
        
        % Get means and sems
        tuning_mean_all = nan(length(unique_heading),length(CP_ts{j}),3);
        tuning_sem_all = tuning_mean_all;
        tuning_mean_correctonly = nan(length(unique_heading_two_zeros),length(CP_ts{j}),3);
        tuning_sem_correctonly = tuning_mean_correctonly;
        tuning_mean_dora = nan(2,length(unique_heading),length(CP_ts{j}),3);
        tuning_sem_dora = tuning_mean_dora;
        tuning_n_dora = tuning_mean_dora;
        tuning_sig_fit_correctonly = nan(3,length(CP_ts{j}),3);
        tuning_linear_fit_correctonly = nan(4,length(CP_ts{j}),3,2);
        
        % For sigmoid fitting. HH20150528
%         sigfunc = @(sig_para, x)(sig_para(1)+ sig_para(2)./ (1 + exp(-x/sig_para(3))));
        sigfunc = @(sig_para, x)(sig_para(1)+ sig_para(2) * normcdf(x,0,sig_para(3)));
        
        function err = cost_function(q,data_cum)
            x = data_cum(:,1);
            y = data_cum(:,2);
            
            z = sigfunc(q,x);
            err = norm(z-y);
        end
        
        % Plotting tuning curves at three time points
        set(figure(145506),'name',['Dora tuning, j = ' num2str(j)]); clf;
        set(gcf,'uni','norm','pos',[0.001       0.078       0.835       0.838]);
        
        for pp = 1:length(CP_ts{j})
            
            for k = 1:3
                % Mean and sem
                % Patch calculation for LEFT & RIGHT choices of 0 headings. HH20160214
                this_tuning_all = tuning_pack{1}{k,pp}(:,~isnan(tuning_pack{1}{k,pp}(1,:)));
                tuning_mean_all(:,pp,k) = nanmean(this_tuning_all,2);
                tuning_sem_all(:,pp,k) = nanstd(this_tuning_all,[],2)./sqrt(sum(~isnan(this_tuning_all),2));
                
                % this_tuning_correctonly =  tuning_pack{2}{k,pp}(:,~isnan(tuning_pack{2}{k,pp}(1,:)));
                this_tuning_correctonly =  tuning_pack{2}{k,pp}(:,  tuning_cell_in_any_cell{k}); % HH20180612
                
                tuning_mean_correctonly(:,pp,k) = nanmean(this_tuning_correctonly,2);
                tuning_sem_correctonly(:,pp,k) = nanstd(this_tuning_correctonly,[],2)./sqrt(sum(~isnan(this_tuning_correctonly),2));;

                this_tuning_dora =  tuning_pack{3}{k,pp}(:,:,~isnan(tuning_pack{2}{k,pp}(1,:))); % Use the same NaN indicator
                tuning_mean_dora(:,:,pp,k) = nanmean(this_tuning_dora,3);
                tuning_n_dora(:,:,pp,k) = sum(~isnan(this_tuning_dora),3);
                tuning_sem_dora(:,:,pp,k) = nanstd(this_tuning_dora,[],3)./sqrt(tuning_n_dora(:,:,pp,k));
                
                
                if linear_or_sigmoid == 2   % Fitting sigmoid function. HH20150528
                    yy = nanmean(this_tuning_correctonly,2);
                    %                 tuning_sig_fit_correctonly(:,pp,k) = nlinfit(unique_heading, yy, sigfunc, [min(yy) range(yy) 1]);
                    
                    % Begin optimization
                    quick = fminsearch(@(q)cost_function(q,[unique_heading_two_zeros yy]),[min(yy) range(yy) 1]);
                    
                    % Output
                    tuning_sig_fit_correctonly(:,pp,k) = quick;
                    
                elseif linear_or_sigmoid == 1    % Linear fitting using group data. HH20161109
                    hh_index = {1:zero_index zero_index+1:size(this_tuning_correctonly,1)};
                    
                    for hh = LEFT:RIGHT % LEFT, RIGHT
                        
                        yy = this_tuning_correctonly(hh_index{hh},:)';
                        xx = unique_heading(hh_index{hh}-(hh==RIGHT))';
                        xx = repmat(xx,size(yy,1),1);
                        xx = xx(~isnan(yy)); yy = yy(~isnan(yy));
                        
                        % [r,p] = corr(xx,yy,'type','pearson'); 
                        [para,~] = polyfit(xx,yy,1);
                        
                        % tuning_linear_fit_correctonly(:,pp,k,hh) = [r p para];
                        tuning_linear_fit_correctonly([3 4],pp,k,hh) = para;
                        
                    end
                    
                    % I use the "delta norm firing as a function of |heading|"
                    % to combine the two sides. It works!   HH20180606
                    this_tuning_diff = this_tuning_correctonly(hh_index{2},:) - ...
                                       flipud(this_tuning_correctonly(hh_index{1},:));
                    yy_diff = this_tuning_diff';
                    xx_abs = repmat(unique_heading(hh_index{2}-1)',size(yy_diff,1),1);
                    xx_abs = xx_abs(~isnan(yy_diff)); yy_diff = yy_diff(~isnan(yy_diff));
                    [r,p] = corr(xx_abs(:),yy_diff(:),'type','pearson');
                    
                    % Share the same p
                    tuning_linear_fit_correctonly([1 2],pp,k,LEFT) = [r, p]';
                    tuning_linear_fit_correctonly([1 2],pp,k,RIGHT) = [r, p]';

                end

            end
            
        end
 
        % ---  Plotting ---
        for pp = 1:length(tuning_time_phase)
            
            for k = 1:3  % For each stim type
                
                % Plotting
                subplot(3,length(tuning_time_phase),pp); hold on; ylabel('All');
                
                % Traditional all heading
                %{
                plot(unique_heading,tuning_mean_all(:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
                h = errorbar(unique_heading,tuning_mean_all(:,tuning_time_phase(pp),k),tuning_sem_all(:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2);
                errorbar_tick(h,10000);
                %}

                % Dora tuning
                
                hh_index = {1:zero_index zero_index:length(unique_heading)};
                cc_marker = {'<','>'};
                
                for cc = 1:2
                    for hh = LEFT:RIGHT
                        %                 plot(unique_heading(1:zero_index)-0.2,tuning_mean_dora(1,1:zero_index,tuning_time_phase(pp),k),'v-','markersize',9,'color',colors(k,:),'markerfacecolor',colors(k,:),'LineWid',2);
    %                 plot(unique_heading(zero_index:end)+0.2,tuning_mean_dora(2,zero_index:end,tuning_time_phase(pp),k),'^-','markersize',9,'color',colors(k,:),'markerfacecolor',colors(k,:),'LineWid',2);
    % 
    %                 plot(unique_heading(1:zero_index)+0.2,tuning_mean_dora(2,1:zero_index,tuning_time_phase(pp),k),'^-','markersize',9,'color',colors(k,:),'markerfacecolor','none','LineWid',2);
    %                 plot(unique_heading(zero_index:end)-0.2,tuning_mean_dora(1,zero_index:end,tuning_time_phase(pp),k),'v-','markersize',9,'color',colors(k,:),'markerfacecolor','none','LineWid',2);

                        h = scatter(unique_heading(hh_index{hh})+0.2*sign(1.5-cc),...
                            tuning_mean_dora(cc,hh_index{hh},tuning_time_phase(pp),k),...
                            tuning_n_dora(cc,hh_index{hh},pp,k),...
                            colors(k,:),[cc_marker{cc}]);
                        
                        if cc == hh
                            set(h,'markerfacecolor',colors(k,:));
                        else
                            set(h,'markerfacecolor','none');
                        end
                    end
                    
                    plot(unique_heading,tuning_mean_dora(cc,:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2)
                end
                
                % plot(unique_heading,tuning_mean_dora(:,:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
                % h = errorbar([unique_heading'+0.2; unique_heading'-0.2],tuning_mean_dora(:,:,tuning_time_phase(pp),k),tuning_sem_dora(:,:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2,'linestyle','none');

                title([tuning_time_phase_title{pp} ', n = ' num2str(size(this_tuning_all,2))]);
                axis tight; xlim(xlim*1.1);
                
                % --------- Tuning all ---------
                subplot(3,length(tuning_time_phase),pp + length(tuning_time_phase));  hold on; ylabel('Correct + wrong');
                
                errorbar(unique_heading,tuning_mean_all(:,tuning_time_phase(pp),k),...
                    tuning_sem_all(:,tuning_time_phase(pp),k),'o-','color',colors(k,:),'LineWid',2,'capsize',0);
                
                
                % --------- Correct only with fitting ----------
                subplot(3,length(tuning_time_phase),pp + 2*length(tuning_time_phase));  hold on; ylabel('Correct only');
                
                plot(unique_heading_two_zeros,tuning_mean_correctonly(:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
                
                errorbar(unique_heading_two_zeros,tuning_mean_correctonly(:,tuning_time_phase(pp),k),...
                    tuning_sem_correctonly(:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2,'linestyle','none','capsize',0);
                
                % try errorbar_tick(h,10000); catch end;
                
                title([tuning_time_phase_title{pp} ', t = ' num2str(CP_ts{j}(tuning_time_phase(pp)))]);
                text(0,0.3+0.02*k,sprintf('n = %g, p = %.2g', size(tuning_cell_in_any_cell{k},1), ...
                           tuning_linear_fit_correctonly(2,tuning_time_phase(pp),k,1)),'color',colors(k,:));
                
                if linear_or_sigmoid == 2        % Plot sigmoid fitting
                    xx = linspace(min(unique_heading),max(unique_heading),100);
                    plot(xx,sigfunc(tuning_sig_fit_correctonly(:,tuning_time_phase(pp),k),xx),'color',colors(k,:),'linew',3);
               
                elseif linear_or_sigmoid == 1   % Plot linear fitting
                
                    line_type = {'-','--'}; % Significant/Insignificant
                    for hh = LEFT:RIGHT % LEFT, RIGHT
                        xx = unique_heading(hh_index{hh})';
                        para = tuning_linear_fit_correctonly(:,tuning_time_phase(pp),k,hh); % [r,p,b,k]
                        plot(xx,para(4)+para(3)*xx,'-','color',colors(k,:),'linew',3,'linestyle',line_type{2-(para(2)<0.05)});
                    end
                end
            end
            
            axis tight; xlim(xlim*1.1);
            SetFigure(15);

        end
        
        % Linear fitting parameters
        set(figure(145525),'name',['Linear fittings, j = ' num2str(j)]); clf;  hold on
        set(gcf,'uni','norm','pos',[0.002       0.532       0.411        0.38]);
        
        for k = 1:3
            % for hh = LEFT:RIGHT
                hh = 2;
                plot(CP_ts{j}, tuning_linear_fit_correctonly(3,:,k,hh),'linew',3,'linestyle',line_type{3-hh},'color',colors(k,:));
                sig = tuning_linear_fit_correctonly(2,:,k,hh) < 0.05;
                plot(CP_ts{j}(sig), mean(tuning_linear_fit_correctonly(3,sig,k,hh),3),'o','linew',2,...
                    'markersize',7,'color',colors(k,:));
                plot(CP_ts{j}, mean(tuning_linear_fit_correctonly(2,:,k,hh),3),[line_type{3-hh}],...
                    'linew',1,'color',colors(k,:));
            % end
        end
        ylabel('Slope of linear fitting');
        ylim([-0.03 0.05])
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j}(1), 0  + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        SetFigure(15);

    end


    dora_tuning_mean_each_cell = []; % Cache the dora tuning for each cell calculated in f1p2p8
    dora_tuning_sem_each_cell = []; % Cache the dora tuning for each cell calculated in f1p2p8
    dora_tuning_n_each_cell = [];
    partial_corr_timewins = {[0 1500],'0 ~ 1500 ms (all stim)';
                             % [-200 0],'-300 ~ 0 ms (before stim)';
                             [0 300],'0 ~ 300 ms (stim start)';
                             [400 700],'400 ~ 700 ms (acc peak)';
                             [700 1000],'700 ~ 1000 ms (intermediate)';
                             [1000 1300],'800 ~ 1100 ms (vel peak)';
                             [1300 1500],'1300 ~ 1500 ms (stim end)';
                             [1500 1800],'1500 ~ 2000 ms (delay)';
                             [2000 2300],'2000 ~ 2300 ms (postsac)'};
    select_for_partial = []; 

    function f1p2p8(debug)    % Rate 2.8. Partial correlation analysis. (Dora's method) 20170327 @UNIGE
        if debug
            dbstack;
            keyboard;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
        select_for_partial = select_tcells;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        find_for_partial = find(select_for_partial);
        unique_heading = group_result(representative_cell).unique_heading;
        partial_corr_coef_all = nan(sum(select_for_partial),size(partial_corr_timewins,1),2,3);  % [cell number, time epoch, (heading coeff, choice coeff), stim_type]
        partial_corr_p_all = nan(sum(select_for_partial),size(partial_corr_timewins,1),2,3);
        partial_corr_coef_all_flipped = nan(sum(select_for_partial),size(partial_corr_timewins,1),2,3);  % For plotting partial corr over time
        
        anova2_p_all = nan(sum(select_for_partial),size(partial_corr_timewins,1),2,3);
        dora_tuning_mean_each_cell = nan(sum(select_for_partial),size(partial_corr_timewins,1),3,2,length(unique_heading));
        dora_tuning_sem_each_cell = nan(sum(select_for_partial),size(partial_corr_timewins,1),3,2,length(unique_heading));
        dora_tuning_n_each_cell = zeros(sum(select_for_partial),size(partial_corr_timewins,1),3,2,length(unique_heading));
        
        % Calculate partial correlation coefficients and p-values for each cell
        progressbar('cell num');
        for i = 1:sum(select_for_partial)  % For cells
            
            this_raw_spike_in_bin = group_result(find_for_partial(i)).mat_raw_PSTH.spike_aligned{1,j};
            this_time = group_result(find_for_partial(i)).mat_raw_PSTH.spike_aligned{2,j};
            this_stim_type = group_result(find_for_partial(i)).mat_raw_PSTH.stim_type_per_trial;
            this_heading = group_result(find_for_partial(i)).mat_raw_PSTH.heading_per_trial;
            this_choice = group_result(find_for_partial(i)).mat_raw_PSTH.choice_per_trial;
            
            for tt = 1:size(partial_corr_timewins,1)
                count_win = partial_corr_timewins{tt,1}(1) <= this_time & this_time <= partial_corr_timewins{tt,1}(2);
                for k = 1:3
                    if isempty(find(this_stim_type==k, 1)); continue; end
                    
                    % --- Partial correlation --- 
                    X=[];
                    X(:,1) = sum(this_raw_spike_in_bin(this_stim_type==k,count_win),2)...
                        /range(partial_corr_timewins{tt,1})*1e3; % Average firing rate in Hz
                    X(:,2) = this_heading(this_stim_type==k);
                    X(:,3) = this_choice(this_stim_type==k);
                    
                    [r,p] = partialcorr(X);
                    partial_corr_coef_all(i,tt,:,k) = r(1,2:3);
                    partial_corr_p_all(i,tt,:,k) = p(1,2:3);
                    
                    % --- Anova2 p values (does not exclude interactions of heading and choice, but Dora used this) ---
                    anova2_p_all(i,tt,:,k) = anovan(X(:,1),{X(:,2),X(:,3)},'display','off')';
                    
                    % --- Dora tuning (putting it here is more flexbile than that in CP_HH calculation) ---
                    real_unique_this_heading = unique(this_heading);
                    withzero_unique_this_heading = unique([this_heading 0]); 
                    for hh = 1:length(real_unique_this_heading)
                        for cc = LEFT:RIGHT
                            this_select = (X(:,2) == real_unique_this_heading(hh))&(X(:,3)==cc);
                            this_mean = mean(X(this_select,1));
                            this_sem = std(X(this_select,1))/sqrt(sum(this_select));
                            hh_shouldbe = withzero_unique_this_heading==real_unique_this_heading(hh); % Deal with cases without zero heading 
                            dora_tuning_mean_each_cell(i,tt,k,cc,hh_shouldbe) = this_mean;
                            dora_tuning_sem_each_cell(i,tt,k,cc,hh_shouldbe) = this_sem;
                            dora_tuning_n_each_cell(i,tt,k,cc,hh_shouldbe) = sum(this_select);
                        end
                    end
                end
                
%                 % Use flipped partial corr (Fig.6 of Zaidel 2017)  
%                 if partial_corr_coef_all(i,tt,2,3) < 0  % Flip according to combined choice partial corr
%                     partial_corr_coef_all_flipped(i,tt,:,:) = - partial_corr_coef_all(i,tt,:,:);
%                 else
%                     partial_corr_coef_all_flipped(i,tt,:,:) = partial_corr_coef_all(i,tt,:,:);
%                 end
                
            end
            progressbar(i/sum(select_for_partial));
        end
        
        % Use R^2 (Fig.4 of Zaidel 2017)
        partial_corr_coef_all_flipped = partial_corr_coef_all.^2;
        
        
        %% Drawing
        set(figure(3099+figN),'name',sprintf('Partial correlation, j = %g, "any significant" out of N = %g, "%s" cells',...
                    j,sum(select_for_partial),t_criterion_txt)); clf; figN = figN+1;
        set(gcf,'uni','norm','pos',[0       0.038       0.994       0.877]);
        [h_sub,~] = tight_subplot(3,size(partial_corr_timewins,1),[0.05 0.02]);

        for tt = 1:size(partial_corr_timewins,1)
            for k = 1:3
                
                axes(h_sub(k+(tt-1)*3));
                % Use partial corr p value
%                 heading_sig = partial_corr_p_all(:,tt,1,k)<0.05;
%                 choice_sig = partial_corr_p_all(:,tt,2,k)<0.05;
                
                % Use ANOVA2 p value (Dora)
                heading_sig = anova2_p_all(:,tt,1,k)<0.05;
                choice_sig = anova2_p_all(:,tt,2,k)<0.05;
                
                h1=plot(partial_corr_coef_all((heading_sig|choice_sig) & ~(heading_sig & choice_sig),tt,1,k),...
                    partial_corr_coef_all((heading_sig|choice_sig) & ~(heading_sig & choice_sig),tt,2,k),...  % "One sig" cells
                    'o','color',colors(k,:));
                hold on;
                h2=plot(partial_corr_coef_all(heading_sig&choice_sig,tt,1,k),...
                    partial_corr_coef_all(heading_sig&choice_sig,tt,2,k),...  % "All sig" cells
                    'o','color','k','markerfacecol','k');
                
                xx = partial_corr_coef_all(heading_sig|choice_sig,tt,1,k);
                yy = partial_corr_coef_all(heading_sig|choice_sig,tt,2,k);
                h_all = plot(xx,yy,'visible','off');
                
                % Draw line if significant
                [rrr,ppp] = corr(xx,yy,'type','Pearson');
                if ppp < inf % 0.05
                    coeff = pca([xx yy]);
                    linPara(1) = coeff(2) / coeff(1);
                    linPara(2) = mean(yy)- linPara(1) *mean(xx);
                    
                    % -- Plotting
                    xxx = linspace(min(xx),max(xx),150);
                    yyy = linPara(1) * xxx + linPara(2);
                    plot(xxx,yyy,'linew',2,'color',colors(k,:));
                end
                
                select_actual_plot = zeros(N,1);
                find_for_partial = find(select_for_partial);
                select_actual_plot(find_for_partial(heading_sig|choice_sig)) = 1;
                set([gca h_all],'ButtonDownFcn',{@Show_individual_cell, h_all, select_actual_plot});

                
                axis([-1 1 -1 1]);
                hold on; plot([-1 1],[0 0],'k--'); plot([0 0],[-1 1],'k--');
                text(-0.9,-0.8,sprintf('%g,%g,%g\nr = %.3g\np= %.3g',sum(heading_sig),sum(choice_sig),sum(heading_sig&choice_sig),rrr,ppp));
                
                if k == 1
                    title(partial_corr_timewins{tt,2});
                end
                if tt == 1 && k == 1
                    xlabel('Partial corr (Heading)');
                    ylabel('Partial corr (Choice)');
                end
                
                                
            end
        end
        
       %% Plot partial correlation over time (Zaidel 2017, Figure 6)  HH20180821
        set(figure(194959),'name',sprintf('Partial correlation over time, j = %g, "any significant" out of N = %g, "%s" cells',...
                    j,sum(select_for_partial),t_criterion_txt)); clf; figN = figN+1;
        set(gcf,'uni','norm','pos',[0       0.038       0.994       0.877]);
        
        time_centers = mean(reshape([partial_corr_timewins{2:end,1}],2,[]),1);
        
        for k = 1:3
            subplot(1,3,k);
            for pp = 1:2
                
               aver_this = mean(partial_corr_coef_all_flipped(:,2:end,pp,k),1); % Note the first time epoch is the whole trial
               sem_this = std(partial_corr_coef_all_flipped(:,2:end,pp,k),[],1)/sqrt(size(partial_corr_coef_all_flipped,1));
               
               if pp == 1
                   col = colors(k,:);
               else
                   col = [0 0 0];
               end
               errorbar(time_centers, aver_this, sem_this,'color',col,'linew',2); hold on;
                
            end
        end
        title('Temporal resolution is too low here for Fisher info. I put the real one in f6p5r2.');
        
    end
        
    function f1p3(debug)      % Rate 3. Correct / Wrong Trials
        if debug
            dbstack;
            keyboard;
        end
        %%
        methods_of_select = {
            %select_bottom_line, 'All cells';
            select_tcells, 'Typical cells';
           % select_no_tcells, 'Non-typical cells'
           };
        
        unique_heading = group_result(representative_cell).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
        unique_heading_for_correct_wrong = unique_heading(unique_heading > 0);
        
        for ms = 1:size(methods_of_select,1)
            set(figure(2254-ms),'name',['Average PSTH (Correct vs Wrong), ' methods_of_select{ms,2}],'pos',[18 67 1645 898]); clf
            h_subplot = tight_subplot(3,1 + length(unique_heading_for_correct_wrong),[0.05 0.02],[0.1 0.1]);
            
            j = 1;
            for k = 1:3
                
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
                    correct_pref_ind = 2 + hh*2-1 ; % Index in PSTH_correct_angles_norm (exclude first two 0 heading)
                    correct_null_ind = 2 + hh*2 ;  % Index in PSTH_correct_angles_norm
                    wrong_pref_ind =  hh*2-1 ; % Index in PSTH_wrong_angles_norm
                    wrong_null_ind =  hh*2 ; % Index in PSTH_wrong_angles_norm
                    
                    for jjj = 1:2
                        PSTH_correct_wrong_this{jjj} = ...
                            cat(3, PSTH_correct_angles_Norm{jjj}(methods_of_select{ms,1},:,[correct_pref_ind correct_null_ind],k),...
                            PSTH_wrong_angles_Norm{jjj}(methods_of_select{ms,1},:,[wrong_pref_ind wrong_null_ind],k));
                    end
                    
                    SeriesComparison({PSTH_correct_wrong_this{1} PSTH_correct_wrong_this{2}},...
                        {rate_ts{1} rate_ts{2} time_markers},...r
                        'Colors',{'k','k','m','m'},'LineStyles',{'-','--'},...
                        'ErrorBar',0,'Xlabel',[],'Ylabel',[],'axes',h_subplot(hh * 3 + k));
                    
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

    function f10p1(debug)   % VarCE and CorCE (Churchland, Neuron, 2011).  HH20161101
        if debug; dbstack; keyboard; end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        select_for_VarCE_find =  find(select_tcells);
        VarCE_win_size = 120; % 60 ms   (Churchland, Neuron, 2011)
        VarCE_step_size = VarCE_win_size/4; % 10 ms
        VarCE_time_periods = {[-300,1800],[-900 100]};
        min_n_for_each_group = 3; % Must have at least min_M_for_each_group trials for each group to ensure the valid Var and Mean of the group
        n_bootstrap = 200;
        down_sample = 1; % As Churchland (4*15-ms bin = 60 ms, non-overlapping windows)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Time centers
        for j = 1:2
            t_centers{j} = [fliplr(VarCE_win_size/2-VarCE_step_size : -VarCE_step_size : min(VarCE_time_periods{j})+VarCE_win_size/2), ...
                VarCE_win_size/2 : VarCE_step_size: max(VarCE_time_periods{j})-VarCE_win_size/2]; % To ensure the first bin >0 aligned with [0,win_size]
        end
        t_centers{3} = t_centers{2}; % Same ts, different choices
        
        % I choose to calculate VarCE for each cell and then average them together instead of *ONLY* using the "residuals"
        % because LIP cells may have heterogeneous VarCE traces. However, when I combine all headings for one cell, I will
        % still use the "residual" method.
        % Notice that stim_type (k) and time alignment (j) would never be combined.
        
        % Preallocation
        residual_per_group = cell(length(select_for_VarCE_find),3,3); % residual_matrix{cell,k,(j=1;j=2 and Tin;j=2 and Tout)}
        n_per_group = residual_per_group;
        mean_per_group = residual_per_group;
        raw_per_group = mean_per_group;
        
        progressbar(['VarCE (n =' num2str(length(select_for_VarCE_find)) ')']);
        
        for nn = length(select_for_VarCE_find):-1:1  % Speed-up without preallocation because I'm lazy :)
            
            % Conditions for grouping
            stim_type_per_trial = group_result(select_for_VarCE_find(nn)).mat_raw_PSTH.stim_type_per_trial';
            heading_per_trial = group_result(select_for_VarCE_find(nn)).mat_raw_PSTH.heading_per_trial';
            unique_heading = unique(heading_per_trial);
            choice_per_trial = group_result(select_for_VarCE_find(nn)).mat_raw_PSTH.choice_per_trial;
            correct_or_zero_per_trial = (choice_per_trial == ((heading_per_trial>0) + 1))  |(heading_per_trial == 0); % Correct or 0 heading;
            pref_null = [group_result(select_for_VarCE_find(nn)).PREF_PSTH LEFT+RIGHT-group_result(select_for_VarCE_find(nn)).PREF_PSTH]; % Pref goes first
            
            % === 1. Calculate the Residual matrix ====
            
            for j = 1:2 % Two epochs
                
                % Get data
                spike_aligned_in_bin = group_result(select_for_VarCE_find(nn)).mat_raw_PSTH.spike_aligned{1,j};
                t_in_bin = group_result(select_for_VarCE_find(nn)).mat_raw_PSTH.spike_aligned{2,j};
                
                % Spike count using sliding window (different from what has been done in HeadingDis_cum_PSTH_HH)
                spike_hist = zeros(length(heading_per_trial),length(t_centers{j})); % Preallocation
                
                for ttt = 1:length(t_centers{j})
                    winBeg = find(t_in_bin>=t_centers{j}(ttt) - VarCE_win_size/2,1);
                    winEnd = find(t_in_bin>=t_centers{j}(ttt) + VarCE_win_size/2,1)-1;
                    spike_hist(:,ttt) = sum(spike_aligned_in_bin(: , winBeg:winEnd),2);  % Just spike counts
                end
                
                % Form groups and calculate Residuals for each group
                for k = 3:-1:1 % Three stim_types that would be never combined
                    select_this_stim_type = (stim_type_per_trial == k);
                    
                    if j == 1   % Aligned to stim on: group according to "Heading"
                        for hh = length(unique_heading):-1:1
                            select_this_group = select_this_stim_type & (heading_per_trial == unique_heading(hh));
                            n_per_group{nn,k,1}(hh,1) = sum(select_this_group); % Number of trials
                            
                            mean_this_group = mean(spike_hist(select_this_group,:),1);
                            residual_this_group = spike_hist(select_this_group,:) - ...
                                           repmat(mean_this_group,sum(select_this_group),1);  % Trials in this group x times
                            % fano_this_group = var(residual_this_group,0,1)./mean_this_group; % Fano = Measured variance / Mean
                                       
                            % - Data saving -
                            mean_per_group{nn,k,1}{hh,1} = repmat(mean_this_group,n_per_group{nn,k,1}(hh,1),1); % M_i x times, corresponding to each residual
                            residual_per_group{nn,k,1}{hh,1} = residual_this_group; % M_i x times
                            raw_per_group{nn,k,1}{hh,1} = spike_hist(select_this_group,:); % M_i x times
                            % fano_each_group{nn,k,j}{hh} = fano_this_group; % 1 x times 
                        end
                    elseif j == 2 % Aligned to saccade on: group according to "Heading x Choice"
                        for hh = length(unique_heading):-1:1
                            for cc = 1:2 % In (pref) & Out (null)
                                select_this_group = select_this_stim_type & (heading_per_trial == unique_heading(hh)) ...
                                    & (choice_per_trial == pref_null(cc)) & correct_or_zero_per_trial; % Pref goes first, and correct only (except 0)!
                                n_per_group{nn,k,1+cc}(hh,1) = sum(select_this_group); % Number of trials
                                
                                if sum(select_this_group) >= min_n_for_each_group % To avoid invalid Var and Means, only needed for j==2
                                    mean_this_group = mean(spike_hist(select_this_group,:),1);
                                    residual_this_group = spike_hist(select_this_group,:) - ...
                                        repmat(mean_this_group,sum(select_this_group),1);  % Trials in this group x times
                                    % fano_this_group = var(residual_this_group,0,1)./mean_this_group; % Fano = Measured variance / Mean

                                    % - Data saving
                                    mean_per_group{nn,k,1+cc}{hh,1} = repmat(mean_this_group,n_per_group{nn,k,1+cc}(hh,1),1); % Corresponding to each residual
                                    residual_per_group{nn,k,1+cc}{hh,1} = residual_this_group; % M_i x times, just for bootstraping
                                    raw_per_group{nn,k,1+cc}{hh,1} = spike_hist(select_this_group,:); % M_i x times
                                    % fano_each_group{nn,k,j}{hh,cc} = fano_this_group; % 1 x times
                                else % Too few trials (if no trials at all, it would be empty)
                                    mean_per_group{nn,k,1+cc}{hh,1} = nan(n_per_group{nn,k,1+cc}(hh,1),length(t_centers{j})); % Nans
                                    residual_per_group{nn,k,1+cc}{hh,1} = nan(n_per_group{nn,k,1+cc}(hh,1),length(t_centers{j})); % Nans
                                    raw_per_group{nn,k,1+cc}{hh,1} = nan(n_per_group{nn,k,1+cc}(hh,1),length(t_centers{j})); % Nans
                                end
                                
                            end
                        end
                    end % of "if j == 1 or 2"
                end % of stim_types
            end % of js
            
            % === 2. Find the phi for this cell ====
            % Note: I use the "grandFano" to find the proper phi instead of using Fano for each group, because 
            % if the trial num is too low for a certain group (the extreme case is 1), the Fano curves is too noisy.
            min_fano = inf;
            for k = 1:3
                for jandcc = 1:3 % To be more concise (when jandcc = 1, j=1; when jandcc >=2, j = 2, cc = jandcc-1)
                    zu_this = cat(1,residual_per_group{nn,k,jandcc}{:}); % Catenate all headings of this cell
                    mean_this = cat(1,mean_per_group{nn,k,jandcc}{:}); % The corresponding means of each trial
                    raw_this = cat(1,raw_per_group{nn,k,jandcc}{:}); % Cache for bootstrapping
                    
                    fano_this = nanvar(zu_this,0,1)./nanmean(mean_this,1);  % Equation(7) in Chruchland paper
                    min_fano = min(min_fano,min(fano_this));

                    fano_grand_per_cell{nn,k,jandcc} = fano_this;
                    zu_per_cell{nn,k,jandcc} = zu_this; % Cache for varCE_grand_all
                    mean_per_cell{nn,k,jandcc} = mean_this;
                    raw_per_cell{nn,k,jandcc} = raw_this;
                end
            end
            phis(nn) = min_fano; % Let phi be the minimal fano factor over epochs and conditions.
            % phi_per_cell(nn) = 0.35; % Hard-coded phis

            % === 3. Calculate VarCE for this cell ====
            for k = 1:3
                for jandcc = 1:3 % To be more concise (when jandcc = 1, j=1; when jandcc >=2, j = 2, cc = jandcc-1)
                    varCE_this_cell = nanvar(zu_per_cell{nn,k,jandcc},0,1)...
                        - phis(nn)*nanmean(mean_per_cell{nn,k,jandcc},1);  % Equation(6) in Chruchland paper
                    varCE_grand_per_cell{nn,k,jandcc} = varCE_this_cell;
                    phi_per_cell{nn,k,jandcc} = repmat(phis(nn),size(zu_per_cell{nn,k,jandcc},1),1); % Just for grand processing
                end                
            end
            
            progressbar((length(select_for_VarCE_find)-nn+1)/length(select_for_VarCE_find));
        end % of cells
        
        % === Plot results using average across cells ===
        % -- Plot phi distributions --
        figure(1839); clf; hold on
        hist(phis,20);
        title(sprintf('median = %g, IQR = %g, %g',median(phis),prctile(phis,25),prctile(phis,75)));
        
        figure(1840); clf; subplot 121; hold on;
        for k = 1:3
            plot(t_centers{1},mean(cat(1,varCE_grand_per_cell{:,k,1}),1),'color',colors(k,:),'linew',2.5);
        end
        
        subplot 122; hold on;
        for k = 1:3
            plot(t_centers{2},mean(cat(1,varCE_grand_per_cell{:,k,2}),1),'color',colors(k,:),'linestyle','-','linew',2.5);
            plot(t_centers{2},mean(cat(1,varCE_grand_per_cell{:,k,3}),1),'color',colors(k,:),'linestyle',':','linew',2.5);
        end

        % === 4. Combined Zu across cells instead of averaging ===
        for k = 1:3
            for jandcc = 1:3  % (when jandcc = 1, j=1; when jandcc >=2, j = 2, cc = jandcc-1)
                zu_grand_all{k,jandcc} = cat(1,zu_per_cell{:,k,jandcc});
                mean_grand_all{k,jandcc} = cat(1,mean_per_cell{:,k,jandcc});
                raw_grand_all{k,jandcc} = cat(1,raw_per_cell{:,k,jandcc});
                phi_grand_all{k,jandcc} = cat(1,phi_per_cell{:,k,jandcc});
                
                % Knock out nans
                select_invalid = any(isnan(zu_grand_all{k,jandcc}),2);
                fprintf('Number of nan trials = %g\n', sum(select_invalid));
                zu_grand_all{k,jandcc}(select_invalid,:) = [];
                mean_grand_all{k,jandcc}(select_invalid,:) = [];
                phi_grand_all{k,jandcc}(select_invalid,:) = [];
                raw_grand_all{k,jandcc}(select_invalid,:) = [];
                
                n_all{k,jandcc} = cat(1,n_per_group{:,k,jandcc});
                n_all{k,jandcc}(n_all{k,jandcc} < min_n_for_each_group) = [];
                
                varCE_grand_all{k,jandcc} = var(zu_grand_all{k,jandcc},0,1)...
                        - mean(mean_grand_all{k,jandcc} .* repmat(phi_grand_all{k,jandcc},1,size(mean_grand_all{k,jandcc},2)),1);  % Equation(6) in Chruchland paper
                n_trials_grand_all(k,jandcc) =  sum(all(~isnan(mean_grand_all{k,jandcc}),2));
            end
        end
        
        figure(1841); clf; subplot 121; hold on;
        for k = 1:3
            plot(t_centers{1},varCE_grand_all{k,1},'color',colors(k,:),'linew',2.5);
        end
        title(sprintf('%g, ',n_trials_grand_all(:,1)));
        
        subplot 122; hold on;
        for k = 1:3
            plot(t_centers{2},varCE_grand_all{k,2},'color',colors(k,:),'linestyle','-','linew',2.5);
            plot(t_centers{2},varCE_grand_all{k,3},'color',colors(k,:),'linestyle','--','linew',2.5);
        end
        title(sprintf('%g, ',n_trials_grand_all(:,2:3)));
        
        % === 4.5 Calculate CorCE === @HH20161220
        set(figure(1842),'unit','norm','pos',[0.0142857142857143 0.0904761904761905 0.857142857142857 0.38952380952381]); clf;
        down_sample_cov = round(VarCE_win_size/VarCE_step_size);
        t_s = t_centers{1}(1:down_sample_cov:end);
        
        for k = 1:3
            % --- (1) Total variance
            cov_grand_all = cov(zu_grand_all{k,1}(:,1:down_sample_cov:end));
            
            % --- (2) Replace the diagonal of cov with VarCEs
            cov_grand_all = cov_grand_all - diag(diag(cov_grand_all)) + diag(varCE_grand_all{k,1}(1:down_sample_cov:end));
            
            % --- (3) Solve corCE (r_ij)
            corCE_grand_all{k,1} = cov_grand_all./sqrt(diag(cov_grand_all)*diag(cov_grand_all)'); % Equation 9 of Ann, 2012
            
            % Ploting
            subplot(1,3,k); imagesc(t_s,t_s,corCE_grand_all{k,1}); axis square;colorbar;
            map = load('CorCE_colormap.mat');
            colormap(map.map);
        end
        drawnow;
        
        % === 5. Bootstrap to get the SE of varCE ===
        if ver_num < 2014
            if matlabpool('size') == 0 ;  matlabpool;  end
        else
            if isempty(gcp('nocreate')); parpool; end
        end

        for k = 1:3
            for jandcc = 1:3
                
                % Speed up the parfor loop
                n_all_tmp = n_all{k,jandcc};
                phi_grand_all_tmp = phi_grand_all{k,jandcc};
                zu_grand_all_tmp = zu_grand_all{k,jandcc}(:,1:down_sample:end);
                raw_grand_all_tmp = raw_grand_all{k,jandcc}(:,1:down_sample:end);
                
                varCE_grand_boot_this = nan(n_bootstrap,size(zu_grand_all_tmp,2));
                parfor_progress(n_bootstrap);
                
                parfor bb = 1:n_bootstrap
                    % Randomly select trials while keeping the ns for each condition
                    zu_grand_boot_this = nan(size(zu_grand_all_tmp));
                    mean_grand_boot_this = zu_grand_boot_this;
                    pointer_first_this = 1; % Indicate the first line of the current group
                    
                    for mm = 1:length(n_all_tmp) % For all groups
                        n_this = n_all_tmp(mm);
                        rand_choose = fix(rand(1,n_this)*n_this) + 1; % Random choose with replacement                        
                        zu_grand_boot_this(pointer_first_this:pointer_first_this+n_this-1,:) = ...
                            zu_grand_all_tmp(pointer_first_this + rand_choose - 1,:); 
                        mean_grand_boot_this(pointer_first_this:pointer_first_this+n_this-1,:) = ...
                            repmat(mean(raw_grand_all_tmp(pointer_first_this:pointer_first_this+n_this-1,:),1),n_this,1);
                        
                        pointer_first_this = pointer_first_this + n_this; % Move to the first line of the next group
                    end
                    
                    % Calculate bootstrapped varCE (phis are not changing)
                    varCE_grand_boot_this(bb,:) = var(zu_grand_boot_this,0,1)...
                        - mean(mean_grand_boot_this .* repmat(phi_grand_all_tmp,1,size(mean_grand_boot_this,2)),1);  % Equation(6) in Chruchland paper
                    
                    parfor_progress; 
                end
                
                varCE_grand_boot_all{k,jandcc} = varCE_grand_boot_this;
            end
        end
        
        %{
        figure(1841); clf; subplot 121; hold on;
        for k = 1:3
            shadedErrorBar(t_centers{1},varCE_grand_boot_all{k,1},{@mean,@std},{'color',colors(k,:),'linew',2.5},0.4);
        end
        title(sprintf('%g, ',n_trials_grand_all(:,1)));
        
        subplot 122; hold on;
        for k = 1:3
            shadedErrorBar(t_centers{2},varCE_grand_boot_all{k,2},{@mean,@std},{'color',colors(k,:),'linew',2.5,'linestyle','-'},0.4);
            shadedErrorBar(t_centers{3},varCE_grand_boot_all{k,3},{@mean,@std},{'color',colors(k,:),'linew',2.5,'linestyle','--'},0.4);
        end
        title(sprintf('%g, ',n_trials_grand_all(:,2:3)));
        %}

        figure(1841); subplot 121; hold on;
        for k = 1:3
            h = errorbar(t_centers{1}(1:down_sample:end),varCE_grand_all{k,1}(1:down_sample:end),...
                std(varCE_grand_boot_all{k,1},[],1), 'color',colors(k,:),'linestyle','none','linew',2.5);
            errorbar_tick(h,inf);
        end
        
        subplot 122; hold on;
        for k = 1:3
            h = errorbar(t_centers{2}(1:down_sample:end),varCE_grand_all{k,2}(1:down_sample:end),...
                 std(varCE_grand_boot_all{k,2},[],1),'color',colors(k,:),'linestyle','none','linew',2.5);
            errorbar_tick(h,inf);
            h = errorbar(t_centers{3}(1:down_sample:end),varCE_grand_all{k,3}(1:down_sample:end),...
                 std(varCE_grand_boot_all{k,3},[],1),'color',colors(k,:),'linestyle','none','linew',2.5);
            errorbar_tick(h,inf);
        end
        
        
        SeriesComparison({reshape([varCE_grand_boot_all{:,1}],n_bootstrap,[],3) reshape([varCE_grand_boot_all{:,2}],n_bootstrap,[],3)},...
            {t_centers{1}(1:down_sample:length(t_centers{1})) t_centers{2}(1:down_sample:length(t_centers{2})) time_markers},...
            'ErrorBar',6,'SEM',0,'Border',[1600, -550],...
            'CompareIndex',[1 2 3;2 3 1],'CompareColor',{colors(1,:),colors(2,:),[0 0.8 0.4]},...
            'Colors',{colors(1,:),colors(2,:),[0 0.8 0.4]},'LineStyles',{'-'});


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
            
            for k = 1:3
                
                %     selectCells_notNaN =  select_for_CP & (~isnan(CP(:,1,k)));
                ys = mean(CP{j}(select_for_CP,:,k));
                errors = std(CP{j}(select_for_CP,:,k))/sqrt(sum(select_for_CP));
                h = shadedErrorBar(CP_ts{j},ys,errors,'lineprops',{'Color',colors(k,:)},'transparent',transparent);
                set(h.mainLine,'LineWidth',2);
                hold on;
                
            end
            
            xlabel(['Center of ' num2str(group_result(representative_cell).mat_raw_PSTH.binSize_CP) ' ms time window']);
            
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
            
            h = shadedErrorBar(rate_ts{j},ys_MD,err_MD,'lineprops',{'Color',[0.5 0.5 0.5]},'transparent',transparent);
            set(h.mainLine,'LineWidth',2);
            
            % Comb v.s. Visual
            ModDiv_to_average = ModDiv_All{j}(select_for_div,:,3);
            ModDiv_to_average = ModDiv_to_average.* repmat(sign(sum(ModDiv_to_average,2)),1,size(ModDiv_to_average,2));
            
            ys_MD = mean(ModDiv_to_average);
            err_MD = std(ModDiv_to_average)/sqrt(size(ModDiv_to_average,1));
            
            h = shadedErrorBar(rate_ts{j},ys_MD,err_MD,'lineprops',{'Color',[0.5 0.5 0.5]},'transparent',transparent);
            set(h.mainLine,'LineWidth',2,'linestyle','-.'); hold on;
            
            
            % Choice divergence
            for k = 1:3
                selectCells_notNaN =  select_for_div & (~isnan(ChoiceDiv_All{j}(:,1,k)));
                
                ys_CD{j}(k,:) = mean(ChoiceDiv_All{j}(selectCells_notNaN,:,k));
                errors = std(ChoiceDiv_All{j}(selectCells_notNaN,:,k))/sqrt(sum(selectCells_notNaN));
                h = shadedErrorBar(rate_ts{j},ys_CD{j}(k,:),errors,'lineprops',{'Color',colors(k,:),'linestyle','-'},'transparent',transparent);
                set(h.mainLine,'LineWidth',2);
            end
            
            % Linear sum of vest and vis
            plot(rate_ts{j}, sum(ys_CD{j}(1:2,:)),'k-');
            
            axis tight;  ylim([-0.1 0.9]);
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
            'Colors',{colors(1,:),colors(2,:),colors(3,:)},'LineStyles',{'-'},...
            'ErrorBar',6,'Xlabel',[],'Ylabel',[],'transparent',0,...
            'CompareIndex',[1:3;1:3],...
            'CompareColor',[mat2cell(colors,ones(3,1))]);
        
        xlim([-300 2300]); ylim([-0.1 0.7 ]); legend off;  plot(xlim,[0 0],'k--'); 
        title(sprintf('Choice divergence (SU, N = %g)',sum(select_for_div)));
        set(gcf,'unit','norm','pos',[0.00952380952380952 0.521904761904762 0.315476190476191 0.385714285714286]);
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
        
        SetFigure(17);
        
        %% A or V? HH20151220
        figure(2247); clf; plot(rate_ts{1},ys_CD{1}');  
        hold on; title('A or V?')
        plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',3,'color',[0.6 0.6 0.6]);
        dt_Gauss = Gauss_vel(2,1)-Gauss_vel(1,1);
        Gauss_acc = diff(Gauss_vel(:,2));
        plot(Gauss_vel(2:end,1),Gauss_acc);
        plot(Gauss_vel(1:end,1),[0 ;cumsum(abs(Gauss_acc))]/10*3,'k','linew',2);
        plot(Gauss_vel(1:end,1),[cumsum(Gauss_vel(:,2))]/100*2/3,'k','linew',2);
        
        
        %% Linear fitting combined trace
        ts = rate_ts{1};
        t_select = rate_ts{1}> 0 & rate_ts{1} <=1500;
   
        ramping1 = ys_CD{1}(1,t_select);
        ramping2 = ys_CD{1}(2,t_select);
        ramping3 = ys_CD{1}(3,t_select);
        
        w = fminsearch(@(w) sum((w(1)*ramping1 + w(2)*ramping2 - ramping3).^2), [.5 .5]);
        
        figure(2308); clf
        plot(ts(t_select),ramping1,'b',ts(t_select),ramping2,'r',...
            ts(t_select),ramping3,'g',ts(t_select),ramping1*w(1)+ramping2*w(2),'k');
        title(sprintf('Fitting Choice Div: w1 = %g, w2 = %g',w(1),w(2)))
        
        
    end
    function f2p2(debug)      % ROC 2. CDiv: Multisensory Enhancement
        if debug
            dbstack;
            keyboard;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         select_for_div =  select_bottom_line;
        select_for_div =  select_tcells;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% CDiv: Multisensory Enhancement
        % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) > 0) & (xls_num{1}(:,header.HD_TargFirst)~=0);% ...
        % & (~isnan(ChoiceDiv_All(:,1,k))) ;% & MemSac_DDI(:,4)<0.55 ;
        
        set(figure(2099+figN),'name','Multisensory Enhancement of CDiv'); clf; figN = figN+1;
        set(gcf,'uni','norm','pos',[0.044        0.49       0.604       0.289]);
        
        for k = 1:3
            h = subplot(1,3,k);
            
            select_this_k = select_for_div & (~isnan(ChoiceDiv_ModDiffer{j}(:,1,k)));
            
            SeriesComparison({ChoiceDiv_ModDiffer{1}(select_this_k,:,k) ChoiceDiv_ModDiffer{2}(select_this_k,:,k)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'ErrorBar',6,...
                'CompareIndex',[1;1],'CompareColor',{colors(k,:)},'PlotPs',1,...
                'Colors',{colors(k,:)},'YLim',[-0.1 0.3],...
                'Transparent',transparent,'LineStyles',{'-'},'axes',h);
            
            legend off; xlim([-400 2100]);
            plot(xlim,[0 0],'k--');
            
            % Gaussian vel
            plot(Gauss_vel(:,1) ,Gauss_vel(:,2)*range(ylim)/5 ,'--','linew',2,'color',[0.6 0.6 0.6]);
            
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
                
                ys = mean(ChoiceDiv_Easy{j}(select_this_k,:,k),1);
                errors = std(ChoiceDiv_Easy{j}(select_this_k,:,k),[],1)/sqrt(sum(select_this_k));
                h = shadedErrorBar(rate_ts{j},ys,errors,'lineprops',{'Color',colors(k,:)},'transparent',transparent);
                set(h.mainLine,'LineWidth',2);  hold on;
                
                ys = mean(ChoiceDiv_Difficult{j}(select_this_k,:,k),1);
                errors = std(ChoiceDiv_Difficult{j}(select_this_k,:,k),[],1)/sqrt(sum(select_this_k));
                h = shadedErrorBar(rate_ts{j},ys,errors,'lineprops',{'Color',colors(k,:)*0.4 + [0.6 0.6 0.6]},'transparent',transparent);
                set(h.mainLine,'LineWidth',2);  hold on;
                
                xlim([time_markers{j}(1,1)-200 time_markers{j}(1,3)+200]); ylim([-0.2 0.7]); 
                
                
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
        
        set(figure(2099+figN),'name','Easy-difficult of CDiv'); clf; figN = figN+1;
        set(gcf,'uni','norm','pos',[0.006       0.582       0.757        0.33]);
        %{
        for j = 1:2
            for k = 1:3
                subplot(2,3,k + (j-1)*3);
                
                select_this_k = select_for_div & (~isnan(ChoiceDiv_All{j}(:,1,k)));
                
                ys = mean(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k),1);
                errors = std(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k),[],1)/sqrt(sum(select_this_k));
                [~,ps] = ttest(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k));  % p_value
                
                h = shadedErrorBar(rate_ts{j},ys,errors,'lineprops',{'Color',colors(k,:)},'transparent',transparent);
                set(h.mainLine,'LineWidth',2); ; hold on;
                
                % xlim([time_markers{j}(1,1)-100 time_markers{j}(1,3)+200]); ylim([-0.2 0.5]); ylim([-0.1 0.2]);
                xlim([  -405        2095]); ylim([-0.05 0.2]);
                
                plot(rate_ts{j},ps,'k');
                
                if sum(ps < p_critical)>0
                    plot(rate_ts{j}(ps < p_critical),max(ylim)*.9,'.','color',modality_diff_colors(k,:));
                    disp('ps = ');
                    ps(ps<p_critical)
                end
                
                plot(xlim,[0 0],'k--');
                for tt = 1:3
                    plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
                end
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/5 + 0,'--','linew',2,'color',[0.6 0.6 0.6]);
                
            end
        end
        %}
        
        for k = 1:3
            h = subplot(1,3,k);
            
            select_this_k = select_for_div & (~isnan(ChoiceDiv_All{j}(:,1,k)));
            
            for j = 1:2
                ys_this{j,k} = mean(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k),1);
                errors_this{j} = std(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k),[],1)/sqrt(sum(select_this_k));
                [~,ps{j}] = ttest(ChoiceDiv_EasyMinusDifficult{j}(select_this_k,:,k));  % p_value
            end
            
            SeriesComparison({shiftdim(ys_this{1,k}',-1) shiftdim(ys_this{2,k}',-1)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'OverrideError',{errors_this{1}', errors_this{2}'},...
                'OverridePs',{ps{1}', ps{2}'},'ErrorBar',6,...
                'CompareIndex',[1;1],'CompareColor',{colors(k,:)},...
                'Colors',{colors(k,:)},'YLim',[-0.1 0.2],...
                'Transparent',transparent,'LineStyles',{'-'},'axes',h);
            
            % xlim([time_markers{j}(1,1)-100 time_markers{j}(1,3)+200]); ylim([-0.2 0.5]); ylim([-0.1 0.2]);
            
            hold on;
            % Gaussian vel
            plot(Gauss_vel(:,1), Gauss_vel(:,2)*range(ylim)/5 + 0,'--','linew',2,'color',[0.6 0.6 0.6]);

            plot(rate_ts{1},ps{1},'k');
            plot(rate_ts{2},ps{2},'c');
            
            plot(xlim,[0 0],'k--');
            
            text(0.1,0,sprintf('min Pvalue = %g',min(ps{1})));
            
            if k == 3 
                plot(rate_ts{1}, ys_this{1,1} + ys_this{1,2},'m--');
            end
            
            xlim([-405 2095]); ylim([-0.1 0.2]);
            legend off;
        end
        
        title(sprintf('Easy - Difficult (N = %g)',sum(select_this_k)));
        SetFigure(15);
        
        %{
        % Use f1p2p7 Dora tuning instead. HH20180612
        
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
%         t_pre_center = mean(time_markers{j}(1,3)) - group_result(representative_cell).mat_raw_PSTH.binSize_CP/2;   % Pre-sac epoch
%         t_post_center = mean(time_markers{j}(1,3)) + group_result(representative_cell).mat_raw_PSTH.binSize_CP/2;  % Post-sac epoch
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
        unique_heading = group_result(representative_cell).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
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
                for k = 1:3
                    this_raw = group_result(find_for_tuning(i)).mat_raw_PSTH.CP{j,k}.raw_CP_result{pp};
                    this_tuning = this_raw.Neu_tuning(:,2);
                    this_tuning_correctonly = this_raw.Neu_tuning_correctonly(:,2);
                    
                    % @HH20160214: Patch calculation for LEFT & RIGHT choices of 0 headings (because I failed to do this
                    % in the original CP_HH and I feel hesitant to redo all the batch files from A to Z right now...).
                    % Note, sadly, that this could be NaN because I also failed to pack the spike counts into
                    % spike_counts_allheadings_grouped if the number of left/right choices were fewer than 3... WTF...

                    % @HH20160907. Today I redo all the batch files finally!
                    
                    if length(this_tuning) == length(unique_heading)
                        zero_spikes_left_right = [mean(this_raw.spike_counts_allheadings_grouped{zero_index, 1}); mean(this_raw.spike_counts_allheadings_grouped{zero_index, 2})];
                        this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); zero_spikes_left_right; this_tuning_correctonly(zero_index+1:end)];
                    else % No zero heading at all in the file 
                        zero_spikes_left_right = [nan;nan];
                        this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); zero_spikes_left_right; this_tuning_correctonly(zero_index:end)];
                        this_tuning = [this_tuning(1:zero_index-1); nan ;this_tuning(zero_index:end)];
                    end
                    
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
                        keyboard;
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
                % errorbar_tick(h,10000); % Obsolete
                title([tuning_time_phase_title{pp} ', n = ' num2str(size(this_tuning_all,2))]);
                axis tight; xlim(xlim*1.1);
                
                subplot(2,length(tuning_time_phase),pp + length(tuning_time_phase));  hold on; ylabel('Correct only');
                plot(unique_heading_two_zeros,tuning_mean_correctonly(:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
                h = errorbar(unique_heading_two_zeros,tuning_mean_correctonly(:,tuning_time_phase(pp),k),tuning_sem_correctonly(:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2);
                % try errorbar_tick(h,10000); catch end;   % Obsolete
                
                title([tuning_time_phase_title{pp} ', t = ' num2str(CP_ts{j}(tuning_time_phase(pp))),...
                    ' n = ' num2str(size(this_tuning_correctonly,2))]);
                axis tight; xlim(xlim*1.1);
                
                SetFigure(15);
                
                % Plot sigmoid fitting
                xx = linspace(min(unique_heading),max(unique_heading),100);
                plot(xx,sigfunc(tuning_sig_fit_correctonly(:,tuning_time_phase(pp),k),xx),'color',colors(k,:),'linew',3);
                
            end
        end
        %}
        
        %{
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
                
%                 load('colormap_for_tuning_hotgram','colormap_for_tuning_hotgram');
%                 colormap(colormap_for_tuning_hotgram);
                
                %                 axis tight;
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
                
            end
            
            SetFigure();
        end
        %}
    end

    marker_size = 11;
    tcell_cross_size = 15;

    function f3p1(debug)      % Correlations 1. Mem-sac v.s. CDiv
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
            'Xlabel',MemSac_indicator_txt,'Ylabel','Pre-sac CDiv','FaceColors',{colors(1,:),colors(2,:),colors(3,:)},'Markers',{'o'},...
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
            'Xlabel',MemSac_indicator_txt,'Ylabel','Post-sac CDiv','FaceColors',{colors(1,:),colors(2,:),colors(3,:)},'Markers',{'o'},...
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
        %         'Xlabel','MemSac Rate Diff (pre)','Ylabel','Pre-sac Cdiv','FaceColors',{colors(1,:),colors(2,:),colors(3,:)},'Markers',{'o'},...
        %         'LineStyles',{'b-','r-','g-'},'MarkerSize',12,'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
        
        % =================== Correlations 2. Mem-sac v.s. CP
        
        j = 2; % Not much time-sensitive, so I use j = 2 here.
        CP_interest = [-400 -100];  % min and max interval. HH20180609
        
        %  CP_interval = CP_ts{j} == CP_interest ;
        
        CP_interval = find(abs(CP_ts{j} - CP_interest(1)) == min(abs(CP_ts{j} - CP_interest(1)))) : ...
                      find(abs(CP_ts{j} - CP_interest(2)) == min(abs(CP_ts{j} - CP_interest(2))));
        
        for k = 1:3
            CP_NS{k} = find_bottom_line(any(CP_p{j}(select_bottom_line,CP_interval,k) > 0.05,2));
            CP_S{k} = find_bottom_line(any(CP_p{j}(select_bottom_line,CP_interval,k) <= 0.05,2));
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
            'FaceColors',{'none',colors(1,:),'none',colors(2,:),'none',colors(3,:)},'Markers',{'o'},...
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
        
        % ========== CD and CP ========        
        h = LinearCorrelation({  mean(ChoiceDiv_All{j}(CP_NS{1},rate_ts{j} > -400 & rate_ts{j} < -100,1),2) / 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_S{1},rate_ts{j} > -400 & rate_ts{j} < -100,1),2) / 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_NS{2},rate_ts{j} > -400 & rate_ts{j} < -100,2),2)/ 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_S{2},rate_ts{j} > -400 & rate_ts{j} < -100,2),2)/ 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_NS{3},rate_ts{j} > -400 & rate_ts{j} < -100,3),2)/ 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_S{3},rate_ts{j} > -400 & rate_ts{j} < -100,3),2)/ 2 + 0.5},...
            {  mean(CP{j}(CP_NS{1},CP_interval,1),2),...
            mean(CP{j}(CP_S{1},CP_interval,1),2),...
            mean(CP{j}(CP_NS{2},CP_interval,2),2),...
            mean(CP{j}(CP_S{2},CP_interval,2),2),...
            mean(CP{j}(CP_NS{3},CP_interval,3),2),...
            mean(CP{j}(CP_S{3},CP_interval,3),2)},...
            'CombinedIndex',[3 12 48],...
            'Xlabel','Choice Div (pre)/2 + 0.5','Ylabel',['CP at Sac ' num2str(CP_interest)],...
            'FaceColors',{'none',colors(1,:),'none',colors(2,:),'none',colors(3,:)},'Markers',{'o'},...
            'LineStyles',{'b:','b:','r:','r:','g:','g:','b-','r-','g-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',nHist,'YHist',nHist,'SameScale',1); figN = figN + 1;
        delete([h.group(1:6).line]);
        plot(xlim,[0.5 0.5],'k--'); plot([0.5 0.5],ylim,'k--');         SetFigure(20);
        
        
        set(gcf,'name',['j = ' num2str(j)]);
        % Annotate tcells
        plot(.5 + .5* mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),mean(CP{j}(select_tcells,CP_interval,1),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(.5 + .5*mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),mean(CP{j}(select_tcells,CP_interval,2),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        plot(.5 + .5*mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,3),2),mean(CP{j}(select_tcells,CP_interval,3),2),...
            '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        
        % ------------ Memsac DDI and Choice Preference
        %
        %         k = 1;
        %         tt = 1;
        %         cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
        %
        %         MemSac_DDI_selected = MemSac_indicator(select_bottom_line);
        %
        %         monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        %         monkey1 = monkeys == 5; monkey1 = monkey1(select_bottom_line)';
        %         monkey2 = monkeys == 10; monkey2 = monkey2(select_bottom_line)';
        %
        %
        %         h = LinearCorrelation({
        %             MemSac_DDI_selected(monkey1 & (~cpref_sig));
        %             MemSac_DDI_selected(monkey1 &(cpref_sig));
        %             MemSac_DDI_selected(monkey1 &(~cpref_sig));
        %             MemSac_DDI_selected(monkey1 &(cpref_sig));
        %             MemSac_DDI_selected(monkey1 &(~cpref_sig));
        %             MemSac_DDI_selected(monkey1 &(cpref_sig));
        %             MemSac_DDI_selected(monkey2 &(~cpref_sig));
        %             MemSac_DDI_selected(monkey2 &(cpref_sig));
        %             MemSac_DDI_selected(monkey2 &(~cpref_sig));
        %             MemSac_DDI_selected(monkey2 &(cpref_sig));
        %             MemSac_DDI_selected(monkey2 &(~cpref_sig));
        %             MemSac_DDI_selected(monkey2 &(cpref_sig));
        %             },...
        %             {
        %             abs(Choice_pref_all(1,monkey1 & ~cpref_sig,tt)) ;
        %             abs(Choice_pref_all(1,monkey1 & cpref_sig,tt)) ;
        %             abs(Choice_pref_all(2,monkey1 & ~cpref_sig,tt)) ;
        %             abs(Choice_pref_all(2,monkey1 & cpref_sig,tt)) ;
        %             abs(Choice_pref_all(3,monkey1 & ~cpref_sig,tt)) ;
        %             abs(Choice_pref_all(3,monkey1 & cpref_sig,tt)) ;
        %             abs(Choice_pref_all(1,monkey2 & ~cpref_sig,tt)) ;
        %             abs(Choice_pref_all(1,monkey2 & cpref_sig,tt)) ;
        %             abs(Choice_pref_all(2,monkey2 & ~cpref_sig,tt)) ;
        %             abs(Choice_pref_all(2,monkey2 & cpref_sig,tt)) ;
        %             abs(Choice_pref_all(3,monkey2 & ~cpref_sig,tt)) ;
        %             abs(Choice_pref_all(3,monkey2 & cpref_sig,tt)) ;
        %             },...
        %             'CombinedIndex',[195 780 3120],...  % (1,2,7,8);(3,4,9,10);(5,6,11,12)
        %             'Xlabel',MemSac_indicator_txt,'Ylabel','abs(Choice preference)',...
        %             'FaceColors',{'none',colors(1,:),'none',colors(2,:),'none',colors(3,:)},'Markers',{'o','o','o','o','o','o','^','^','^','^','^','^',},...
        %             'LineStyles',{'b:','b:','r:','r:','g:','g:','b:','b:','r:','r:','g:','g:','b-','r-','g-'},'MarkerSize',marker_size,...
        %             'figN',figN,'XHist',20,'YHist',20,...
        %             'XHistStyle','stacked','YHistStyle','grouped','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        
    end

    function f3p1p2(debug)      % Mem-sac v.s. CPref
        if debug
            dbstack;
            keyboard;
        end
        % ------------ Memsac DDI and Choice Preference (single modality)
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 5; monkey1 = monkey1(select_bottom_line)';
        monkey2 = monkeys == 10; monkey2 = monkey2(select_bottom_line)';
        
        memsac_sig = MemSac_indicator_p' < 0.05;
        
        for k = 1:3
            
            tt = 3; % From stim on to stim off
            
            cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
            
            MemSac_DDI_selected = MemSac_indicator(select_bottom_line);
            
            h = LinearCorrelation({
                MemSac_DDI_selected(monkey1 & ~or(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey1 & xor(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey1 & and(cpref_sig, memsac_sig));
                
                MemSac_DDI_selected(monkey2 & ~or(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey2 & xor(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey2 & and(cpref_sig, memsac_sig));
                },...
                {
                abs(Choice_pref_all(k,monkey1 & ~or(cpref_sig, memsac_sig),tt)) ;
                abs(Choice_pref_all(k,monkey1 & xor(cpref_sig, memsac_sig),tt));
                abs(Choice_pref_all(k,monkey1 & and(cpref_sig, memsac_sig),tt));
                
                abs(Choice_pref_all(k,monkey2 & ~or(cpref_sig, memsac_sig),tt)) ;
                abs(Choice_pref_all(k,monkey2 & xor(cpref_sig, memsac_sig),tt)) ;
                abs(Choice_pref_all(k,monkey2 & and(cpref_sig, memsac_sig),tt));
                
                },...
                'CombinedIndex',[48 63],...
                'Xlabel',MemSac_indicator_txt,'Ylabel','abs(Choice preference)',...
                'FaceColors',{'none',colors(k,:)*0.2 + [0.8 0.8 0.8],colors(k,:)},'Markers',{'o','o','o','^','^','^'},...
                'LineStyles',{'k:','k:','k:','k:','k:','k:','k:','k-'},'MarkerSize',15,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,...
                'Method','Pearson','FittingMethod',2); figN = figN + 1;
                        

            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
            delete([h.group(1:4).line]);
            % Annotate tcells
            plot(MemSac_DDI_selected(select_tcells(select_bottom_line)),abs(Choice_pref_all(k,select_tcells(select_bottom_line),tt)),...
                '+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(MemSac_DDI_selected(:),abs(Choice_pref_all(k,:,tt)),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_bottom_line});
            
            % Override histogram
            cla(h.ax_xhist);
            HistComparison( {MemSac_DDI_selected(memsac_sig), MemSac_DDI_selected(~memsac_sig)},...
                'EdgeColors',{'k','k'}, 'FaceColors',{'k', 'none'},...
                'MeanType','Median','Axes',h.ax_xhist,'TTest',0);
            h.ax_xhist.XLim = xlim(h.ax_raw);
            
            cla(h.ax_yhist);
            HistComparison( {abs(Choice_pref_all(k,cpref_sig,tt))',abs(Choice_pref_all(k, ~cpref_sig,tt))'},...
                           'EdgeColors',{'k','k'}, 'FaceColors',{colors(k,:), 'none'},...
                           'MeanType','Median','Axes',h.ax_yhist,'TTest',0);
            h.ax_yhist.XLim = ylim(h.ax_raw);
            
        end
    end
    function f3p2(debug)      % Correlations 2. Choice preference v.s. modality preference (Anne Fig.2f-h; Fig.3a)
        if debug
            dbstack;
            keyboard;
        end
        
        %% ---------------- 1 Choice pref v.s. modality pref (Anne Fig.3a)-------------------
        tt = 3;
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 5; monkey1 = monkey1(select_bottom_line)';
        monkey2 = monkeys == 10; monkey2 = monkey2(select_bottom_line)';
        
        for k = 1:3
            set(figure(figN),'name','CDiv vs. MDiv','pos',[17 514 1151 449]);
            cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
            mpref_sig = Modality_pref_p_value_all(1,:,tt) < 0.05;
            
            h = LinearCorrelation({
                (Modality_pref_all(1, monkey1 & ~cpref_sig & ~mpref_sig,tt)) ;
                (Modality_pref_all(1, monkey2 & ~cpref_sig & ~mpref_sig,tt)) ;
                (Modality_pref_all(1, monkey1 & xor(cpref_sig , mpref_sig),tt));
                (Modality_pref_all(1, monkey2 & xor(cpref_sig , mpref_sig),tt));
                (Modality_pref_all(1, monkey1 & cpref_sig & mpref_sig,tt));...
                (Modality_pref_all(1, monkey2 & cpref_sig & mpref_sig,tt))},...
                {
                (Choice_pref_all(k,monkey1 & ~cpref_sig & ~mpref_sig,tt)) ;
                (Choice_pref_all(k,monkey2 & ~cpref_sig & ~mpref_sig,tt)) ;
                (Choice_pref_all(k,monkey1 & xor(cpref_sig , mpref_sig),tt)) ;
                (Choice_pref_all(k,monkey2 & xor(cpref_sig , mpref_sig),tt)) ;
                (Choice_pref_all(k,monkey1 & cpref_sig & mpref_sig,tt)) ;...
                (Choice_pref_all(k,monkey2 & cpref_sig & mpref_sig,tt)) },...
                'CombinedIndex',[48 63],...
                'Ylabel','Choice preference (pre)','Xlabel','Modality preference (pre)',...
                'FaceColors',{'w','w',colors(k,:)*0.2 + [0.8 0.8 0.8],colors(k,:)*0.2 + [0.8 0.8 0.8],colors(k,:),colors(k,:)},'Markers',{'o','^'},...
                'LineStyles',{'k:','k:','k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,...
                'Method','Pearson','FittingMethod',2); figN = figN + 1;
            
            delete([h.group(1:6).line h.diag]); 
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       SetFigure(15);
            axis([-1 1 -1 1])
            set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1);
        
           
            % Annotate tcells
            plot(Modality_pref_all(1,select_tcells(select_bottom_line),tt),...
                Choice_pref_all(k,select_tcells(select_bottom_line),tt),'+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(Modality_pref_all(1,:,tt),Choice_pref_all(k,:,tt),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
            
            % Override histogram
            cla(h.ax_xhist);
            HistComparison( {Modality_pref_all(1,mpref_sig,tt)', Modality_pref_all(1,~mpref_sig,tt)'},...
                'EdgeColors',{'k','k'}, 'FaceColors',{'k', 'none'},...
                'MeanType','Median','Axes',h.ax_xhist,'TTest',0);
            h.ax_xhist.XLim = xlim(h.ax_raw);
            
            %     for i = 1:3
            %         set(h.group(i).dots,'color',colors(k,:));
            %     end
        end
    end
        
    function f3p2p2(debug)      % Correlations 2. Choice preference between modalities
        if debug
            dbstack;
            keyboard;
        end
       %% ---------------- 2 Choice pref (vest and visual) -------------------
        
        tt = 3;
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 5; monkey1 = monkey1(select_bottom_line)';
        monkey2 = monkeys == 10; monkey2 = monkey2(select_bottom_line)';

        
        set(figure(figN),'name','CDiv (visual) vs. CDiv (vest)','pos',[17 514 1151 449]);
        
        cpref_sig_1 = Choice_pref_p_value_all(1,:,tt) < 0.05;
        cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
        
        % Abs() or not?
        Choice_pref_all_temp = abs(Choice_pref_all);
        
        h = LinearCorrelation({
                (Choice_pref_all_temp(2, monkey1 & ~cpref_sig_1 & ~cpref_sig_2,tt)) ;
                (Choice_pref_all_temp(2, monkey2 & ~cpref_sig_1 & ~cpref_sig_2,tt)) ;
                (Choice_pref_all_temp(2, monkey1 & xor(cpref_sig_1 , cpref_sig_2),tt));
                (Choice_pref_all_temp(2, monkey2 & xor(cpref_sig_1 , cpref_sig_2),tt));
                (Choice_pref_all_temp(2, monkey1 & cpref_sig_1 & cpref_sig_2,tt));...
                (Choice_pref_all_temp(2, monkey2 & cpref_sig_1 & cpref_sig_2,tt))},...
                {
                (Choice_pref_all_temp(1,monkey1 & ~cpref_sig_1 & ~cpref_sig_2,tt)) ;
                (Choice_pref_all_temp(1,monkey2 & ~cpref_sig_1 & ~cpref_sig_2,tt)) ;
                (Choice_pref_all_temp(1,monkey1 & xor(cpref_sig_1 , cpref_sig_2),tt)) ;
                (Choice_pref_all_temp(1,monkey2 & xor(cpref_sig_1 , cpref_sig_2),tt)) ;
                (Choice_pref_all_temp(1,monkey1 & cpref_sig_1 & cpref_sig_2,tt)) ;...
                (Choice_pref_all_temp(1,monkey2 & cpref_sig_1 & cpref_sig_2,tt)) },...
                'CombinedIndex',[63],'PlotCombinedOnly', 1, ...
                'Ylabel','Vestibular choice preference','Xlabel','Visual choice preference',...
                'FaceColors',{'none','none',[0.8 0.8 0.8],[0.8 0.8 0.8],'k','k'},'Markers',{'o','^'},...
                'LineStyles',{'k:','k:','k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,'Method','Pearson','FittingMethod',2, ...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1); figN = figN + 1;
        
        delete([h.group(1:6).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       %  SetFigure(20);
        
        % Annotate tcells
        plot(Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt),Choice_pref_all_temp(1,select_tcells(select_bottom_line),tt),...
            '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        
        axis([-1 1 -1 1])
        set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(Choice_pref_all(2,:,tt),Choice_pref_all(1,:,tt),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
        
        clear Choice_pref_all_temp;
        
        
       %% ---------------- 3 Choice pref (single vs comb) -------------------
        two_face_colors = fliplr({'none','none',[0.8 0.8 1],[0.8 0.8 1],colors(1,:),colors(1,:);
                                  'none','none',[1 0.8 0.8],[1 0.8 0.8],colors(2,:),colors(2,:)});
        
        for k = 1:2  % Plot it separately
            
            
            set(figure(figN),'name','CDiv (single) vs. CDiv (comb)','pos',[17 514 1151 449]);
            
            cpref_sig_k = Choice_pref_p_value_all(k,:,tt) < 0.05;
%             cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
            cpref_sig_3 = Choice_pref_p_value_all(3,:,tt) < 0.05;
            
            % Abs() or not?
            Choice_pref_all_temp = (Choice_pref_all);
            
           
            h = LinearCorrelation({
                (Choice_pref_all_temp(k, monkey1 & cpref_sig_k & cpref_sig_3,tt));
                (Choice_pref_all_temp(k, monkey2 & cpref_sig_k & cpref_sig_3,tt));
                (Choice_pref_all_temp(k, monkey1 & xor(cpref_sig_k , cpref_sig_3),tt));
                (Choice_pref_all_temp(k, monkey2 & xor(cpref_sig_k , cpref_sig_3),tt));
                (Choice_pref_all_temp(k, monkey1 & ~cpref_sig_k & ~cpref_sig_3,tt)) ;
                (Choice_pref_all_temp(k, monkey2 & ~cpref_sig_k & ~cpref_sig_3,tt)) ;
                },...
                {
                (Choice_pref_all_temp(3,monkey1 & cpref_sig_k & cpref_sig_3,tt)) ;
                (Choice_pref_all_temp(3,monkey2 & cpref_sig_k & cpref_sig_3,tt));
                (Choice_pref_all_temp(3,monkey1 & xor(cpref_sig_k, cpref_sig_3),tt)) ;
                (Choice_pref_all_temp(3,monkey2 & xor(cpref_sig_k, cpref_sig_3),tt));
                (Choice_pref_all_temp(3,monkey1 & ~cpref_sig_k & ~cpref_sig_3,tt)) ;
                (Choice_pref_all_temp(3,monkey2 & ~cpref_sig_k & ~cpref_sig_3,tt)) ;
                },...
                'CombinedIndex',[63],'PlotCombinedOnly', 1, ...
                'Ylabel','Combined choice preference','Xlabel','Single choice preference',...
                'FaceColors',two_face_colors(k,:),'Markers',{'o','^'},...
                'LineStyles',{'k:','k:','k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,...
                'Method','Pearson','FittingMethod',2); figN = figN + 1;
            
            % delete([h.group(1:6).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       %   SetFigure(20);
            axis([-1 1 -1 1]); set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1); axis square;
            
            % Annotate tcells
            plot((Choice_pref_all_temp(k,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
                '+','markersize',tcell_cross_size,'color','k','linew',2);
%             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
%                 '+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot([Choice_pref_all_temp(k,:,tt)],[Choice_pref_all_temp(3,:,tt)],'visible','off'); hold on;
            set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_cpref_mpref],1});
            
            % Plotting histograms manually because we have to deal with different significant standards here. HH20180613
            cla(h.ax_xhist);
            HistComparison( {Choice_pref_all_temp(k,cpref_sig_k,tt)', Choice_pref_all_temp(k,~cpref_sig_k,tt)'},...
                'EdgeColors',{'k','k'}, 'FaceColors',{colors(k,:), 'none'},...
                'MeanType','Median','Axes',h.ax_xhist,'TTest',0);
            xlim([-1 1])
            cla(h.ax_yhist);
            HistComparison( {Choice_pref_all_temp(3,cpref_sig_3,tt)', Choice_pref_all_temp(3,~cpref_sig_3,tt)'},...
                'EdgeColors',{'k','k'}, 'FaceColors',{colors(3,:), 'none'},...
                'MeanType','Median','Axes',h.ax_yhist,'TTest',0);
            xlim([-1 1])
                       
            clear Choice_pref_all_temp;
            axis tight;
        end
        %             for i = 1:3
        %                 set(h.group(i).dots,'color',colors(k,:));
        %             end
        
        %% ===  (Comb - visual) VS (Comb - vest) %HH20160830 ===
        set(figure(figN),'name','Cdiv(Comb - visual) VS Cdiv(Comb - vest)','pos',[17 514 1151 449]);
                
        % Abs() or not?
        Choice_pref_all_temp_comb_minus_vest = abs(Choice_pref_all(3,:,:))-abs(Choice_pref_all(1,:,:));
        Choice_pref_all_temp_comb_minus_vis = abs(Choice_pref_all(3,:,:))-abs(Choice_pref_all(2,:,:));
        
        cellTypes = [group_result.Waveform_broad];
        
        h = LinearCorrelation({
            (Choice_pref_all_temp_comb_minus_vest(1,monkey1 & ~cellTypes,tt)) ;
            (Choice_pref_all_temp_comb_minus_vest(1,monkey2 & ~cellTypes,tt)) ;
            (Choice_pref_all_temp_comb_minus_vest(1,monkey1 & cellTypes,tt)) ;
            (Choice_pref_all_temp_comb_minus_vest(1,monkey2 & cellTypes,tt)) ;
            },...
            {
            (Choice_pref_all_temp_comb_minus_vis(1,monkey1 & ~cellTypes,tt)) ;
            (Choice_pref_all_temp_comb_minus_vis(1,monkey2 & ~cellTypes,tt)) ;
            (Choice_pref_all_temp_comb_minus_vis(1,monkey1 & cellTypes,tt)) ;
            (Choice_pref_all_temp_comb_minus_vis(1,monkey2 & cellTypes,tt)) ;
            },...
            'CombinedIndex',15,'PlotCombinedOnly', 1, ...
            'Xlabel','Combined - Vest (Choice preference)','Ylabel','Combined - Visual',...
            'FaceColors',{'none','none','k','k'},'Markers',{'o','^'},...
            'LineStyles',{'k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
            'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
            'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        axis([-1 1 -1 1]); set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1); axis square;
        
        % Annotate tcells
        plot((Choice_pref_all_temp_comb_minus_vest(1,select_tcells(select_bottom_line),tt)),...
            (Choice_pref_all_temp_comb_minus_vis(1,select_tcells(select_bottom_line),tt)),...
            '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        %             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
        %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot([Choice_pref_all_temp_comb_minus_vest(1,:,tt)],[Choice_pref_all_temp_comb_minus_vis(1,:,tt)],'visible','off'); hold on;
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_cpref_mpref],1});
        
        clear Choice_pref_all_temp;
        axis tight;
        
        %% ===  Mean raw delta firing (Comb - visual) VS (Comb - vest) %HH20170719 ===
        set(figure(figN),'name','Cdiv(Comb - visual) VS Cdiv(Comb - vest)','pos',[17 514 1151 449]);
                
        % Abs() or not?
        mean_raw_delta_firing = mean(PSTH_all_raw_PrefminusNull{1}(:,0<=rate_ts{1} & rate_ts{1}<=1500,:),2);
        mean_raw_delta_firing_comb_minus_vest = abs(mean_raw_delta_firing(:,3)) - abs(mean_raw_delta_firing(:,1));
        mean_raw_delta_firing_comb_minus_vis = abs(mean_raw_delta_firing(:,3)) - abs(mean_raw_delta_firing(:,2));
        
        h = LinearCorrelation({
            (mean_raw_delta_firing_comb_minus_vest(monkey1)) ;
            (mean_raw_delta_firing_comb_minus_vest(monkey2)) ;
            },...
            {
            (mean_raw_delta_firing_comb_minus_vis(monkey1)) ;
            (mean_raw_delta_firing_comb_minus_vis(monkey2)) ;
            },...
            'CombinedIndex',[3],'PlotCombinedOnly', 1, ...
            'Xlabel','Combined - Vest (raw PSTH)','Ylabel','Combined - Visual',...
            'FaceColors',{'k'},'Markers',{'o','^'},...
            'LineStyles',{'k:','k:','k-'},'MarkerSize',marker_size,...
            'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
            'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        axis([-1 1 -1 1]); axis square;
        
        % Annotate tcells
        plot((mean_raw_delta_firing_comb_minus_vest(select_tcells(select_bottom_line))),...
            (mean_raw_delta_firing_comb_minus_vis(select_tcells(select_bottom_line))),...
            '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        %             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
        %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot([mean_raw_delta_firing_comb_minus_vest],[mean_raw_delta_firing_comb_minus_vis],'visible','off'); hold on;
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_cpref_mpref],1});
        
        clear Choice_pref_all_temp;
        axis tight;
        
    end

    function f3p2p3(debug)      %  Choice pref (pre vs post)
        if debug
            dbstack;
            keyboard;
        end
        
        %% ---------------- Choice pref (pre vs post) -------------------
        
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
            'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;         SetFigure(20);
        
        
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

    function f3p4(debug)      % Correlations 4. Cell position and cell type HH20170715
        if debug
            dbstack;
            keyboard;
        end
        
        all_monkey_hemis = Position_all(:,1:2);
        all_APs = Position_all(:,3);
        all_VDs = Position_all(:,4);
        all_depths = Position_all(:,5);
        all_spikewid = [group_result(select_bottom_line).Waveform_peakToTrough]';
        all_celltype = [group_result(select_bottom_line).Waveform_broad]';
        allmonkey = xls_num{1}(select_bottom_line,header.Monkey); 
        
        % === Do some analyses ===
        %% 1. Draw depths
        
        figure(7161901);  clf
        set(gcf,'uni','norm','pos',[0.429       0.061       0.177       0.486]);
        [counts,centers] = hist(all_depths,29);
        
        layerKmeans = kmeans(all_depths,3,'Replicates',10,'Display','final');
        layerColors = repmat([0.9;0.5;0],1,3);
        
        for ll = 1:3
            countsThis = hist(all_depths(layerKmeans == ll), centers);
            set(bar(centers, countsThis,1),'facecolor',layerColors(ll,:),'linew',0.1); hold on;
        end
        
        
        fitOptions = fitoptions('gauss3','StartPoint',[5 700 100 5 1100 50 5 1700 50]);
        fit_model = fit(centers',counts','gauss3',fitOptions);
        gauss_f = @(x,a,b,c)a*exp(-((x-b)/c).^2);
        
       %  bar(centers,counts,1); 
        xx = linspace(min(xlim),max(xlim),500);
        hold on; plot(xx,fit_model(xx)','r'); xlabel('depth');
        plot(xx,gauss_f(xx,fit_model.a1,fit_model.b1,fit_model.c1)','b');
        plot(xx,gauss_f(xx,fit_model.a2,fit_model.b2,fit_model.c2)','b');
        plot(xx,gauss_f(xx,fit_model.a3,fit_model.b3,fit_model.c3)','b');
        view(90,90)
        
       %% 2. Draw waveform width 
        
       % Align all spikes
       figure(616); clf; hold on
       set(gcf,'uni','norm','pos',[0.426        0.62        0.25       0.304]);
       
       waveformPlotBefore = 0.6; %ms
       waveformPlotAfter = 1; %ms
       waveformPlotTs = - waveformPlotBefore : interp_dt: waveformPlotAfter;
       
       allAlignedWaveform = nan(N, (waveformPlotBefore + waveformPlotAfter)/interp_dt + 1);
       allBackgroundFiring = nan(N,1);
       
       for cc = 1:length(group_result)
           
           % -- This waveform --
           waveform_this = group_result(cc).Waveform;           
           waveform_this = [nan(1,1000) waveform_this nan(1,1000)]; % Padding of nans
           
           alignTime = group_result(cc).Waveform_peak + 1000;
           
           waveformAligned = waveform_this(alignTime - waveformPlotBefore/interp_dt : ...
               alignTime + waveformPlotAfter/interp_dt);
           
           if group_result(cc).Waveform_broad
               plot( waveformPlotTs , waveformAligned,'color', [0.8 0.8 0.8]);
           else
               plot( waveformPlotTs , waveformAligned,'color',hsv2rgb(rgb2hsv([1 0 0])-[0 0.7 0]));
           end
           %             plot((waveform_t1-waveform_t_left)*dt,waveform_this(waveform_t1),'or');
           %             plot((waveform_t_right-waveform_t_left)*dt,waveform_this(waveform_t_right),'og');
           
           allAlignedWaveform(cc,:) = waveformAligned;
           
           % -- This mean firing rate --
           fromStartToWhen = 100; % in ms
           dt = group_result(cc).mat_raw_PSTH.spike_aligned{2,1}(2) - group_result(cc).mat_raw_PSTH.spike_aligned{2,1}(1);
           allBackgroundFiring(cc) = mean(sum(group_result(cc).mat_raw_PSTH.spike_aligned{1,1}(:,1:fromStartToWhen),2)/fromStartToWhen*1000); % in Hz
           

       end
       
       meanBroadWaveform = nanmean(allAlignedWaveform ([group_result.Waveform_broad],:),1);
       seBroadWaveform = nanstd(allAlignedWaveform ([group_result.Waveform_broad],:),[],1) ./ 1;%...
                         %sqrt(sum(~isnan(allAlignedWaveform([group_result.Waveform_broad],:)),1));

       meanNarrowWaveform = nanmean(allAlignedWaveform (~ [group_result.Waveform_broad],:),1);
       seNarrowWaveform = nanstd(allAlignedWaveform (~ [group_result.Waveform_broad],:),[],1) ./1; % ...
                       %  sqrt(sum(~isnan(allAlignedWaveform(~ [group_result.Waveform_broad],:)),1));

       shadedErrorBar(waveformPlotTs, meanBroadWaveform, seBroadWaveform,'lineprops',{'color','k'});
       shadedErrorBar(waveformPlotTs, meanNarrowWaveform, seNarrowWaveform,'lineprops',{'color','r'});
                     
       xlabel('ms');
       xlim([-0.6 1.1])
       

        
       %% 3. Waveform width distribution
        figure(18616);  clf; hold on;
        set(gcf,'uni','norm','pos',[0.033       0.631       0.389       0.289]);

        [counts,centers] = hist(all_spikewid,30);
        
        cellKmeans = kmeans(all_spikewid,2,'Replicates',10,'Display','final');
        cellColors = repmat([0.8;0.2],1,3);
        
        for ll = 1:2
            countsThis = hist(all_spikewid(cellKmeans == ll), centers);
            set(bar(centers, countsThis,1),'facecolor',cellColors(ll,:),'linew',0.01); hold on;
        end
        
        fitOptions = fitoptions('gauss2','StartPoint',[5 0.2 0.1 5 0.4 0.1]);
        fit_model = fit(centers',counts','gauss2',fitOptions);
        gauss_f = @(x,a,b,c)a*exp(-((x-b)/c).^2);
        
       %  bar(centers,counts,1); 
        xx = linspace(min(xlim),max(xlim),500);
        hold on; plot(xx,fit_model(xx)','r'); xlabel('depth');
        plot(xx,gauss_f(xx,fit_model.a1,fit_model.b1,fit_model.c1)','b');
        plot(xx,gauss_f(xx,fit_model.a2,fit_model.b2,fit_model.c2)','b');
        xlim([0.2 0.6])
        
        % TTest of two cell types
        meanNarrow = mean(all_spikewid(~all_celltype));
        stdNarrow = std(all_spikewid(~all_celltype));
        meanBroad = mean(all_spikewid(all_celltype));
        stdBroad = std(all_spikewid(all_celltype));
        [~,p] = ttest2(all_spikewid(~all_celltype),all_spikewid(all_celltype));
        
        plot([meanNarrow meanBroad],max(ylim)*0.9*[1 1],'o');
        
        herrorbar(meanNarrow,max(ylim)*0.9,stdNarrow); text(meanNarrow,max(ylim),sprintf('n = %g',sum(~all_celltype)));
        herrorbar(meanBroad,max(ylim)*0.9,stdBroad); text(meanBroad,max(ylim),sprintf('n = %g, p = %.2g',sum(all_celltype),p));
        
       
       %% 4. Width and depth: 2-D

        colorsMaxMod = [1 1 1; colors];   
        monkeys = [5 10];
        monkeyMarkers = {'o','^'};
        
        anySig = any(Choice_pref_p_value_all(:,:,3) < 0.05);
        allMaxMod = zeros(sum(select_bottom_line),1); % None sig = 0
        [~,MaxMod] = max(abs(Choice_pref_all(:,:,3)),[],1); % Find the max modality
        allMaxMod(anySig) = MaxMod(anySig);
        
        % 
        figure(1903); clf; hold on;
        set(gcf,'uni','norm','pos',[0.035       0.059       0.392        0.49]);
        
        for mm = 1:2
            for kk = 0:3
                if ~isempty(allmonkey == monkeys(mm))
                    thisToPlot = (allmonkey == monkeys(mm)) & (allMaxMod == kk);
                    plot(all_spikewid(thisToPlot), all_depths(thisToPlot),...
                        monkeyMarkers{mm},'markerfacecol',colorsMaxMod(kk+1,:),'markeredgecol','k','markersize',10);
                end
            end
        end
        axis ij
        
        %{
        h = LinearCorrelation({
            all_spikewid
            },...
            {
            all_depths
            },...
            'CombinedIndex',[],...
            'Xlabel','Spike width','Ylabel','Depth',...
            'FaceColors',{'none','k'},'Markers',{'o'},...
            'LineStyles',{'k:','k-.','k-'},'MarkerSize',12,...
            'figN',1948,'XHist',30,'YHist',34,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman','FittingMethod',2);
        axis ij
        view(h.ax_yhist,90,90)
        %}
        
        %% AP and depth
        %{
            LinearCorrelation({
            all_APs
            },...
            {
            all_depths
            },...
            'CombinedIndex',[],...
            'Xlabel','AP','Ylabel','Depth',...
            'FaceColors',{'none','k'},'Markers',{'o'},...
            'LineStyles',{'k:','k-.','k-'},'MarkerSize',12,...
            'figN',194852,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman');
        %}
        
        %% 5. Spike waveform and baseline firing rate (HH20190330)
        figure(19329); clf;
        cache = nan(max(sum(all_celltype),sum(~all_celltype)),2);
        cache(1:sum(~all_celltype),1) = allBackgroundFiring(~all_celltype);  % Narrow
        cache(1:sum(all_celltype),2) = allBackgroundFiring(all_celltype); % Broad
        result = BarComparison(cache,'Fig',19329);
        text(1.5, max(ylim), sprintf('p = %g',result.ps_ttest(2,end)));
        set(gca,'xticklabel',{'narrow','broad'});
        ylabel('Background activity');
        
        figure(193291); clf;
        LinearCorrelation(all_depths, allBackgroundFiring,'Fig',193291,'Xlabel', 'Depth','Ylabel','Background Firing');
       
        %% 6. Spike waveform and contralateral (HH20190330)
        all_if_contra = [group_result.if_contralateral]';
        narrow_contra = sum(~all_celltype & all_if_contra)/sum(~all_celltype);
        broad_contra = sum(all_celltype & all_if_contra)/sum(all_celltype);
        figure(190331); clf; 
        bar([1 2],[narrow_contra broad_contra]);
        set(gca,'xticklabel',{'narrow','broad'});
        ylabel('contra %');
        
        %% Draw AP,DV
%         %{
        figure(7162221); 
        subplot(2,2,1);
        plot(all_VDs(all(bsxfun(@eq,all_monkey_hemis,[5 1]),2)),all_APs(all(bsxfun(@eq,all_monkey_hemis,[5 1]),2)),'o');
        axis equal; axis([10 18 -7 6]);
        title('Polo L');
        subplot(2,2,2);
        plot(all_VDs(all(bsxfun(@eq,all_monkey_hemis,[5 2]),2)),all_APs(all(bsxfun(@eq,all_monkey_hemis,[5 2]),2)),'o');
        axis equal;axis([10 18 -7 6]);
        title('Polo R');
        subplot(2,2,3);
        plot(all_VDs(all(bsxfun(@eq,all_monkey_hemis,[10 1]),2)),all_APs(all(bsxfun(@eq,all_monkey_hemis,[10 1]),2)),'o');
        axis equal;axis([10 18 -7 6]);
        title('Messi L');
        subplot(2,2,4);
        plot(all_VDs(all(bsxfun(@eq,all_monkey_hemis,[10 2]),2)),all_APs(all(bsxfun(@eq,all_monkey_hemis,[10 2]),2)),'o');
        axis equal;axis([10 18 -7 6]);
        title('Messi R');
        %}

    end

    function f3p4p1(debug)      % Cell type and position
        if debug
            dbstack;
            keyboard;
        end
        
        all_monkey_hemis = Position_all(:,1:2);
        all_APs = Position_all(:,3);
        all_VDs = Position_all(:,4);
        all_depths = Position_all(:,5);
        all_spikewid = [group_result(select_bottom_line).Waveform_peakToTrough]';
        all_celltype = [group_result(select_bottom_line).Waveform_broad]';
        allmonkey = xls_num{1}(select_bottom_line,header.Monkey);

       %% 5. Multi-linear regression
        if isempty( firstDivergenceTime )
            f1p1p6p1(0)
        end
        
        XsName = {'monkey','celltype','depth','AP','MidLateral','typeXdepth'};        
        Xs = [ allmonkey, all_celltype , all_depths, all_APs, all_VDs ,  all_celltype.* all_depths];
        norm_Xs = (Xs - nanmean(Xs)) ./ nanstd(Xs);
        
%         XsInUse = [ 1 2 3 4 5];
        XsInUse = [ 2 3 4 5];
        
        Ys = [MemSac_indicator(select_bottom_line,:), ...  % 1
              group_MemSac_AI(select_bottom_line,3), ...   % 2
              (Modality_pref_all(:,:,3))', ...               % 3-5 
              abs(Choice_pref_all(:,:,3))', ...            % 6-8 
              firstDivergenceTime(select_bottom_line,:)];  % 9-11
        YsName = {'1.MemsacDDI','2.MemsacAI','3.ModPref1','4.ModPref2','5.ModPref3',...
                  '6.CPref1','7.CPref2','8.CPref3','9.divTime1','10.divTime2','11.divTime3'};  
        
        ps = []; ps2 = []; 
        fitlmModel = {}; fitlmModelSimplified = {};
        
        for yy = 1:size(Ys,2)
            % --- Linear regression ---
            
            % Use glmfit
            [~,~,STATS] = glmfit( norm_Xs(:,XsInUse) , Ys(:,yy));
            ps(yy,:) = STATS.p;
            
            % Use fitlm for comparison
            mdl = fitlm(norm_Xs(:,XsInUse), Ys(:,yy));
            fitlmModel{yy} = mdl;
            ps2(yy,:) = mdl.anova{1:end-1,end}';
            
            % Try simplification
            mdlSimple = step(mdl,'Nsteps',10,'Criterion','BIC');
            fitlmModelSimplified{yy} = mdlSimple;
            
            %             [B,~,STATS] = glmfit( [XsANOVA(:,1:3)] , Ys(:,yy));
            %             ps(yy,:) =  STATS.p;
            
            % --- ANOVA ---
            %anovan            
        end
        ps(:,1) = []; % Remove baseline
        tmp = array2table(ps,'VariableNames',XsName(XsInUse),'RowNames',YsName)

        % -- show simplified --
        for yy = 1:size(Ys,2)
            nametmp = fitlmModelSimplified{yy}.CoefficientNames(2:end);
            fprintf('%g: ', yy);
            fprintf('%s, ', nametmp{:});  fprintf(':');
            fprintf('%g\t',fitlmModelSimplified{yy}.Coefficients.pValue(2:end)');
            fprintf('\n');
        end
        
        % Plot significants
        significants = find(ps<0.05);
        nSig = length(significants);
        figure(4603); clf ; set(gcf,'uni','norm','pos',[0.009       0.045       0.552       0.861]);
        hs = tight_subplot(ceil(nSig/ceil(nSig/sqrt(nSig))),ceil(nSig/sqrt(nSig)),[.1 .07],[.03 .07],[.07 .03]);
        
        for ss = 1:nSig
            [thisY, thisX] = ind2sub(size(ps),significants(ss));
            h = LinearCorrelation({ Xs(allmonkey == 5,XsInUse(thisX)), ...
                                    Xs(allmonkey == 10,XsInUse(thisX))},...
                                  { Ys(allmonkey == 5,thisY),...
                                    Ys(allmonkey == 10, thisY)},...
                            'Combined',3, 'FittingMethod',1,'LineStyles',{'k--','r--','k-'},'Markers',{'o','^'},...
                            'PlotCombinedOnly', XsInUse(thisX) == 1 ,... % If monkey is X axis, no fitting for each monkey
                            'xlabel',XsName{XsInUse(thisX)},'ylabel',YsName{thisY},'Axes',hs(ss));
            title(sprintf('(%.2g): %s',ps(significants(ss)),num2str([h.group(:).p],'%.2g, ')));
            legend off;
        end
        
       %% 6. ANOVA
        
        
        
        
       %% 7. Pair-wise Correlations 
        
        % 1. Depth and Choice_pref
        figure(7162036); clf;
        
        % Position = [monkey hemi AP VD depth whichlayer]
        which_mod = 1; which_pos = 5;
        cpref_sig = Choice_pref_p_value_all(which_mod,:,3) < 0.05;

        LinearCorrelation({
            Position_all(~cpref_sig,which_pos);
            Position_all(cpref_sig,which_pos);
            },...
            {
            (abs(Choice_pref_all(which_mod,~cpref_sig,3)));
            (abs(Choice_pref_all(which_mod,cpref_sig,3)));
            },...
            'CombinedIndex',3,'PlotCombinedOnly',0,...
            'Xlabel','','Ylabel','',...
            'FaceColors',{'none','k'},'Markers',{'o'},...
            'LineStyles',{'k:','k-.','k-'},'MarkerSize',12,...
            'figN',7162036,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman','FittingMethod',2);
  
        %% 8. Vest divergence time comparison
        figure(162015); clf; hold on;
        for k = 1:3
            vestDivBroad = Ys(all_celltype,8 + k);
            nBroad = sum(~isnan(vestDivBroad));
            vestDivNarrow = Ys(~all_celltype, 8 + k);
            nNarrow = sum(~isnan(vestDivNarrow));
            
            means = [nanmean(vestDivBroad), nanmean(vestDivNarrow)];
            sems = [nanstd(vestDivBroad)/sqrt(nBroad), nanstd(vestDivNarrow)/sqrt(nNarrow)];
            [~,p] = ttest2(vestDivBroad,vestDivNarrow);
            
            bar([(k-1)*2+1 k*2], means); 
            errorbar([(k-1)*2+1 k*2], means, sems);
            text(k*2-0.5,mean(means)*1.1,sprintf('%g, %g\np = %g',nBroad,nNarrow,p));
        end
          
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
        A_memSac = group_MemSac_DDI(select_bottom_line,2:end);  % Read memsac from .mat file. HH20150413
        A_memSac = (A_memSac-0.5)*2; % Map to -1 ~ 1 for better visualization
        
        % Now I add mod_div. HH20150415
        % A_moddiv = reshape(ModDiv_All{j_PCA_A}(select_bottom_line,:,:),sum(select_bottom_line),[]);
        A_moddiv = reshape(ModDiv_All{j_PCA_A}(select_bottom_line,:,1),sum(select_bottom_line),[]); % Only vis-vest
        
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
        
        for sort_ind = length(sort_method):-1:1
            
            sort_begin = sort_method{sort_ind}(1);
            sort_end = sort_method{sort_ind}(end);
            
            mean_for_sort = mean(PCA_A(:,sort_begin:sort_end),2);
            mean_for_sort(isnan(mean_for_sort)) = -inf;  % Nan goes last
            
            [~,sort_order] = sort(mean_for_sort,'descend');
            
            % Enlarge mem-sac part for clarity
            A_forplot = [reshape(repmat(A_memSac,enlarge_factor,1),size(A_memSac,1),[]) A_choicediv A_moddiv A_moddiv(:,end)];% A_CP A_CP(:,end)]; % To ensure the last value to be shown
            A_forplot = A_forplot(sort_order,:);
            
            A_forplot = [A_forplot; A_forplot(end,:)]; % To ensure the last cell to be shown
            
            A_forplot(isnan(A_forplot)) = -0.5; % This workaround is to avoid that cells which are next the NaNs cannot be displayed.
            
            % [nn, tt] = meshgrid(1:size(A_forplot,1), 1:size(A_forplot,2));
            
            set(figure(2099+figN),'name',['CDiv and MDiv, hot-gram, j_PCA_A = ' num2str(j_PCA_A)],...
                'unit','norm','pos',[0.581547619047619 0.177142857142857 0.414880952380952 0.738095238095238]); 
            clf; figN = figN+1;
            
            % h1 = surf(tt,nn,A_forplot','Edgecolor','none');
            h1 = imagesc(A_forplot); colormap jet
            % axis xy;
            
            % Annotate time separations
            hold on;
            CD_begin_at = enlarge_factor * 5 + 1;
            
            % Stim on / stim off / sac on
            for ttt = 1:3
                for tt = 0: 3% 5
                    plot(CD_begin_at + tt*length(rate_ts{j_PCA_A}) + find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,ttt),1) * [1 1],...
                        [1 size(A_forplot,1)],'k','linesty',marker_for_time_markers{j_PCA_A}{ttt},'linewid',1.5);
                end
            end
            
            plot3([CD_begin_at CD_begin_at],[1 size(A_forplot,1)],[1 1],'k','linewid',3);
            for tt = 1: 4 % 6
                plot(CD_begin_at + tt*[length(rate_ts{j_PCA_A})  length(rate_ts{j_PCA_A})],[1 size(A_forplot,1)],'k','linewid',3);
            end
            
            % Annotate tcells
            cell_loc = select_tcells(select_bottom_line);
            h2 = plot(0,find(cell_loc(sort_order))+.5,'+','markersize',5,'color','k','linew',1.5);
            
            % Annotate sort methods
            temp_beg = @(x)((x<=5)*((x-1)*enlarge_factor) + 1 + (x>5)*(CD_begin_at+(x-5)-1));
            temp_end = @(x)((x<=5)*((x)*enlarge_factor) + 1 + (x>5)*(CD_begin_at+(x-5)));
            
            sort_begin_forplot = temp_beg(sort_begin);
            sort_end_forplot = temp_end(sort_end);
            
            plot([sort_begin_forplot sort_end_forplot],[size(A_forplot,1)+1 size(A_forplot,1)+1],'color',colors(2,:),'linewid',5);
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
            f4p0(debug);
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
        % Pack all data (PSTH) into matrix B: [vest_pref, vest_null, vis_pref, vis_null, comb_pref, comb_null]
        % Only decicion period is included
        
        find_PCA_B = find(select_for_PCA_B);
        
        PCA_B = nan(sum(select_for_PCA_B),6 * sum(PCA_B_time_range));
        
        for i = 1:sum(select_for_PCA_B)
            raw_PSTH = group_result(find_PCA_B(i)).mat_raw_PSTH.PSTH{j_PCA_B,ALL_CHOICE,1}.ys;
            if size(raw_PSTH,1)==6
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
            PCA_B_projPC{dim} = (reshape(weights_PCA_B_PC(:,dim)' * PCA_B,[],6))';
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
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %
        %         % Annotate tcells
        %         plot(PCA_B_PC(select_tcells(select_for_PCA_B),1),PCA_B_PC(select_tcells(select_for_PCA_B),2),'+','markersize',16,'color',colors(2,:),'linew',2);
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
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        %         delete([h.group(1:2).line]);
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %
        %         % Annotate tcells
        %         plot(PCA_B_PC(select_tcells(select_for_PCA_B),1),abs(Choice_pref_all(k,select_tcells(select_for_PCA_B),tt)),'+','markersize',16,'color',colors(2,:),'linew',2);
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
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        %
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %
        %         % Annotate tcells
        %         plot(mean(MemSac_DDI(select_tcells,MemSac_DDI_phase),2),PCA_B_PC(select_tcells(select_for_PCA_B),1),...
        %             '+','markersize',16,'color',colors(2,:),'linew',2);
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
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        %
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %         delete([h.group(1:2).line]);
        %         % Annotate tcells
        %         plot(mean(MemSac_DDI_selected(select_tcells(select_for_PCA_B),MemSac_DDI_phase),2),abs(Choice_pref_all(k,select_tcells(select_for_PCA_B),tt)),...
        %             '+','markersize',16,'color',colors(2,:),'linew',2);
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
        %             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        %         delete([h.group(1:2).line]);
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        %
        %         % Annotate tcells
        %         plot((Modality_pref_all(k,select_tcells(select_bottom_line),tt)), PCA_B_PC(select_tcells(select_bottom_line),2),...
        %             '+','markersize',16,'color',colors(2,:),'linew',2);
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
        plot((1:denoised_dim)',cumsum(PCA_B_explained(1:denoised_dim)),'ro-','markersize',8,'linew',1.5,'markerfacecol',colors(2,:));
        plot([0 1],[0 PCA_B_explained(1)],'r-','linew',1.5);
        plot(xlim,[1 1]*sum(PCA_B_explained(1:denoised_dim)),'r--');
        text(denoised_dim,sum(PCA_B_explained(1:denoised_dim))*0.9,[num2str(sum(PCA_B_explained(1:denoised_dim))) '%'],'color',colors(2,:));
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
            SeriesComparison(shiftdim(norm_proj_PC_this',-1),PCA_B_times,'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'axes',h);
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
        
        % ======== 2D ========= %  % Increase to 3d
        set(figure(2099+figN),'pos',[18 170 898 786],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        
        which_two_dimension = [1,2,3];
        
        % -- For plotting time markers
        plotInt = 200; % ms
        plotPerTimeBin = fix(plotInt/(PCA_B_times(2)-PCA_B_times(1))); % Should be 100/10 = 10
        plotMinInd = fix(-min(PCA_B_times)/plotInt)*plotPerTimeBin;
        plotInd = plotMinInd : plotPerTimeBin : length(PCA_B_times);
        %%
        cla
        for k = 1:3
            
            % Time markers
            start_time = 1; % Start point
            
            % -- Big ball
            % Pref
            
            % 3-d. Not good for plotting
            %{
            h_pref( k) = plot3(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,start_time),...
                              PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,start_time),...
                              PCA_B_projPC{which_two_dimension(3)}((k-1)*2+1,start_time),...
                             'o','color',colors(k,:),'markersize',20,'markerfacecol',colors(k,:));
            % Null
            h_null(k) = plot3(PCA_B_projPC{which_two_dimension(1)}(k*2,start_time),...
                             PCA_B_projPC{which_two_dimension(2)}(k*2,start_time),...
                             PCA_B_projPC{which_two_dimension(3)}(k*2,start_time),...
                             'o','color',colors(k,:),'markersize',20,'linew',3);
            
            % -- Tranjectory
            % Pref
            plot3(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,:),...
                 PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,:),...
                 PCA_B_projPC{which_two_dimension(3)}((k-1)*2+1,:),...
                 '-','color',colors(k,:),'linew',3);
            % Null
            plot3(PCA_B_projPC{which_two_dimension(1)}(k*2,:),...
                PCA_B_projPC{which_two_dimension(2)}(k*2,:),...
                PCA_B_projPC{which_two_dimension(3)}(k*2,:),...
                '--','color',colors(k,:),'linew',3);
            
            % -- Time markers
            colorsHsv = repmat(rgb2hsv(colors(k,:)),length(plotInd),1);
            colorsHsv(:,2) = linspace(0.3,1,length(plotInd));
            colorsHsv(:,3) = 1;
            colorsRGB = hsv2rgb(colorsHsv);
            for pp = 1:length(plotInd)
                plot3(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,plotInd(pp)),...
                      PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,plotInd(pp)),...
                      PCA_B_projPC{which_two_dimension(3)}((k-1)*2+1,plotInd(pp)),...
                      'o','color',colorsRGB(pp,:),'markerfacecol',colorsRGB(pp,:),'linew',2,'markersize',13);
                  
                plot3(PCA_B_projPC{which_two_dimension(1)}(k*2,plotInd(pp)),...
                      PCA_B_projPC{which_two_dimension(2)}(k*2,plotInd(pp)),...
                      PCA_B_projPC{which_two_dimension(3)}(k*2,plotInd(pp)),...
                      'o','color',colorsRGB(pp,:),'markerfacecol','none','linew',2,'markersize',13);
            end
            %}
            
            h_pref( k) = plot(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,start_time),...
                              PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,start_time),...
                             'o','color',colors(k,:),'markersize',20,'markerfacecol',colors(k,:));
            % Null
            h_null(k) = plot(PCA_B_projPC{which_two_dimension(1)}(k*2,start_time),...
                             PCA_B_projPC{which_two_dimension(2)}(k*2,start_time),...
                             'o','color',colors(k,:),'markersize',20,'linew',3);
            
            % -- Tranjectory
            % Pref
            plot(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,:),...
                 PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,:),...
                 '-','color',colors(k,:),'linew',3);
            % Null
            plot(PCA_B_projPC{which_two_dimension(1)}(k*2,:),...
                PCA_B_projPC{which_two_dimension(2)}(k*2,:),...
                '--','color',colors(k,:),'linew',3);
            
            % -- Time markers
            colorsHsv = repmat(rgb2hsv(colors(k,:)),length(plotInd),1);
            colorsHsv(:,2) = linspace(0.2,1,length(plotInd));
            
            if k == 3
                colorsHsv(:,3) = linspace(.9,colorsHsv(1,3),length(plotInd));
            end
            
            colorsRGB = hsv2rgb(colorsHsv);
            
            for pp = 1:length(plotInd)
                plot(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,plotInd(pp)),...
                      PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,plotInd(pp)),...
                      'o','color',colorsRGB(pp,:),'markerfacecol',colorsRGB(pp,:),'linew',0.1,'markersize',15);
                  
                plot(PCA_B_projPC{which_two_dimension(1)}(k*2,plotInd(pp)),...
                      PCA_B_projPC{which_two_dimension(2)}(k*2,plotInd(pp)),...
                      'o','color',colorsRGB(pp,:),'markerfacecol','none','linew',2,'markersize',15);
            end            
            
        end
        %%
        axis tight;  grid off;
        axis off
        % xlabel('PC1'); ylabel('PC2');
        
        % xlims = xlim; ylims = ylim;
        % h_text = text(xlims(1),ylims(2),'');
        
        
        % ======= Euclidean distance =========
        h_timeaxis = axes('pos', [0.026 0.616 0.344 0.363] ,'color','none');  hold on
        
        prefs = [1 3 5];
        nulls = [2 4 6];
        
        % distance = sqrt((PCA_B_projPC{1}(prefs,:) - PCA_B_projPC{1}(nulls,:)).^2 + (PCA_B_projPC{2}(prefs,:) - PCA_B_projPC{2}(nulls,:)).^2 + (PCA_B_projPC{3}(prefs,:) - PCA_B_projPC{3}(nulls,:)).^2);
        % distance = distance';
        
        distance = sum(reshape(cell2mat((cellfun(@(x)(x(prefs,:)-x(nulls,:)).^2',PCA_B_projPC,'uni',0))),[],3,denoised_dim),3);
        
        % Normalize distance
        % distance = (distance - repmat(min(distance,[],1),size(distance,1),1))./repmat(range(distance),size(distance,1),1);
        % set(figure(2099+figN),'pos',[18 355 766 601],'name',['Euclidean distance, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        
        for k = 1:3
            plot(repmat(PCA_B_times',1,3),distance(:,k),'color',colors(k,:),'linew',3);
        end
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_PCA_B}(1),Gauss_vel(:,2)*max(ylim)/4,'--','linew',3,'color',[0.6 0.6 0.6]);
        axis tight
        
        plot(repmat(PCA_B_times',1,3),sum(distance(:,[1 2]),2),'k--','linew',2);
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
        
        plot(repmat(PCA_B_times(1:end-1)',1,3), diff(distance),'linew',3); hold on;
        axis([-50 1800 -100 950]);
        
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
                            set(h_pref(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}((kk-1)*2+1,t_ind),...
                                           'ydata',PCA_B_projPC{which_two_dimension(2)}((kk-1)*2+1,t_ind),...
                                           'zdata',PCA_B_projPC{which_two_dimension(3)}((kk-1)*2+1,t_ind));
                            set(h_null(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}(kk*2,t_ind),...
                                           'ydata',PCA_B_projPC{which_two_dimension(2)}(kk*2,t_ind),...
                                           'zdata',PCA_B_projPC{which_two_dimension(3)}(kk*2,t_ind));
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
                        set(h_pref(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}((kk-1)*2+1,t_ind),...
                                       'ydata',PCA_B_projPC{which_two_dimension(2)}((kk-1)*2+1,t_ind),...
                                       'zdata',PCA_B_projPC{which_two_dimension(3)}((kk-1)*2+1,t_ind));
                        set(h_null(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}(kk*2,t_ind),...
                                       'ydata',PCA_B_projPC{which_two_dimension(2)}(kk*2,t_ind),...
                                       'zdata',PCA_B_projPC{which_two_dimension(3)}(kk*2,t_ind));
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

    function f51p1(debug)     % 1. PCA over heading and choice for each modality (Dora's method) %HH20160701
        if debug;  dbstack;  keyboard;  end
        
        find_PCA_B = find(select_for_PCA_B);
        
%         PCA_B_heading = nan(sum(select_for_PCA_B),6 * sum(PCA_B_time_range));
%         
%         for i = 1:sum(select_for_PCA_B)
%             raw_PSTH = group_result(find_PCA_B(i)).mat_raw_PSTH.PSTH{j_PCA_B,ALL_CHOICE,1}.ys;
%             if size(raw_PSTH,1)==6
%                 PCA_B_heading(i,:) = reshape(raw_PSTH(:,PCA_B_time_range)',[],1)';
%             end
%         end
        for k = 1:3
            PCA_B_heading_eachk{k} = PSTH_correct_angles_raw{j_PCA_B}(select_for_PCA_B,PCA_B_time_range,:,k);
            PCA_B_heading_eachk{k} = reshape(PCA_B_heading_eachk{k},size(PCA_B_heading_eachk{k},1),[]); % Cell * (0-pref,0-null,1-pref,1-null, ..., 8-pref, 8-null)
        end
        
        PCA_B_heading = [PCA_B_heading_eachk{1} PCA_B_heading_eachk{2} PCA_B_heading_eachk{3}];
        PCA_B_heading(sum(isnan(PCA_B_heading),2) > 0,:) = [];
        
        % Do PCA
        [weights_PCA_B_heading_PC, score, latent, ~, PCA_B_explained_heading] = pca(PCA_B_heading');
        
        % Manually rotate PC1 and PC2 (they're also orthogonal, so there's no
        % change in PCA results). @HH20150424
%         theta = 55 / 180 * pi; % Rotate counterclockwise the two PC bases
%         THETA = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
%         weights_PCA_B_PC(:,1:2) = weights_PCA_B_PC(:,1:2) * THETA;
%         weights_PCA_B_PC(:,2) = -weights_PCA_B_PC(:,2);
        
        % Projecting the raw data onto the first several eigenvectors
        for dim = 1:denoised_dim
            PCA_B_heading_projPC{dim} = (reshape(weights_PCA_B_heading_PC(:,dim)' * PCA_B_heading,[],3*size(PSTH_correct_angles_raw{j_PCA_B},3)))';
            %     projPC{dim} = reshape(score(:,dim),[],6)'
        end       
        
        % ======== 2D ========= %
        set(figure(2099+figN),'pos',[18 170 898 786],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        
        which_two_dimension = [1,2];
        
        for k = 1:3
%             k = 1;
            colors_angles{k} = colormap(gray);
            colors_angles{k} = ones(size(colors_angles{k},1),3) - colors_angles{k} .* repmat([1 1 1]-colors(k,:),size(colors_angles{k},1),1);
            colors_angles{k} = colors_angles{k}(round(linspace(20,length(colors_angles{k}),5)),:);
            
            for hh = 1:5 % 5 headings
                %             % Time markers
                start_time = 1; % Start point
                % Pref
                plot(PCA_B_heading_projPC{which_two_dimension(1)}(2*hh-1+(k-1)*10,start_time),...
                     PCA_B_heading_projPC{which_two_dimension(2)}(2*hh-1+(k-1)*10,start_time),...
                     'o','color',colors_angles{k}(hh,:),'markersize',20,'markerfacecol',colors_angles{k}(hh,:));
                % Null
                plot(PCA_B_heading_projPC{which_two_dimension(1)}(2*hh+(k-1)*10,start_time),...
                     PCA_B_heading_projPC{which_two_dimension(2)}(2*hh+(k-1)*10,start_time),...
                     'o','color',colors_angles{k}(hh,:),'markersize',20,'linew',3);
                
                % Pref
                plot(PCA_B_heading_projPC{which_two_dimension(1)}(2*hh-1+(k-1)*10,:),...
                     PCA_B_heading_projPC{which_two_dimension(2)}(2*hh-1+(k-1)*10,:),'-','color',colors_angles{k}(hh,:),'linew',3);
                % Null
                plot(PCA_B_heading_projPC{which_two_dimension(1)}(2*hh+(k-1)*10,:),...
                     PCA_B_heading_projPC{which_two_dimension(2)}(2*hh+(k-1)*10,:),'--','color',colors_angles{k}(hh,:),'linew',3);
            end
            
            axis tight;  grid off;
            axis off
        end
        
        
        % ============  1. Variance explained =================
        set(figure(2099+figN),'pos',[936 491 733 463],'name',['Variance explained, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        plot((1:length(PCA_B_explained_heading))', cumsum(PCA_B_explained_heading),'o-','markersize',8,'linew',1.5);
        plot((1:denoised_dim)',cumsum(PCA_B_explained_heading(1:denoised_dim)),'ro-','markersize',8,'linew',1.5,'markerfacecol',colors(2,:));
        plot([0 1],[0 PCA_B_explained_heading(1)],'r-','linew',1.5);
        plot(xlim,[1 1]*sum(PCA_B_explained_heading(1:denoised_dim)),'r--');
        text(denoised_dim,sum(PCA_B_explained_heading(1:denoised_dim))*0.9,[num2str(sum(PCA_B_explained_heading(1:denoised_dim))) '%'],'color',colors(2,:));
        SetFigure(); xlabel('Num of principal components'); ylabel('Explained variability (%)'); ylim([0 100]);
        
        % ============  2. 1-D Trajectory of Eigen-neurons =================
        
        set(figure(2099+figN),'pos',[462 67 1207 889],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        ds = [1:9];
        
        for d = 1:length(ds)
            ranges(d) = range(PCA_B_heading_projPC{ds(d)}(:));
        end
        
        % Plotting
        for d = 1:length(ds)
            % Normalize to [-1,1]
            gain = max(ranges) / (2 * 0.8) ;
            offset = mean(PCA_B_heading_projPC{ds(d)}(:));
            norm_proj_PC_this = (PCA_B_heading_projPC{ds(d)}-offset)/gain;
            
            h = subplot(fix(sqrt(length(ds))),ceil(length(ds)/fix(sqrt(length(ds)))),d);
            
            colors_angles_all = [colors_angles{1};colors_angles{2};colors_angles{3}];
            colors_angles_all = reshape(repmat(colors_angles_all,1,2)',3,[])';
            colors_angles_all = mat2cell(colors_angles_all,ones(30,1));
            
            SeriesComparison(shiftdim(norm_proj_PC_this',-1),PCA_B_times,'Colors',...
               colors_angles_all,'LineStyles',{'-','--'},'axes',h);
            axis tight; ylim([-1 1]); xlim([-200 1860])
            for tt = 1:3
                plot([1 1] * time_markers{j_PCA_B}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_PCA_B}{tt},'linew',1.5);
            end
            set(gca,'ytick',[-1 0 1]);
            xlabel([]); ylabel('Amplitude (a.u.)'); title(['Eigen-neuron ' num2str(ds(d))]); legend off;
        end
        % set(get(gcf,'children'),'ylim',[min_ylim max_ylim]);
        SetFigure(15);
        
    end

    %% =========== Demixed PCA, preparation ================  HH20190628
        j_dPCA = 2; % Align to onset
        
        unique_abs_heading = unique(abs(group_result(representative_cell).unique_heading)); 

        dPCA1_PSTHAver_NHCMT = [];
        dPCA1_PSTHTrial_NHCMTK = [];
        dPCA1_TrialN_NHCM = [];
        
        dPCA2_PSTHAver_NMCT = [];
        dPCA2_PSTHTrial_NMCTK = [];
        dPCA2_TrialN_NMC = [];
        
    
    function f55p0()     % Demixed PCA, preparation

        time_length = length(rate_ts{j_dPCA});
        
        max_trialNum = 90; % This is known post hoc (the actual number is 87)
        dPCA1_PSTHAver_NHCMT = NaN(N,length(unique_abs_heading),2,3,time_length);
        dPCA1_PSTHTrial_NHCMTK = NaN(N,length(unique_abs_heading),2,3,time_length,max_trialNum);
        dPCA1_TrialN_NHCM = NaN(N,length(unique_abs_heading),2,3);
        
        max_trialNum = 250; % This is known post hoc (the actual number is 233)
        dPCA2_PSTHAver_NMCT = NaN(N,3,2,time_length);
        dPCA2_PSTHTrial_NMCTK = NaN(N,3,2,time_length,max_trialNum);
        dPCA2_TrialN_NMC = NaN(N,3,2);
        
       %% Organize my data
        
        
        % Should put in Batch file
        for nn = 1:N
            this_unique_abs_heading = unique(abs(group_result(nn).unique_heading)); % Group abs(headings) to the same levels
            this_unique_abs_heading_add0 = union(this_unique_abs_heading,0); % Deal with sessions without 0 heading

            this_raw = group_result(nn).mat_raw_PSTH;
            this_PSTH_per_trial = this_raw.spike_hist{j_dPCA};
            this_stim_type_per_trial = this_raw.stim_type_per_trial';
            this_heading_per_trial = this_raw.heading_per_trial';
            this_choice_per_trial = this_raw.choice_per_trial;
            this_correctOr0_per_trial = ((this_choice_per_trial > 1.5) == (this_heading_per_trial >= 0)) | ...
                                        (this_heading_per_trial == 0); % Correct or 0 heading trials
            
            for stim_type = 1:3
               kk = find(stim_type == group_result(i).mat_raw_PSTH.unique_stim_type);
               if ~isempty(kk)   % We have this condition
                   for cc = 1:2  % 1: pref; 2: Null
                       
                    %% Different ways of dealing with headings

                       % 1) f55p1: PSTH: [N, abs(Heading), choice, modality], for each modality separately

                       for hh = 1:length(unique_abs_heading)
                           this_selected = (this_stim_type_per_trial == kk) & ...
                               (abs(this_heading_per_trial) == this_unique_abs_heading_add0(hh)) & ...
                               (this_choice_per_trial == (cc == 1) * this_raw.PREF + (cc == 2) * (3 - this_raw.PREF)) & ... % Align to [Pref, Null]
                               (this_correctOr0_per_trial); % Correct only
                           
                           dPCA1_PSTHTrial_NHCMTK(nn, hh, cc, stim_type, :, 1:sum(this_selected)) = this_PSTH_per_trial(this_selected,:)';
                           dPCA1_PSTHAver_NHCMT(nn, hh, cc, stim_type, :) = mean(this_PSTH_per_trial(this_selected,:),1);
                           dPCA1_TrialN_NHCM (nn, hh, cc, stim_type) = sum(this_selected);
                       end
                       
                       % 2) f55p2: PSTH: [N, modality, choice], different headings grouped together, exactly what the 2nd reviewer wants, like TDR
                       % Note that here 'choice' is still in the third dimension to keep compatible with the dPCA_plot function.
                       this_selected = (this_stim_type_per_trial == kk) & ...
                           (this_choice_per_trial == (cc == 1) * this_raw.PREF + (cc == 2) * (3 - this_raw.PREF)) & ...  % Align to [Pref, Null]
                           (this_correctOr0_per_trial); % Correct only
                       
                       dPCA2_PSTHTrial_NMCTK(nn, stim_type, cc, :, 1:sum(this_selected)) = this_PSTH_per_trial(this_selected,:)';
                       dPCA2_PSTHAver_NMCT(nn, stim_type, cc, :) = mean(this_PSTH_per_trial(this_selected,:),1);
                       dPCA2_TrialN_NMC (nn, stim_type, cc) = sum(this_selected);

                       
                   end
               end
            end
        end

    end

   %% 1. dPCA over heading and choice, separately for each stim_type
    % Use abs(heading), correct only. Because lots of wrong trials are missing (at the largest heading).
    % In fact, if wrong trials are included, there will be only less than 20 cells that have data for ALL parameter combinations.
    % But in the meantime, this will make choice and heading highly correlated, so the task-difficulty dependence trace is still 
    % in the choice's second component.
    
    function f55p1(debug)     % Demixed PCA. %HH20190625
        if debug;  dbstack;  keyboard;  end
        
        if isempty(dPCA1_PSTHAver_NHCMT)
            f55p0();
        end
        
        set(0,'defaultAxesColorOrder','remove');

        for kk = 1:3
            
            min_trials = 2; % Follows Kobak's paper
            
            % Select cells that have at least min_trials trials for ALL parameter combinations
            this_full_data = all(all(dPCA1_TrialN_NHCM(:,:,:,kk) >= min_trials,2),3); 
            
            firingRates = squeeze(dPCA1_PSTHTrial_NHCMTK(this_full_data,:,:,kk,:,:));
            firingRatesAverage = squeeze(dPCA1_PSTHAver_NHCMT(this_full_data,:,:,kk,:));
            trialNum = dPCA1_TrialN_NHCM(this_full_data,:,:,kk);
            ifSimultaneousRecording = 0;
            
            combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}; 
            margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
            margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
            time = rate_ts{j_dPCA};
            timeEvents = time_markers{j_dPCA}(1,:);
           
            % -- Adapted from Kobak et al. 2016 --
            %{
            
            %% --- dPCA without regularization and ignoring noise covariance ---

            % This is the core function.
            % W is the decoder, V is the encoder (ordered by explained variance),
            % whichMarg is an array that tells you which component comes from which
            % marginalization
            
            tic
            [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
                'combinedParams', combinedParams);
            toc
            
            explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
                'combinedParams', combinedParams);
            
            dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
                'explainedVar', explVar, ...
                'marginalizationNames', margNames, ...
                'marginalizationColours', margColours, ...
                'whichMarg', whichMarg,                 ...
                'time', time,                        ...
                'timeEvents', timeEvents,               ...
                'timeMarginalization', 3, ...
                'legendSubplot', 16);
            
            set(gcf,'Name',sprintf('Without Regulation, stim_type = %g, cell_number = %g',kk, sum(this_full_data)));
           %}
            
           %% --- dPCA with regularization and Noise (Co)variance ---
%             %{
            % This function takes some minutes to run. It will save the computations
            % in a .mat file with a given name. Once computed, you can simply load
            % lambdas out of this file:
            %   load('tmp_optimalLambdas.mat', 'optimalLambda')
            
            % Please note that this now includes noise covariance matrix Cnoise which
            % tends to provide substantial regularization by itself (even with lambda set
            % to zero).
            
            optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
                'combinedParams', combinedParams, ...
                'simultaneous', ifSimultaneousRecording, ...
                'numRep', 2, ...  % increase this number to ~10 for better accuracy
                'filename', 'tmp_optimalLambdas.mat');
            
            Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
                firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);
            
            [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
                'combinedParams', combinedParams, ...
                'lambda', optimalLambda, ...
                'Cnoise', Cnoise);
            
            explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
                'combinedParams', combinedParams);
            
            dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
                'explainedVar', explVar, ...
                'marginalizationNames', margNames, ...
                'marginalizationColours', margColours, ...
                'whichMarg', whichMarg,                 ...
                'time', time,                        ...
                'timeEvents', timeEvents,               ...
                'timeMarginalization', 3,           ...
                'legendSubplot', 16);
            
            set(gcf,'Name',sprintf('With Regulation, stim_type = %g, cell_number = %g',kk, sum(this_full_data)));

            %}
            
    
        end
        
        set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255); % Restore color order

    end

    %% 2. dPCA over choice and modality (exactly what the 2nd reviewer asked)
     % Also ignore the wrong trials. Everything is like the TDR analysis in the original manuscript
    
    function f55p2(debug)     % Demixed PCA. %HH20190625
        if debug;  dbstack;  keyboard;  end
        
        if isempty(dPCA1_PSTHAver_NHCMT)
            f55p0();
        end
        
        set(0,'defaultAxesColorOrder','remove');
                
        
        % --- Select cells that have at least min_trials trials for ALL parameter combinations ---
        min_trials = 2; % Follows Kobak's paper
        this_full_data = all(all(dPCA2_TrialN_NMC(:,:,:) >= min_trials,2),3);
        
        firingRates = dPCA2_PSTHTrial_NMCTK(this_full_data,:,:,:,:);
        firingRatesAverage = dPCA2_PSTHAver_NMCT(this_full_data,:,:,:);
        trialNum = dPCA2_TrialN_NMC(this_full_data,:,:);
        ifSimultaneousRecording = 0;
        
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'Modality', 'Decision', 'Condition-independent', 'D/M Interaction'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
        time = rate_ts{j_dPCA};
        timeEvents = time_markers{j_dPCA}(1,:);
        
        % -- Adapted from Kobak et al. 2016 --
    
        %% --- dPCA with regularization and Noise (Co)variance ---
        % This function takes some minutes to run. It will save the computations
        % in a .mat file with a given name. Once computed, you can simply load
        % lambdas out of this file:
        %   load('tmp_optimalLambdas.mat', 'optimalLambda')
        
        % Please note that this now includes noise covariance matrix Cnoise which
        % tends to provide substantial regularization by itself (even with lambda set
        % to zero).
        
        optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
            'combinedParams', combinedParams, ...
            'simultaneous', ifSimultaneousRecording, ...
            'numRep', 10, ...  % increase this number to ~10 for better accuracy
            'filename', 'tmp_optimalLambdas.mat');
        
        Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
            firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);
        
        [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
            'combinedParams', combinedParams, ...
            'lambda', optimalLambda, ...
            'Cnoise', Cnoise);
        
        explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
            'combinedParams', combinedParams);
        
%         dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%             'explainedVar', explVar, ...
%             'marginalizationNames', margNames, ...
%             'marginalizationColours', margColours, ...
%             'whichMarg', whichMarg,                 ...
%             'time', time,                        ...
%             'timeEvents', timeEvents,               ...
%             'timeMarginalization', 3,           ...
%             'legendSubplot', 16);
        
        
        % === dPCA with Decoding ===
        
        decodingClasses = {[(1:3)' (1:3)'], repmat([1:2], [3 1]), [], [(1:3)' (4:6)']};
        
        accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
            'lambda', optimalLambda, ...
            'combinedParams', combinedParams, ...
            'decodingClasses', decodingClasses, ...
            'simultaneous', ifSimultaneousRecording, ...
            'numRep', 20, ...        % increase to 100
            'filename', 'tmp_classification_accuracy.mat');
        
        dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
        
        accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
            'lambda', optimalLambda, ...
            'combinedParams', combinedParams, ...
            'decodingClasses', decodingClasses, ...
            'simultaneous', ifSimultaneousRecording, ...
            'numRep', 20, ...        % increase to 100
            'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
            'filename', 'tmp_classification_accuracy.mat');
        
        dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
        
        componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
        
        dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'marginalizationNames', margNames, ...
            'marginalizationColours', margColours, ...
            'whichMarg', whichMarg,                 ...
            'time', time,                        ...
            'timeEvents', timeEvents,               ...
            'timeMarginalization', 3,           ...
            'legendSubplot', 16,                ...
            'componentsSignif', componentsSignif);
        
        set(gcf,'Name',sprintf('With Regulation, Choice amd Modality, cell_number = %g', sum(this_full_data)));

        set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255); % Restore color order


    end

    %% 3. dPCA over choice and modality and heading
     % Also ignore the wrong trials. Everything is like the TDR analysis in the original manuscript
    
    function f55p3(debug)     % Demixed PCA. %HH20190711
        if debug;  dbstack;  keyboard;  end
        
        if isempty(dPCA1_PSTHAver_NHCMT)
            f55p0();
        end
        
        set(0,'defaultAxesColorOrder','remove');
        
        % --- Select cells that have at least min_trials trials for ALL parameter combinations ---
        min_trials = 2; % Follows Kobak's paper
        this_full_data = all(all(dPCA1_TrialN_NHCM(:,:,:) >= min_trials,2),3);
        
        firingRates = dPCA1_PSTHTrial_NHCMTK(this_full_data,:,:,:,:,:);
        firingRatesAverage = dPCA1_PSTHAver_NHCMT(this_full_data,:,:,:,:);
        trialNum = dPCA1_TrialN_NHCM(this_full_data,:,:,:);
        ifSimultaneousRecording = 0;
        
        % 1: abs(Heading), 2: Choice, 3: Modality, 4: Time
        combinedParams = {{1, [1 4]}, {2, [2 4]}, {3, [3 4]}, {4}};  % No interaction. Just like Rossi-Pool
        margNames = {'Heading', 'Decision', 'Modality', 'Time'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
        time = rate_ts{j_dPCA};
        timeEvents = time_markers{j_dPCA}(1,:);
        
        % -- Adapted from Kobak et al. 2016 --
    
        %% --- dPCA with regularization and Noise (Co)variance ---
        % This function takes some minutes to run. It will save the computations
        % in a .mat file with a given name. Once computed, you can simply load
        % lambdas out of this file:
        %   load('tmp_optimalLambdas.mat', 'optimalLambda')
        
        % Please note that this now includes noise covariance matrix Cnoise which
        % tends to provide substantial regularization by itself (even with lambda set
        % to zero).
        
        optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
            'combinedParams', combinedParams, ...
            'simultaneous', ifSimultaneousRecording, ...
            'numRep', 2, ...  % increase this number to ~10 for better accuracy
            'filename', 'tmp_optimalLambdas.mat');
        
        Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
            firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);
        
        [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
            'combinedParams', combinedParams, ...
            'lambda', optimalLambda, ...
            'Cnoise', Cnoise);
        
        explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
            'combinedParams', combinedParams);
        
        dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'marginalizationNames', margNames, ...
            'marginalizationColours', margColours, ...
            'whichMarg', whichMarg,                 ...
            'time', time,                        ...
            'timeEvents', timeEvents,               ...
            'timeMarginalization', 4,           ...  % Note here
            'legendSubplot', 16, ...
            'showNonsignificantComponents', 1);
        
        set(gcf,'Name',sprintf('With Regulation, Choice amd Modality, cell_number = %g', sum(this_full_data)));
                
        set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255); % Restore color order

    end

    function f55p4(debug)     % Significant cells using choice divergence / PSTH p test (2nd reviewer)
        if debug;  dbstack;  keyboard;  end
        
        for jj = 1:2
            CD_p_all{jj} = nan(N,length(rate_ts{jj}),3);
            
            for nn = 1:N
%                 CD_p_all{j}(nn,:,:) = group_result(nn).mat_raw_PSTH.ChoiceDivergence_ALL_perm{jj}.p';  % From choice divergence, permutation test (a caveat
                                                                                                         % is that this will overestimate the significance, because
                                                                                                         % of limitation in the number of permutations)
                CD_p_all{jj}(nn,:,:) = group_result(nn).mat_raw_PSTH.PSTH{jj,1,1}.ps';  % j = 1  % From PSTH, t-test
            end
            
            percent_sign{jj} = squeeze(sum(CD_p_all{jj}(:,:,:) < p_critical,1) / N * 100);
            
        end
        
        figure(155010); set(gcf,'uni','norm','pos',[0.023       0.127       0.936       0.549]); clf
        
        if isempty(firstDivergenceTime)
            f1p1p6p1(0); % Calculate divergence time
        end
        
        for k = 1:3
            
            % [divergence_time_ordered, plot_order] = sort(firstDivergenceTime(:,k)); % Not good because some cell does not have divergence time
            [~, plot_order] =sort(abs(group_ChoicePreference(3,:,3)),'descend');
%             [~, plot_order] =sort(group_ChoicePreference_pvalue(k,:,3),'ascend');
            p_ordered = group_ChoicePreference_pvalue(k,plot_order,3);

            subplot(1,5,k); hold on;
            imagesc('XData', rate_ts{1}, 'CData', CD_p_all{1}(plot_order,:,k) < p_critical); 
            
            sign_ind = find(p_ordered < p_critical);
            plot(0*sign_ind,sign_ind,'rs','linew',1,'markerfacecol','r');
            
            axis tight;
            
%             for nn = 1:N
%                  plot(divergence_time_ordered(nn),nn,'ro');
%             end
            
            % Time markers
            for tt = 1:3
                plot([1 1] * time_markers{j_for_SVM}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_for_SVM}{tt},'linew',1);
            end
            plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/3 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
            
        end
        
        h_subplot = subplot(1,5,[4 5]);
        SeriesComparison({shiftdim(percent_sign{1}, -1) shiftdim(percent_sign{2}, -1)},...
            {rate_ts{1} rate_ts{2} time_markers},...
            'Colors',{colors(1,:),colors(2,:),colors(3,:)},'LineStyles',{'-'},...
            'Xlabel',[],'Ylabel','Fraction of cells (%)','axes',h_subplot);
        plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/3 + p_critical*100,'--','linew',2.5,'color',[0.6 0.6 0.6]);
        
        plot(xlim,[p_critical p_critical]*100,'k--'); % By definition, false alarm level (baseline) = p_critical
        legend off; SetFigure();

    end

    ANOVAN_1_p_all = [];
    ANOVAN_2_p_all = [];
    ANOVAN_ts = [];
    
    function f55p5(debug)     % Significant cells using ANOVA (2nd reviewer)  HH20190701
        if debug;  dbstack;  keyboard;  end
        
        ANOVAN_p_tRange = {[-300 1800],[-300 300]};
        
        if isempty(ANOVAN_1_p_all)
            for jj = 1:2
                
                this_t_ind = find(ANOVAN_p_tRange{jj}(1) < rate_ts{jj} & rate_ts{jj} < ANOVAN_p_tRange{jj}(2));
                this_t_ind = this_t_ind(1:2:end); % Downsample a little bit
                
                ANOVAN_ts{jj} = rate_ts{jj}(this_t_ind);
                ANOVAN_1_p_all{jj} = nan(N,length(ANOVAN_ts{jj}),3);
                ANOVAN_2_p_all{jj} = nan(N,length(ANOVAN_ts{jj}),2,3);
                
                %                 parfor_progress(N);
                progressbar(['j = ' num2str(jj) ', cell num']);
                
                for nn = 1:N
                    
                    this_PSTH = group_result(nn).mat_raw_PSTH.spike_hist{jj}(:,this_t_ind);
                    this_stim_type = group_result(nn).mat_raw_PSTH.stim_type_per_trial;
                    this_heading = group_result(nn).mat_raw_PSTH.heading_per_trial;
                    this_choice = group_result(nn).mat_raw_PSTH.choice_per_trial';
                    
                    for ttt = 1 : length(this_t_ind)
                        % --- ANOVA_1: on heading, stim_type, and choice ---
                        [this_p, ~, ~] = anovan(this_PSTH(:,ttt),{this_heading, this_stim_type, this_choice},'continuous',1,'display','off');
                        ANOVAN_1_p_all{jj}(nn,ttt,:) = this_p';
                        
                        % --- ANOVA_2: on heading and choice, for each stim_type separately ---
                        for kk = 1:3
                            select_kk = this_stim_type == kk;
                            [this_p, ~, ~] = anovan(this_PSTH(select_kk,ttt),{this_heading(select_kk), this_choice(select_kk)},'continuous',1,'display','off');
                            ANOVAN_2_p_all{jj}(nn,ttt,:,kk) = this_p';
                        end
                    end
                    
                    progressbar(nn/N)
                    %                     parfor_progress();
                end
                %                 parfor_progress(0);
                
            end
        end
        
        for jj = 1:2
            percent_significant_1{jj} = squeeze(sum(ANOVAN_1_p_all{jj}(:,:,:) < p_critical,1) / N * 100);
            percent_significant_2{jj} = squeeze(sum(ANOVAN_2_p_all{jj}(:,:,:,:) < p_critical,1) / N * 100);
        end
        
        figure(191948); set(gcf,'Name','ANCOVA: heading, stim_type, choice','uni','norm','pos',[ 0.036905         0.1     0.93869     0.37619]); clf
        
        h_subplot = subplot(1,3,1);
        SeriesComparison({shiftdim(percent_significant_1{1}, -1) shiftdim(percent_significant_1{2}, -1)},...
            {ANOVAN_ts{1} ANOVAN_ts{2} time_markers},...
            'Colors',{'k','c','m'},'LineStyles',{'-'},...
            'Xlabel',[],'Ylabel','Fraction of cells (%)','figN',191948,'axes',h_subplot);
        legend off; ylim([0 90]); title('Heading, Choice, Modality')
        plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/3 + p_critical*100,'--','linew',2.5,'color',[0.6 0.6 0.6]);
        
        plot(xlim,[p_critical p_critical]*100,'k--'); % By definition, false alarm level (baseline) = p_critical
        
        for hc = 1:2 % Heading, choice
            h_subplot = subplot(1,3,1+hc);
            SeriesComparison({shiftdim(squeeze(percent_significant_2{1}(:,hc,:)), -1) shiftdim(squeeze(percent_significant_2{2}(:,hc,:)), -1)},...
                {ANOVAN_ts{1} ANOVAN_ts{2} time_markers},...
                'Colors',{colors(1,:),colors(2,:),colors(3,:)},'LineStyles',{'-'},...
                'Xlabel',[],'Ylabel','Fraction of cells (%)','figN',191948,'axes',h_subplot);
            legend off;            ylim([0 90])
            
            if hc == 1; title('Heading'); else title('Choice'); end
            
            plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/3 + p_critical*100,'--','linew',2.5,'color',[0.6 0.6 0.6]);
            plot(xlim,[p_critical p_critical]*100,'k--'); % By definition, false alarm level (baseline) = p_critical
        end
        
        SetFigure();
        
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
%         SVM_training_epoch = [1500 1700];    % This make it comparable with PCA results

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
        if ver_num < 2014
            if matlabpool('size') == 0 ;  matlabpool;  end
        else
            if isempty(gcp('nocreate')); parpool; end
        end
        
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
            svm_training_choice(nn) = svmtrain(firing_for_this_training,answers_choice,'box',1e-5,'tol',1e-7,'autoscale',0);
            svm_training_modality(nn) = svmtrain(firing_for_this_training,answers_modality,'box',1e-5,'tol',1e-7,'autoscale',0);
            
            parfor_progress;
        end
        parfor_progress(0);
        toc; disp('Done!');
        
        %% Averaged weights (bagging)
%         weights_svm_choice_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_choice,'uniform',0));
%         weights_svm_modality_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_modality,'uniform',0));

        % By default, 'autoscale' = true, so here I should rescale back the weight! HH20170613
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
            weights_svm_choice_mean(sortt(who_are_tcells(sortt))),0.5),'facecolor',colors(2,:),'edgecolor','none');
        
        xlim([-2 sum(select_for_SVM_actual)+1]); ylim([-1 1]); title('Choice');
        set(gca,'xtick',[1 sum(select_for_SVM_actual)]);
        % errorbar(svm_weights_choice_mean(sortt),svm_weights_choice_sem(sortt),'.');
        
        % Weights for modality decoder
        subplot(2,2,2);
        [~,sortt] = sort(abs(weights_svm_modality_mean),'descend');
        set(bar(weights_svm_modality_mean(sortt),0.5),'facecolor','k','edgecolor','none'); hold on;
        set(bar(find(who_are_tcells(sortt)),...
            weights_svm_modality_mean(sortt(who_are_tcells(sortt))),0.5),'facecolor',colors(2,:),'edgecolor','none');
        
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
            'SameScale',1,'Method','Pearson','FittingMethod',2);
        delete([h.diag h.group(1:2).line]); legend off;
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
        text(min(xlim)+0.2,min(ylim)+0.2,sprintf('r^2 = %g, p = %g',h.group(3).r_square,h.group(3).p),'fontsize',11);
        
        % Annotate tcells
        h_line = plot(weights_svm_modality_mean(who_are_tcells),...
            weights_svm_choice_mean(who_are_tcells),...
            '+','markersize',8,'color',colors(2,:),'linew',2);
        
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
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2);
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),weights_svm_choice_mean(select_tcells(select_for_SVM),1),...
            '+','markersize',16,'color',colors(2,:),'linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(MemSac_indicator(select_for_SVM),(weights_svm_choice_mean(:,1)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_SVM});
        
        
        % ---------  Weights for eigen-neuron 1 and weights for SVM choice decoder
        if isempty(weights_PCA_B_PC)
            f5p0(0);
        end
        
        set(figure(613),'pos',[89 26 761 606]); clf
                
        h = LinearCorrelation({
            Modality_pref_all(1,:,3);
            % weights_PCA_B_PC(:,1);
            },...
            {
            % (weights_svm_choice_mean);
            (weights_svm_modality_mean);
            },...
            'CombinedIndex',[],... % 'Xlabel','Weight for eigen-neuron 1','Ylabel','Weight in choice decoder',...
             'Xlabel','Modality_pref_all (vis - vest)','Ylabel','Weight in modality decoder',...
            'FaceColors',{'k'},'Markers',{'o'},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',613,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman','FittingMethod',2);
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(weights_PCA_B_PC(select_tcells(select_for_SVM),1),weights_svm_choice_mean(select_tcells(select_for_SVM),1),...
            '+','markersize',16,'color',colors(2,:),'linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        set([gca h.group.dots],'ButtonDownFcn',{@Show_individual_cell, h.group.dots, select_for_SVM});
        
        % ------------ Weights for SVM choice decoder VS choice divergence and multisensory enhancement 20170609 ------------------
        Choice_pref_all_comb = abs(Choice_pref_all(3,:,3)); % stim-on to stim-off
        Choice_pref_all_vest = abs(Choice_pref_all(1,:,3));
        Choice_pref_all_vis = abs(Choice_pref_all(2,:,3));
        
        figure(1234); clf;
        hs = tight_subplot(2,2,0.1,0.1,0.1);
        
        % -- SVM weights vs Div_vest --
        xxs = {Choice_pref_all_vest(:), 'Divergence vest';
            Choice_pref_all_vis(:),  'Divergence vis';
            Choice_pref_all_comb(:), 'Divergence comb';
           % Choice_pref_all_comb(:)./(Choice_pref_all_vest(:)+Choice_pref_all_vis(:)), '(comb-vest)+(comb-vis)'};
           % 2*Choice_pref_all_comb(:)-(Choice_pref_all_vest(:)+Choice_pref_all_vis(:)), '(comb-vest)+(comb-vis)'};
            Choice_pref_all_comb(:)-(Choice_pref_all_vest(:)), '(comb-vest)+(comb-vis)'};
   
        for xxxx = 1:length(xxs)
            xx = xxs{xxxx,1};
            yy = (weights_svm_choice_mean);
            hl = LinearCorrelation(xx,yy,'Axes',hs(xxxx),'MethodOfCorr','Spearman','FittingMethod',2);
            hold on;
            delete([hl.leg]);
            
            xlabel(xxs{xxxx,2}); ylabel('abs(svm weight)');
            text(min(xlim),min(ylim)+range(ylim)*0.9,...
                sprintf('r^2=%g\n p=%g\n k = %g\\pm%g',hl.group.r_square,hl.group.p,hl.group.para(1),hl.group.paraSE(1)))
            
            % Annotate tcells
            plot(xx(select_tcells),yy(select_tcells(select_for_SVM),1),...
                '+','markersize',16,'color',colors(2,:),'linew',2);

        end

        
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
        shadedErrorBar(SVM_testing_tcenters,nanmean(corr_rate_choice_by_choice),...
            nanstd(corr_rate_choice_by_choice),'lineprops',{'color','m','linew',2}); hold on;
        shadedErrorBar(SVM_testing_tcenters,nanmean(corr_rate_modality_by_choice),...
            nanstd(corr_rate_modality_by_choice),'lineprops',{'color',[.87 .49 0],'linew',2});
        legend('Decodes choice','Decodes modality','location','best');
        
        axis tight;  ylim([0.35 1.05]); 
        plot(xlim,[0.5 0.5],'k--'); title('Choice decoder');
        
        % Time markers
        for tt = 1:3
            plot([1 1] * time_markers{j_for_SVM}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_for_SVM}{tt},'linew',1);
        end
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_for_SVM}(1), .5+ Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
        subplot(1,2,2);
        shadedErrorBar(SVM_testing_tcenters,nanmean(corr_rate_modality_by_modality),...
            nanstd(corr_rate_modality_by_modality),'lineprops',{'color',[.87 .49 0],'linew',2}); hold on;
        shadedErrorBar(SVM_testing_tcenters,nanmean(corr_rate_choice_by_modality),...
            nanstd(corr_rate_choice_by_modality),'lineprops',{'color','m','linew',2});
        
        axis tight; ylim([0.35 1.05]); 
        plot(xlim,[0.5 0.5],'k--'); title('Modality decoder');
        
        % Time markers
        for tt = 1:3
            plot([1 1] * time_markers{j_for_SVM}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_for_SVM}{tt},'linew',1);
        end
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_for_SVM}(1),.5+ Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
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
        
        % projs = {weights_svm_choice_mean,weights_svm_modality_mean};
        projs = {weights_svm_choice_allbootstrap,weights_svm_modality_allbootstrap}; % With bootstrap
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

%         t_span = 200;    
        t_span = 1000;
        t_step = 50; % For time evolution
        toi = [950 1500] +  time_markers{j_for_reg}(1,1);   % Time of interests for weights illustration
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        find_for_regression = find(select_for_regression);
        t_centers = ts(1) + t_span / 2 : t_step : ts(end)-t_span / 2;
        
        weight_vest_vis = nan(length(find_for_regression),2 * 2 + 2,length(t_centers)); % Weights (3), p_values (3), r^2, p-value of r^2
        measure_pred_all = nan(length(find_for_regression) * 10,2,length(t_centers));
        
%         weight_vest_vis_lsqlin = nan(length(find_for_regression),2,length(t_centers));
        
        progressbar(['Number of cells (all ' num2str(length(find_for_regression)) ')']);
        
        for i = 1:length(find_for_regression)  % For each cell
            
            for tcc = 1:length(t_centers) % For each time bin
                
                t_range = t_centers(tcc) - t_span / 2 <= ts & ts <= t_centers(tcc) + t_span / 2;
                
                r = [];
                for k = 1:3
                    r(:,k) = mean(group_result(find_for_regression(i)).mat_raw_PSTH.PSTH{j_for_reg,2,k}.ys(:,t_range),2);
                end
                
                
                % GLM fit
                r = r - repmat(nanmean(r,1),size(r,1),1); % Mean removed
                r(any(isnan(r),2),:) = [];
                [b,~,stat] = glmfit(r(:,1:2) ,r(:,3),'normal','link','identity','constant','off');

%                 r(any(isnan(r),2),:) = [];
%                 [b,~,stat] = glmfit(r(:,1:2) ,r(:,3),'normal','link','identity','constant','on');
%                 b = b(2:3); stat.p = stat.p(2:3);
                
                % Override by ridge regression
                X = r(:,1:2); Y = r(:,3);
                lamda = 0;
                b = inv(X'*X + lamda*eye(size(X,2))) * X' * Y;

                [~,~,~,~,stat_reg] = regress(r(:,3),[ones(size(r,1),1) r(:,1:2)]);
                weight_vest_vis(i,:,tcc) = [b' stat.p' stat_reg([1 3])]; % Weights, r^2 and p-value

%                 b_lsqlin = lsqlin([r(:,1:2)],r(:,3),[],[],[],[],[0 0]');
%                 b_lsqlin = lsqlin([ones(size(r(:,1))) r(:,1:2)],r(:,3),[],[],[],[],[0 0 0]'); b_lsqlin = b_lsqlin(2:3);

%                 weight_vest_vis_lsqlin(i,:,tcc) = b_lsqlin;
                
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
        
        temp_col = {colors(1,:),colors(2,:),'k'};
        temp_marker = {'s','o','o'}; % {'','',''};
        for pp = 1 : 3
            h = shadedErrorBar(repmat(t_centers',1,1),mean_paras(:,pp),sem_paras(:,pp),...
                'lineprops',{[temp_marker{pp} temp_col{pp} '-'],'linew',2,'markersize',10},'transparent',transparent);
            hold on;
        end
        
        plot(repmat(t_centers',1,1),mean_paras(:,1)+mean_paras(:,2),'color',colors(3,:),'linew',2);
        plot(xlim,[1 1],'g:','linew',3);
        
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
            
            % Predicted v.s. 
            ax0 = subplot(2,3,1 + (toii-1)*3);
            h = LinearCorrelation( measure_pred_all(:,1,toi_ind), measure_pred_all(:,2,toi_ind),...
                'CombinedIndex',[],...
                'Xlabel','Measured response','Ylabel','Predicted response',...
                'FaceColors',{[0.9 0.9 0.9]},'Markers',{'.'},...
                'LineStyles',{'k-'},'MarkerSize',10,...
                'XHist',0,'YHist',0,...
                'XHistStyle','grouped','YHistStyle','grouped','SameScale',1,'Method','Pearson','FittingMethod',2,'axes',ax0);
            text(min(xlim)+2,min(ylim)+2,sprintf('r^2 = %g, p = %g',h.group(1).r_square,h.group(1).p),'FontSize',13);
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
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,'Method','Pearson','FittingMethod',2,'axes',ax2);
            
            if toii == 1 ; xlabel(''); end
            delete([h.group(1).line h.diag]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--'); plot(xlim,[1 1],'k--'); plot([1 1],ylim,'k--');
            text(min(xlim)+0.2,min(ylim)+0.2,sprintf('r^2 = %g, p = %g',h.group(2).r_square,h.group(2).p),'FontSize',13);
            legend off;
            set(gca,'xtick',-10:1:10,'ytick',-10:1:10);
            title(sprintf('Ridge \\lamda = %g', lamda))

            %{
            % Annotate tcells
            h_t = plot(weight_vest_vis(select_tcells(select_for_regression),1,toi_ind),...
                weight_vest_vis(select_tcells(select_for_regression),2,toi_ind),'+','markersize',10,'color',colors(2,:),'linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(weight_vest_vis(:,1,toi_ind), weight_vest_vis(:,2,toi_ind),'visible','off'); hold on;
            set([gca [h.group.dots] h_line h_t],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_regression});
            %}
            
            %{
            ax3 = subplot(2,4,4 + (toii-1)*4);
            plot(ax3,weight_vest_vis_lsqlin(:,1,toi_ind),weight_vest_vis_lsqlin(:,2,toi_ind),'o'); hold on
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--'); plot(xlim,[1 1],'k--'); plot([1 1],ylim,'k--');
            set(gca,'xtick',-10:1:10,'ytick',-10:1:10);
            %}
        end
        
        SetFigure(18);
    end

    %% Regression 2. Fitting model traces with real data. (We decided to do this when Alex came to Shanghai) @HH20170808
    %%%%%%%%%%%%%%%%%%% for f7p2: load model traces %%%%%%%%%%%%%%%
    % Using "result = lip_HH({},{'mean_diff_PSTH_correct_allheading','diff_PSTH_correct_mean_allheading','ts'});"
    
    % ==== Load model data (Fit target) from the trace without heterogeneity ====
    % Generated from calling "result = lip_HH({},{'ts','mean_diff_PSTH_correct_allheading'});"
    model_mean_PSTH_trace_without_heter = load('mean_trace_without_heter.mat'); % Gamma = 0
    % model_mean_PSTH_trace_without_heter = load('mean__trace_without_heter_gamma=1.mat'); % Gamma = 1
    % model_mean_PSTH_trace_without_heter = load('mean__trace_without_heter_gamma=8.mat'); % Gamma = 8
%     model_mean_PSTH_trace_without_heter = load('diff_PSTH_trace_without_heter_vestvelocity.mat'); % Vestibular velocity

    model_ts = model_mean_PSTH_trace_without_heter.result.ts*1000;
    model_PSTH_optimal_M1 = squeeze(model_mean_PSTH_trace_without_heter.result.mean_diff_PSTH_correct_allheading);
    
    % ==== Single cells from the model ====
    model_PSTH_trace_with_heter_optimal = load('diff_PSTH_trace_with_heter.mat');
    model_PSTH_trace_with_heter_vest10_vis1 = load('diff_PSTH_trace_with_heter_vest10_vis1.mat');
    
    % model_PSTH_trace_with_heter_shortTau = load('diff_PSTH_trace_with_heter_shortTau.mat')
    % model_PSTH_trace_with_heter_shortTau = load('diff_PSTH_trace_with_heter_shortTau_fixed2000ms.mat');
    model_PSTH_trace_with_heter_shortTau = load('diff_PSTH_trace_with_heter_shortTau_fixed100ms.mat'); % 20181010 Fixed real dynamics with 100ms tau
    
%     real_diff_PSTH_trace_MST = load('diff_PSTH_trace_MST_scaling0.78.mat'); % HH20180916  % Scale of sigma
    real_diff_PSTH_trace_MST = load('diff_PSTH_trace_MST_scaling1.mat'); % HH20180916
    
    increase_step_trick = 0; alpha = [];
    common_ts = []; fit_time_range = []; test_time_range = [];  valid_time_range = [];
    fit_deltaPSTHs_per_cell = [];  test_deltaPSTHs_per_cell = []; valid_deltaPSTHs_per_cell = [];
    fit_model_PSTH_interp = [];  test_model_PSTH_interp = []; valid_model_PSTH_interp = [];
    data_to_fit = []; use_data_to_fit = []; swap_visual_vest = []; data_to_fit_time = [];  data_to_fit_PSTH = []; 
    
    function f7p2(debug)
        if debug;  dbstack;  keyboard;  end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        use_data_to_fit = 5; % 
        swap_visual_vest = 0; % Two fellows of Alex suggested us try flip the vest and vis labels HH20180908
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        data_to_fit = {%Times  %Data   %Name 
                       rate_ts{1}, PSTH_all_raw_PrefminusNull{1},'Real LIP data';
                       model_ts, model_PSTH_trace_with_heter_optimal.a,'Optimal model with heterogeneity';
                       model_ts, model_PSTH_trace_with_heter_vest10_vis1.a,'vest10_vis1 with heterogeneity';
                       model_ts, model_PSTH_trace_with_heter_shortTau.a, 'ShortTau with heterogeneity';
                       real_diff_PSTH_trace_MST.PSTH_ts_scaled, permute(real_diff_PSTH_trace_MST.group_PSTH_prefMnull_smoothed,[3 1 2]),'Real MST data';
                       };
                   
        data_to_fit_time = data_to_fit{use_data_to_fit,1};
        data_to_fit_PSTH = data_to_fit{use_data_to_fit,2};
        
        % Times
        find_common_ts_in_rate_ts = find(min(model_ts)<=data_to_fit_time & data_to_fit_time<=max(model_ts));
        common_ts = data_to_fit_time(find_common_ts_in_rate_ts);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % alpha = 10;
        % fit_time_range = 500 < common_ts & common_ts <= 1200;
        % fit_time_range = rand(1,length(common_ts)) < 0.2;

        % Interleaved training and testing
        fit_time_range = false(1,length(common_ts)); 
        test_time_range = fit_time_range;
        valid_time_range = fit_time_range;
        
        fit_time_range(1:3:end) = true;
        test_time_range(2:3:end) = true;
        valid_time_range(3:3:end) = true; % Use other (time points) to test
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fit_ts = common_ts(fit_time_range);
        test_ts = common_ts(test_time_range);
        valid_ts = common_ts(valid_time_range);
        
        
        % Optimal model data
        fit_model_PSTH_interp = interp1(model_ts,model_PSTH_optimal_M1(:,:),fit_ts);
        test_model_PSTH_interp = interp1(model_ts,model_PSTH_optimal_M1(:,:),test_ts);
        valid_model_PSTH_interp = interp1(model_ts,model_PSTH_optimal_M1(:,:),valid_ts);

        % Heter data (basis functions)
        fit_deltaPSTHs_per_cell =  data_to_fit_PSTH(:,find_common_ts_in_rate_ts(fit_time_range),:); % Should be indices in original rate_ts{1}!
        test_deltaPSTHs_per_cell =  data_to_fit_PSTH(:,find_common_ts_in_rate_ts(test_time_range),:);
        valid_deltaPSTHs_per_cell =  data_to_fit_PSTH(:,find_common_ts_in_rate_ts(valid_time_range),:);
        
        
       %% Do fitting
        alphas = [1];
        
        global fitting_dynamics;
        for aa = 1:length(alphas) % Scan alpha
            
            fitting_dynamics = [];
            alpha = alphas(aa);
            fprintf('alpha = %g\n',alpha);
            
            %         w0 = increase_step_trick + ones(1,sum(select_bottom_line))/sum(select_bottom_line); % Start from straight average
            w0 = increase_step_trick + (1+randn(1,size(data_to_fit_PSTH,1)))/size(data_to_fit_PSTH,1); % Start from straight average
            proj_PSTH = reshape(w0 * reshape(fit_deltaPSTHs_per_cell,size(fit_deltaPSTHs_per_cell,1),[]),[],3);
            
            opt = optimset('MaxFunEvals',50000,'MaxIter',300,'PlotFcns',@f7p2_fitting_deltaPSTH_plot_func);
            set(figure(1458),'name','Optimization PlotFcns'); clf;
            %         fitted_w = fminsearch(@(w)f7p2_fitting_deltaPSTH_cost_func(fit_model_PSTH_interp,fit_deltaPSTHs_per_cell,w),w0,opt);
            
            if swap_visual_vest == 0 % Normal labeling
                fitted_w = fmincon(@(w)f7p2_fitting_deltaPSTH_cost_func(fit_model_PSTH_interp,fit_deltaPSTHs_per_cell,w),w0,[],[],[],[],0*w0,[],[],opt);
            else % Total flip, and only fit vest and vis
                ys = fit_model_PSTH_interp;
                ys(:,3) = nan; % Don't fit combined
                
                xs = fit_deltaPSTHs_per_cell(:,:,[2,1,3]);
                xs(:,:,3) = nan;
                
                fitted_w = fmincon(@(w)f7p2_fitting_deltaPSTH_cost_func(ys,xs,w),w0,[],[],[],[],0*w0,[],[],opt);
            end
            
            scan_alpha_costs{aa} = fitting_dynamics{1};
            scan_alpha_weight{aa} = fitted_w;
        end

%         save('Z:\Labtools\HH_Tools\DataHub\scan_alpha.mat','scan_alpha_costs','scan_alpha_weight','alphas');
%         Weighted_sum_PSTH(fitted_w',{'fit'},select_bottom_line)
        
        %% Plotting Scanning alpha result
        %{
        load('Z:\Labtools\HH_Tools\DataHub\scan_alpha.mat');
            
        figure(9041914); clf; subplot(2,2,1);
        for i = 1:5:length(scan_alpha_costs)
            hold on; 
            plot(log(scan_alpha_costs{i}(:,1)),'k'); 
            plot(log(scan_alpha_costs{i}(:,2)),colors(2,:)); 
            plot(log(scan_alpha_costs{i}(:,3)),colors(1,:));
            text(length(scan_alpha_costs{i}(:,3)),log(scan_alpha_costs{i}(end,3)),num2str(alphas(i)));
        end
        xlabel('Steps'); ylabel('log(MSE)');
        
        subplot(2,2,2);
        weights = cell2mat(scan_alpha_weight');
        [coeff,score,latent,~,explained] = pca(weights);
        for i=1:length(scan_alpha_costs)
            hold on; scatter(score(i,1),score(i,2),i*3,'ob');
        end
        xlabel('PC1'); ylabel('PC2');
        
        subplot(2,2,3);
        min_cost = cellfun(@(x)min(x(:,3)),scan_alpha_costs)';
        for i=1:length(scan_alpha_costs)
            hold on; scatter(alphas(i),log(min_cost(i)),i*3,'ob');
        end
        xlabel('alpha'); ylabel('log(validating MSE)'); 
        SetFigure(15);
        %}
    end
    function stop = f7p2_fitting_deltaPSTH_plot_func(weight,optimValue,~)
        global fitting_dynamics;
        weight = weight-increase_step_trick;

        persistent costs;
        if optimValue.iteration == 0
            costs = [];
        end
        set(findobj(gcf,'type','axes'),'visible','off');
        h = tight_subplot(2,2,[0.1 0.1],0.1,0.1);
        
        % Update cost
        proj_test_PSTH = reshape((weight) * reshape(test_deltaPSTHs_per_cell,size(test_deltaPSTHs_per_cell,1),[]),[],3);
        test_cost = mean((proj_test_PSTH(:) - test_model_PSTH_interp(:)).^2); % Mean Squared error
        proj_valid_PSTH = reshape((weight) * reshape(valid_deltaPSTHs_per_cell,size(valid_deltaPSTHs_per_cell,1),[]),[],3);
        valid_cost = mean((proj_valid_PSTH(:) - valid_model_PSTH_interp(:)).^2); % Mean Squared error
        costs = [costs; optimValue.fval-alpha*norm(weight) test_cost valid_cost];
        
        if mod(optimValue.iteration,1) == 0
            % Fitted delta PSTH
            %             plot(h(1),real_ts,model_PSTH_interp,'linew',2); hold(h(1),'on');
            %             plot(h(1),real_ts,reshape((weight) * reshape(real_deltaPSTHs_per_cell,size(real_deltaPSTHs_per_cell,1),[]),[],3),'--','linew',2);
            plot(h(1),data_to_fit_time,interp1(model_ts,model_PSTH_optimal_M1(:,:),data_to_fit_time),'linew',2); hold(h(1),'on');
            plot(h(1),data_to_fit_time,reshape((weight) * reshape(data_to_fit_PSTH,size(data_to_fit_PSTH,1),[]),[],3),'--','linew',2);
            xlim(h(1),[-400 2400]);
            ylim(h(1),[-5 20]);
            ylims = ylim((h(1))); indicator_pos = min(ylims)*ones(1,length(common_ts));
            plot(h(1),[min(model_ts) min(model_ts)],[min(ylims) max(ylims)],'k-');
            plot(h(1),[max(model_ts) max(model_ts)],[min(ylims) max(ylims)],'k-');
            plot(h(1),common_ts(fit_time_range),indicator_pos(fit_time_range),'ok','markerfacecol','k','markersize',5);
            plot(h(1),common_ts(test_time_range),indicator_pos(test_time_range),'or','markerfacecol',colors(2,:),'markersize',5);
            plot(h(1),common_ts(valid_time_range),indicator_pos(valid_time_range),'ob','markerfacecol',colors(1,:),'markersize',5);
            % plot(h(1),[min(real_ts) max(real_ts)],min(ylims)*ones(1,2),'k-','linew',5);
            title(h(1),sprintf('Solid: model; Dashed: %s; Swap = %g',data_to_fit{use_data_to_fit,3}, swap_visual_vest));
            ylabel(h(1),'delta firing');
            
            % Projected PSTH
            if use_data_to_fit == 1 % Real LIP data
                for kk = 1:size(PSTH_all_raw{1},3)
                    PSTH_projected(1,:,kk) = weight*(PSTH_all_raw{1}(select_bottom_line,:,kk));
                end
                
                % Note here the errorbars should be STD instead of SEM. (Not independent sampling, but bootstrapping)
                SeriesComparison(PSTH_projected,data_to_fit_time,...
                    'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),colors(3,:),colors(3,:)},'LineStyles',{'-','--','-','--','-','--'},...
                    'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Weighted sum of firing','axes',h(2));
                hold on;    legend off;
                xlim(h(2),[-400 2200]);
            end
            
            % Fitting error
            plot(h(3), 1:size(costs,1), log(costs(:,1)),'ok-');  % Training error
            hold(h(3),'on');
            plot(h(3), 1:size(costs,1), log(costs(:,2)),'or-');  % Testing error
            plot(h(3), 1:size(costs,1), log(costs(:,3)),'ob-');  % Testing error

            legend(h(3),{'Training','Testing','Validating'});
            title(h(3),sprintf('\\alpha = %g, step = %g',alpha,optimValue.iteration));
            ylabel(h(3),'log(mean of squared error)');
            xlabel(h(3),'# Iteration');
            
            % Weight vs choice signal
            if use_data_to_fit == 1 % Real LIP data
                xx =  abs(group_ChoicePreference(3,select_bottom_line,3));
                plot(h(4),xx,weight,'o'); hold(h(4),'on');
                
                [r,p] = corr(xx',weight');
                [linPara,S] = polyfit(xx,weight,1);
                xxx = linspace(min(xx),max(xx),2);
                Y = polyval(linPara,xxx);
                plot(h(4),xxx,Y,'k');
                title(h(4),sprintf('p = %3.3g, r = %3.3g',p,r));
                
                xlabel(h(4),'Choice strength of each cell');
                ylabel(h(4),'Fitting weight');
            else % Modeled data
                plot(h(4),weight,'o'); hold(h(4),'on');
                xlabel(h(4),'# cell');
                ylabel(h(4),'Fitting weight');
                
            end
%             keyboard
        end
        stop = false;
        
        % Save data
        fitting_dynamics = {costs};
        
    end
    function cost = f7p2_fitting_deltaPSTH_cost_func(model_y,real_ys,weight)
        weight = weight-increase_step_trick;
        proj_PSTH = reshape((weight) * reshape(real_ys,size(real_ys,1),[]),[],3);
        cost = nanmean((proj_PSTH(:) - model_y(:)).^2)+alpha*norm(weight); % Mean Squared error
    end


    %% Regression 3. Fitting model traces with real data (New version): Use traces of each |heading| separately
    %% HH20180623: D2206 Shanghai -> Chengdu
    %%%%%%%%%%%%%%%%%%% for f7p3: load model traces %%%%%%%%%%%%%%%
    
    % ==== Load model data (Fit target) from the trace without heterogeneity ====
    model_mean_PSTH_trace_without_heter_EachHeading = load('diff_PSTH_eachHeadingTest.mat');
    model_PSTH_EachHeading = squeeze(mean(model_mean_PSTH_trace_without_heter_EachHeading.diff_PSTH_correct_mean_headings,1));
    model_PSTH_EachHeading = permute(model_PSTH_EachHeading,[1 3 2]); % [Time, |Heading|, modality]
    model_ts_EachHeading = linspace(0, 1500, size(model_PSTH_EachHeading,1)); % Temporally
    
    % ==== Single cells from the model ====
    model_PSTH_trace_with_heter_optimal_EachHeading = [];% load('diff_PSTH_trace_with_heter.mat');
    model_PSTH_trace_with_heter_vest10_vis1_EachHeading = []; %load('diff_PSTH_trace_with_heter_vest10_vis1.mat');
    model_PSTH_trace_with_heter_shortTau_EachHeading = []; %load('diff_PSTH_trace_with_heter_shortTau.mat');
    
    data_to_fit_PSTH_EachHeading=[];
    
    function f7p3(debug)  % Fit each Heading
        if debug;  dbstack;  keyboard;  end
        
        use_data_to_fit = 1; % 
                
        data_to_fit = {%Times  %Data   %Name 
                       % rate_ts{1}, PSTH_all_raw_PrefminusNull{1},'Real LIP data';  
                       rate_ts{1},  PSTH_correct_angles_raw{1}(:,:,1:2:end,:) - PSTH_correct_angles_raw{1}(:,:,2:2:end,:),'Real LIP data'; 
                       model_ts_EachHeading, model_PSTH_trace_with_heter_optimal.a,'Optimal model with heterogeneity';
                       model_ts_EachHeading, model_PSTH_trace_with_heter_vest10_vis1.a,'vest10_vis1 with heterogeneity';
                       model_ts_EachHeading, model_PSTH_trace_with_heter_shortTau.a, 'ShortTau with heterogeneity';
                       };
        data_to_fit_time = data_to_fit{use_data_to_fit,1};
        data_to_fit_PSTH_EachHeading = data_to_fit{use_data_to_fit,2};
        
        % Times
        find_common_ts_in_rate_ts = find(min(model_ts_EachHeading)<=data_to_fit_time & data_to_fit_time<=max(model_ts_EachHeading));
        common_ts = data_to_fit_time(find_common_ts_in_rate_ts);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % alpha = 10;
        % fit_time_range = 500 < common_ts & common_ts <= 1200;
        % fit_time_range = rand(1,length(common_ts)) < 0.2;

        % Interleaved training and testing
        fit_time_range = false(1,length(common_ts)); 
        test_time_range = fit_time_range;
        valid_time_range = fit_time_range;
        
        fit_time_range(1:3:end) = true;
        test_time_range(2:3:end) = true;
        valid_time_range(3:3:end) = true; % Use other (time points) to test
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fit_ts = common_ts(fit_time_range);
        test_ts = common_ts(test_time_range);
        valid_ts = common_ts(valid_time_range);
        
        % Real data, now: [NumOfCells, Time, |Headings|, Modality]

        fit_deltaPSTHs_per_cell =  data_to_fit_PSTH_EachHeading(:,find_common_ts_in_rate_ts(fit_time_range),:,:); % Should be indices in original rate_ts{1}!
        test_deltaPSTHs_per_cell =  data_to_fit_PSTH_EachHeading(:,find_common_ts_in_rate_ts(test_time_range),:);
        valid_deltaPSTHs_per_cell =  data_to_fit_PSTH_EachHeading(:,find_common_ts_in_rate_ts(valid_time_range),:);
        
        % Model data,  now: [Time, |Headings|, Modality]
        fit_model_PSTH_interp = interp1(model_ts_EachHeading,model_PSTH_EachHeading(:,:,:),fit_ts);
        test_model_PSTH_interp = interp1(model_ts_EachHeading,model_PSTH_EachHeading(:,:,:),test_ts);
        valid_model_PSTH_interp = interp1(model_ts_EachHeading,model_PSTH_EachHeading(:,:,:),valid_ts);
        
       %% Do fitting  
        alphas = [1];
        
        global fitting_dynamics;
        for aa = 1:length(alphas)% Scan alpha
            
            fitting_dynamics = [];
            alpha = alphas(aa);
            fprintf('alpha = %g\n',alpha);
            
            %         w0 = increase_step_trick + ones(1,sum(select_bottom_line))/sum(select_bottom_line); % Start from straight average
            w0 = increase_step_trick + (1+randn(1,size(data_to_fit_PSTH_EachHeading,1)))/size(data_to_fit_PSTH_EachHeading,1); % Start from straight average
            
            opt = optimset('MaxFunEvals',50000,'MaxIter',300,'PlotFcns',@f7p3_fitting_deltaPSTH_plot_func);
            set(figure(1458),'name','Optimization PlotFcns'); clf;
            %         fitted_w = fminsearch(@(w)f7p2_fitting_deltaPSTH_cost_func(fit_model_PSTH_interp,fit_deltaPSTHs_per_cell,w),w0,opt);
            
            fitted_w = fmincon(@(w)f7p3_fitting_deltaPSTH_cost_func(fit_model_PSTH_interp,fit_deltaPSTHs_per_cell,w),w0,...
                                                                                                      [],[],[],[],0*w0,[],[],opt);
        end

    end

    function stop = f7p3_fitting_deltaPSTH_plot_func(weight,optimValue,~)  % Plotting function
        
        global fitting_dynamics;
        weight = weight-increase_step_trick;

        persistent costs;
        
        if optimValue.iteration == 0
            costs = [];
        end
        
        
        % Update cost
        test_cost = f7p3_fitting_deltaPSTH_cost_func(test_model_PSTH_interp , test_deltaPSTHs_per_cell, weight);
        valid_cost = f7p3_fitting_deltaPSTH_cost_func(valid_model_PSTH_interp , valid_deltaPSTHs_per_cell, weight);
        
%         proj_test_PSTH = reshape((weight) * reshape(test_deltaPSTHs_per_cell,size(test_deltaPSTHs_per_cell,1),[]),[],3);
%         test_cost = mean((proj_test_PSTH(:) - test_model_PSTH_interp(:)).^2); % Mean Squared error
%         
%         proj_valid_PSTH = reshape((weight) * reshape(valid_deltaPSTHs_per_cell,size(valid_deltaPSTHs_per_cell,1),[]),[],3);
%         valid_cost = mean((proj_valid_PSTH(:) - valid_model_PSTH_interp(:)).^2); % Mean Squared error
        
        costs = [costs; optimValue.fval-alpha*norm(weight) test_cost valid_cost];
        
        if mod(optimValue.iteration,1) == 0
        
            set(findobj(gcf,'type','axes'),'visible','off');
            set(gcf,'uni','norm','pos',[0.017       0.208       0.947       0.675]);
            h = tight_subplot(2,4,[0.1 0.05],0.1,0.1);
            
            % ------ Fitted delta PSTH for each heading
            headings = [0 1 2 4 8];
            for hh = 1:5
                plot(h(hh),data_to_fit_time,interp1(model_ts_EachHeading,squeeze(model_PSTH_EachHeading(:,hh,:)),...
                    data_to_fit_time),'linew',2); hold(h(hh),'on');
                
                
                allY_this = reshape(data_to_fit_PSTH_EachHeading(:,:,hh,:),size(data_to_fit_PSTH_EachHeading,1),[]);
                proj_PSTH_this = reshape(nansum( weight' .* allY_this),[],3);
                
                plot(h(hh),data_to_fit_time, proj_PSTH_this,'--','linew',2);
                
                xlim(h(hh),[-400 2400]);
                ylim(h(hh),[-5 40]);
                
                plot(h(hh),[min(model_ts_EachHeading) min(model_ts_EachHeading)],[min(ylim) max(ylim)],'k-');
                plot(h(hh),[max(model_ts_EachHeading) max(model_ts_EachHeading)],[min(ylim) max(ylim)],'k-');
                
                title(h(hh),sprintf('|Heading| = %g',headings(hh)));
            end
            
            %             plot(h(1),common_ts(fit_time_range),indicator_pos(fit_time_range),'ok','markerfacecol','k','markersize',5);
            %             plot(h(1),common_ts(test_time_range),indicator_pos(test_time_range),'or','markerfacecol',colors(2,:),'markersize',5);
            %             plot(h(1),common_ts(valid_time_range),indicator_pos(valid_time_range),'ob','markerfacecol',colors(1,:),'markersize',5);
            %             % plot(h(1),[min(real_ts) max(real_ts)],min(ylims)*ones(1,2),'k-','linew',5);
            title(h(1),sprintf('Solid: model; Dashed: %s',data_to_fit{use_data_to_fit,3}));
            ylabel(h(1),'delta firing');
            
            % ------- Projected PSTH
            if use_data_to_fit == 1 % Real LIP data
                for kk = 1:size(PSTH_all_raw{1},3)
                    PSTH_projected(1,:,kk) = weight*(PSTH_all_raw{1}(select_bottom_line,:,kk));
                end
                
                % Note here the errorbars should be STD instead of SEM. (Not independent sampling, but bootstrapping)
                SeriesComparison(PSTH_projected,data_to_fit_time,...
                    'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),colors(3,:),colors(3,:)},'LineStyles',{'-','--','-','--','-','--'},...
                    'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Weighted sum of firing','axes',h(8));
                hold on;    legend off;
                xlim(h(2),[-400 2200]);
            end
            
            % ------- Fitting error
            plot(h(6), 1:size(costs,1), log(costs(:,1)),'ok-');  % Training error
            hold(h(6),'on');
            plot(h(6), 1:size(costs,1), log(costs(:,2)),'or-');  % Testing error
            plot(h(6), 1:size(costs,1), log(costs(:,3)),'ob-');  % Testing error

            legend(h(6),{'Training','Testing','Validating'});
            title(h(6),sprintf('\\alpha = %g, step = %g',alpha,optimValue.iteration));
            ylabel(h(6),'log(mean of squared error)');
            xlabel(h(6),'# Iteration');
            
            % ------- Weight vs choice signal
            if use_data_to_fit == 1 % Real LIP data
                xx =  abs(group_ChoicePreference(3,select_bottom_line,3));
                plot(h(7),xx,weight,'o'); hold(h(7),'on');
                
                [r,p] = corr(xx',weight');
                [linPara,S] = polyfit(xx,weight,1);
                xxx = linspace(min(xx),max(xx),2);
                Y = polyval(linPara,xxx);
                plot(h(7),xxx,Y,'k');
                title(h(7),sprintf('p = %3.3g, r = %3.3g',p,r));
                
                xlabel(h(7),'Choice strength of each cell');
                ylabel(h(7),'Fitting weight');
            else % Modeled data
                plot(h(7),weight,'o'); hold(h(7),'on');
                xlabel(h(7),'# cell');
                ylabel(h(7),'Fitting weight');
                
            end
%             keyboard
        end
        stop = false;
        
        % Save data
        fitting_dynamics = {costs};
        
    end
    function cost = f7p3_fitting_deltaPSTH_cost_func(model_y,real_ys,weight)
        weight = weight-increase_step_trick;
        
        %      proj_PSTH = reshape((weight) * reshape(real_ys,size(real_ys,1),[]),[],3);
        
        % Some cells didn't have 0 headings, so I use nansum but not matrix product.
        allY = reshape(real_ys,size(real_ys,1),[]);
        proj_PSTH = reshape( nansum(repmat(weight',1,size(allY,2)) .* allY), ...
                            size(model_y,1), size(model_y,2), size(model_y,3));
        
        cost = mean((proj_PSTH(:) - model_y(:)).^2) + alpha*norm(weight); % Mean Squared error
        
    end


    fisherSimpleGu = [];
    findForFisher = [];

    %% ====== Calculate Fisher information of heading ======== HH20180619
    % 1. Simplest method like Gu 2010: Sum over cells (slope / mean)
    function f6p5p1(debug)  
        if debug;  dbstack;  keyboard;  end
        %%
        
        j = 1;
                
        findForFisher = find(select_bottom_line);

        fisherSize = size(group_result(1).mat_raw_PSTH.fisherSimpleGu);
        fisherSimpleGu = nan(length(findForFisher),fisherSize(1),fisherSize(2),fisherSize(3),fisherSize(4));  % Last one: varPoisson/varReal
        
        for ii = 1:length(findForFisher)
            
            thisCell = findForFisher(ii);
            fisherSimpleGu (ii,:,:,:,:) = group_result(thisCell).mat_raw_PSTH.fisherSimpleGu;

            % Moved to batch file
            %{ 
            for k = 1:3
               
                tmpForPar = group_result(thisCell).mat_raw_PSTH.CP{j,k}.raw_CP_result;
              
                % parfor tt = 1:length(CP_ts{j})   % Parfor is even slower...
                for tt = 1:length(CP_ts{j})
                   
                    % Get data
                    thisTuning = tmpForPar{tt}.Neu_tuning;
                    if thisPref == LEFT % To keep PREF = RIGHT (Actually no effect because we have slope^2 afterwards in Fisher information)
                        thisTuning(:,2:3) = flipud(thisTuning(:,2:3));
                    end
                    linearFit = polyfit(thisTuning(:,1),thisTuning(:,2),1);
                    slopeInRad = linearFit(1)*(180/pi);
                    varPoisson = mean(thisTuning(:,2)); % Assuming Poisson with fano = 1
                    varReal = mean((thisTuning(:,3) * sqrt(thisN)).^2); % Using real variance (usually much large)
                    
                    % Compute simple Fisher
                    thisFisherVarPoisson = slopeInRad^2/varPoisson;
                    thisFisherVarReal = slopeInRad^2/varReal;
                    
                    % Save
                    tmpOut(tt,1,:) = [thisFisherVarPoisson, thisFisherVarReal];
                end
                
                fisherSimpleGu (ii,:,k,:) = tmpOut;
            end
            %}
            
        end
        
        % -------- Plotting FI(t) ---------
        figure(4232); clf
        set(gcf,'uni','norm','pos',[0.442       0.166       0.557       0.739]);
        plotRange = 1:find(CP_ts{1}>=1600,1);
        bootN = 2000;
        
        rangeTitles = {'+/-8 degree','+/-2 degrees'};
        titles = {'Poisson assump.','Real variance'};
        
        for ranges = 1:2
            for varMethod = 1:2  %  1: Gu's Poisson assumption   2: Real variance (10x larger than Poisson because of varCE)
                h = subplot(2,2,varMethod + (ranges-1)*2);
                
                % Get sum of Fisher and std by bootstrap (Gu 2010)
                sumBoots = bootstrp(bootN,@(x)sum(x,1),squeeze(fisherSimpleGu(:,plotRange,:,ranges,varMethod)));
                sumBoots = reshape(sumBoots,bootN,[],3);
                
                SeriesComparison(sumBoots, CP_ts{1}(plotRange),...
                    'SEM',0,'Errorbar',2,'Axes',h,'Transparent',transparent);
                
                sumFisherMean = squeeze(mean(sumBoots,1));
                % sumFisherVestPlusVIs = sumFisherMean(:,1) + sumFisherMean(:,2);
                sumFisherVestPlusVIs = sumFisherMean(:,1) - nanmean(sumFisherMean(CP_ts{1}<0,1)) ...
                                     + sumFisherMean(:,2) - nanmean(sumFisherMean(CP_ts{1}<0,2)) ...
                                     + nanmean(sumFisherMean(CP_ts{1}<0,3)); % Remove baseline
                                 
                plot(CP_ts{1}(plotRange),sumFisherVestPlusVIs,'m','linew',2);
                
                % Gaussian vel
                axis tight; plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
                xlim([-100 1600]); ylim([0 max(ylim)*1.1])
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
                legend off;
                title(titles{varMethod})
                ylabel(rangeTitles{ranges});
                text(min(xlim),min(ylim),sprintf('n = %g',length(findForFisher)))
            end
        end
        
        SetFigure(15)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FIScatterTimeRange = 500 < CP_ts{1} & CP_ts{1} <= 1500;
    FIScatterHeadingRange = 1; % +/-8 degrees
    FIScatterVariance = 2; % Real variance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1.1 Plot scatter FI of different modalities
    function f6p5p1p1(debug)  
        if debug;  dbstack;  keyboard;  end
        
        findForFisher = find(select_tcells);

        if isempty(fisherSimpleGu), f6p5p1(0); end
        
        % -------- Plotting cell-by-cell FI (Yong Gu wants this) ---------
 
        % Which info?
        FI_all_temp = squeeze(mean(fisherSimpleGu(:,FIScatterTimeRange,:,FIScatterHeadingRange,FIScatterVariance),2))';
        FI_temp_comb_minus_vest = FI_all_temp(3,:)-FI_all_temp(1,:);
        FI_temp_comb_minus_vis = FI_all_temp(3,:)-FI_all_temp(2,:);
        FI_temp_vest_plus_vis = FI_all_temp(1,:) + FI_all_temp(2,:);
        
        %%  1. ====== vest and visual  ========
        tt = 3;
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 5; monkey1 = monkey1(findForFisher)';
        monkey2 = monkeys == 10; monkey2 = monkey2(findForFisher)';

        set(figure(figN),'name','FI (visual) vs. FI (vest)','pos',[17 514 1151 449]);
        
        cpref_sig_1 = Choice_pref_p_value_all(1,findForFisher,tt) < 0.05;  % Use Choice pref p value as indicators
        cpref_sig_2 = Choice_pref_p_value_all(2,findForFisher,tt) < 0.05;
        
        
        h = LinearCorrelation({
                (FI_all_temp(2, monkey1 & ~cpref_sig_1 & ~cpref_sig_2)) ;
                (FI_all_temp(2, monkey2 & ~cpref_sig_1 & ~cpref_sig_2)) ;
                (FI_all_temp(2, monkey1 & xor(cpref_sig_1 , cpref_sig_2)));
                (FI_all_temp(2, monkey2 & xor(cpref_sig_1 , cpref_sig_2)));
                (FI_all_temp(2, monkey1 & cpref_sig_1 & cpref_sig_2));...
                (FI_all_temp(2, monkey2 & cpref_sig_1 & cpref_sig_2))},...
                {
                (FI_all_temp(1,monkey1 & ~cpref_sig_1 & ~cpref_sig_2)) ;
                (FI_all_temp(1,monkey2 & ~cpref_sig_1 & ~cpref_sig_2)) ;
                (FI_all_temp(1,monkey1 & xor(cpref_sig_1 , cpref_sig_2))) ;
                (FI_all_temp(1,monkey2 & xor(cpref_sig_1 , cpref_sig_2))) ;
                (FI_all_temp(1,monkey1 & cpref_sig_1 & cpref_sig_2)) ;...
                (FI_all_temp(1,monkey2 & cpref_sig_1 & cpref_sig_2)) },...
                'CombinedIndex',[63],'PlotCombinedOnly', 1,...
                'Ylabel','Vestibular FI','Xlabel','Visual FI',...
                'FaceColors',{'none','none',[0.8 0.8 0.8],[0.8 0.8 0.8],'k','k'},'Markers',{'o','^'},...
                'LineStyles',{'k:','k:','k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,'Method','Pearson','FittingMethod',2, ...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1); figN = figN + 1;
        
        delete([h.group(1:6).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       %  SetFigure(20);
        
%         % Annotate tcells
%         plot(FI_all_temp(2,select_tcells(select_tcells)),FI_all_temp(1,select_tcells(select_tcells)),...
%             '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
                
        % Show individual cell selected from the figure. HH20150424
        h_line = plot((FI_all_temp(2,:)),(FI_all_temp(1,:)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_tcells});
        
        
       %% 2. ====== Single vs Comb =======
        two_face_colors = fliplr({'none','none',[0.8 0.8 1],[0.8 0.8 1],colors(1,:),colors(1,:);
                                  'none','none',[1 0.8 0.8],[1 0.8 0.8],colors(2,:),colors(2,:)});
        
        for k = 1:2  % Plot it separately
            
            set(figure(figN),'name','FI (single) vs. FI (comb)','pos',[17 514 1151 449]);
            
            cpref_sig_k = Choice_pref_p_value_all(k,findForFisher,tt) < 0.05;
%             cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
            cpref_sig_3 = Choice_pref_p_value_all(3,findForFisher,tt) < 0.05;
           
           
            h = LinearCorrelation({
                (FI_all_temp(k, monkey1 & cpref_sig_k & cpref_sig_3));
                (FI_all_temp(k, monkey2 & cpref_sig_k & cpref_sig_3));
                (FI_all_temp(k, monkey1 & xor(cpref_sig_k , cpref_sig_3)));
                (FI_all_temp(k, monkey2 & xor(cpref_sig_k , cpref_sig_3)));
                (FI_all_temp(k, monkey1 & ~cpref_sig_k & ~cpref_sig_3)) ;
                (FI_all_temp(k, monkey2 & ~cpref_sig_k & ~cpref_sig_3)) ;
                },...
                {
                (FI_all_temp(3,monkey1 & cpref_sig_k & cpref_sig_3)) ;
                (FI_all_temp(3,monkey2 & cpref_sig_k & cpref_sig_3));
                (FI_all_temp(3,monkey1 & xor(cpref_sig_k, cpref_sig_3))) ;
                (FI_all_temp(3,monkey2 & xor(cpref_sig_k, cpref_sig_3)));
                (FI_all_temp(3,monkey1 & ~cpref_sig_k & ~cpref_sig_3)) ;
                (FI_all_temp(3,monkey2 & ~cpref_sig_k & ~cpref_sig_3)) ;
                },...
                'CombinedIndex',[63],'PlotCombinedOnly', 1,...
                'Ylabel','Combined FI','Xlabel','Single FI',...
                'FaceColors',two_face_colors(k,:),'Markers',{'o','^'},...
                'LineStyles',{'k:','k:','k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,...
                'Method','Pearson','FittingMethod',2); figN = figN + 1;
            
            % delete([h.group(1:6).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       %   SetFigure(20);
            axis square;
            
            % Annotate tcells
%             plot((FI_all_temp(k,select_tcells(select_tcells))),(FI_all_temp(3,select_tcells(select_tcells))),...
%                 '+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot((FI_all_temp(k,:)),(FI_all_temp(3,:)),'visible','off'); hold on;
            set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_tcells],1});
                                   
            axis tight;
        end
        %             for i = 1:3
        %                 set(h.group(i).dots,'color',colors(k,:));
        %             end
        
        %% ===  3. (Comb - visual) VS (Comb - vest)  ===
        set(figure(figN),'name','FI(Comb - visual) VS FI(Comb - vest)','pos',[17 514 1151 449]);
                
        cellTypes = [group_result.Waveform_broad];
        cellTypes = cellTypes(findForFisher);
        cellTypes(:) = 0;

        
        h = LinearCorrelation({
            (FI_temp_comb_minus_vest(1,monkey1 & ~cellTypes)) ;
            (FI_temp_comb_minus_vest(1,monkey2 & ~cellTypes)) ;
            (FI_temp_comb_minus_vest(1,monkey1 & cellTypes)) ;
            (FI_temp_comb_minus_vest(1,monkey2 & cellTypes)) ;
            },...
            {
            (FI_temp_comb_minus_vis(1,monkey1 & ~cellTypes)) ;
            (FI_temp_comb_minus_vis(1,monkey2 & ~cellTypes)) ;
            (FI_temp_comb_minus_vis(1,monkey1 & cellTypes)) ;
            (FI_temp_comb_minus_vis(1,monkey2 & cellTypes)) ;
            },...
            'CombinedIndex',15,'PlotCombinedOnly', 1, ...
            'Xlabel','FI (Combined - Vest)','Ylabel','FI (Combined - Visual)',...
            'FaceColors',{'none','none','k','k'},'Markers',{'o','^'},...
            'LineStyles',{'k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
            'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
            'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        axis square;
        
        % Annotate tcells
%         plot((FI_temp_comb_minus_vest(1,select_tcells(select_bottom_line))),...
%             (FI_temp_comb_minus_vis(1,select_tcells(select_bottom_line))),...
%             '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        %             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
        %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(FI_temp_comb_minus_vest(1,:),FI_temp_comb_minus_vis(1,:),'visible','off'); hold on;
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_tcells],1});
        
        axis tight;   
        
        
       %% ===  4. Comb VS (Vest + Vis)  ===
        set(figure(figN),'name','FI(Comb) VS FI(Vest + Vis)','pos',[17 514 1151 449]);
                
        cellTypes = [group_result.Waveform_broad];
        cellTypes(:) = 0;
        
        cellTypes = cellTypes(findForFisher);
        
        h = LinearCorrelation({
            (FI_temp_vest_plus_vis(1,monkey1 & ~cellTypes)) ;
            (FI_temp_vest_plus_vis(1,monkey2 & ~cellTypes)) ;
            (FI_temp_vest_plus_vis(1,monkey1 & cellTypes)) ;
            (FI_temp_vest_plus_vis(1,monkey2 & cellTypes)) ;
            },...
            {
            (FI_all_temp(3,monkey1 & ~cellTypes)) ;
            (FI_all_temp(3,monkey2 & ~cellTypes)) ;
            (FI_all_temp(3,monkey1 & cellTypes)) ;
            (FI_all_temp(3,monkey2 & cellTypes)) ;
            },...
            'CombinedIndex',15,'PlotCombinedOnly', 1, ...
            'Xlabel','FI (Vest + Vis)','Ylabel','FI (Combined)',...
            'FaceColors',{'none','none','k','k'},'Markers',{'o','^'},...
            'LineStyles',{'k:','k:','k:','k:','k-'},'MarkerSize',marker_size,'SameScale',1,...
            'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
            'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        axis square;
        
        % Annotate tcells
%         plot((FI_temp_comb_minus_vest(1,select_tcells(select_bottom_line))),...
%             (FI_temp_comb_minus_vis(1,select_tcells(select_bottom_line))),...
%             '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        %             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
        %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(FI_temp_vest_plus_vis(1,:),FI_all_temp(3,:),'visible','off'); hold on;
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_tcells],1});
        
        clear Choice_pref_all_temp;
        axis tight;     
        
        
    end
    
    % 1.2. FI VS CPref ====
    function f6p5p1p2(debug)
        if debug;  dbstack;  keyboard;  end
        if isempty(fisherSimpleGu), f6p5p1(0); end
        
        monkeys = xls_num{1}(:,header.Monkey);  
        monkey1 = monkeys == 5; monkey1 = monkey1(findForFisher)';
        monkey2 = monkeys == 10; monkey2 = monkey2(findForFisher)';

        % Which info?
        FI_all_temp = squeeze(mean(fisherSimpleGu(:,FIScatterTimeRange,:,FIScatterHeadingRange,FIScatterVariance),2))';
        tt = 3;
        
        for k = 1:3  % Plot it separately
            
            set(figure(figN),'name','FI vs. Cpref','pos',[17 514 1151 449]);
            
            cpref_sig_k = Choice_pref_p_value_all(k,findForFisher,tt) < 0.05;
            
            Choice_pref_all_temp = abs(Choice_pref_all(:,findForFisher,tt));
           
            h = LinearCorrelation({
                (Choice_pref_all_temp(k, monkey1 & cpref_sig_k ));
                (Choice_pref_all_temp(k, monkey2 & cpref_sig_k ));
                (Choice_pref_all_temp(k, monkey1 & ~cpref_sig_k)) ;
                (Choice_pref_all_temp(k, monkey2 & ~cpref_sig_k )) ;
                },...
                {
                (FI_all_temp(k,monkey1 & cpref_sig_k)) ;
                (FI_all_temp(k,monkey2 & cpref_sig_k ));
                (FI_all_temp(k,monkey1 & ~cpref_sig_k)) ;
                (FI_all_temp(k,monkey2 & ~cpref_sig_k)) ;
                },...
                'CombinedIndex',15,'PlotCombinedOnly', 1, ...
                'Ylabel','FI','Xlabel','ChoicePref',...
                'FaceColors',{colors(k,:),colors(k,:),'none','none'},'Markers',{'o','^'},...
                'LineStyles',{'k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked',...
                'Method','Spearman','FittingMethod',2); figN = figN + 1;
            
            % delete([h.group(1:6).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       %   SetFigure(20);
            
            % Annotate tcells
%             plot((FI_all_temp(k,select_tcells(select_tcells))),(FI_all_temp(3,select_tcells(select_tcells))),...
%                 '+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(Choice_pref_all_temp(k,:),FI_all_temp(k,:),'visible','off'); hold on;
            set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_tcells],1});
                                   
            axis tight;
        end
        %             for i = 1:3
        %                 set(h.group(i).dots,'color',colors(k,:));
        %             end
        
    end
        
    FIDora_partial_corr_all = []; FIDora_partial_corr_all_flipped = []; FIDora_partial_corr_all_Rsquare = []; 
    FIDora_conditional_slope2overVar = []; FIDora_conditional_choiceslope2overVar = [];
    
    function f6p5p2(debug)  % Dora's partial corr / multivariate linear regression slope
        if debug;  dbstack;  keyboard;  end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
        select_for_partial = select_tcells;
        
        FIDora_binSize = 200; % in ms
        FIDora_stepSize = 10;
        FIDora_tCenters = FIDora_binSize/2 : FIDora_stepSize : 1500 - FIDora_binSize/2 ;
        conditional_variance_min_reps = 10;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        find_for_partial = find(select_for_partial);
        
        % Calculate partial correlation coefficients for each cell
        
        if isempty(FIDora_partial_corr_all)
            
            FIDora_partial_corr_all = nan(sum(select_for_partial),length(FIDora_tCenters),2,3);  % [cell number, time epoch, (heading coeff, choice coeff), stim_type]
            FIDora_partial_corr_all_flipped = nan(sum(select_for_partial),length(FIDora_tCenters),2,3);  % For plotting partial corr over time
            FIDora_conditional_slope2overVar = nan(sum(select_for_partial),length(FIDora_tCenters),1,3);  
            
            progressbar('cell num');
            
            for i = 1:sum(select_for_partial)  % For cells
                
                this_raw_spike_in_bin = group_result(find_for_partial(i)).mat_raw_PSTH.spike_aligned{1,j};
                this_time = group_result(find_for_partial(i)).mat_raw_PSTH.spike_aligned{2,j};
                this_stim_type = group_result(find_for_partial(i)).mat_raw_PSTH.stim_type_per_trial;
                this_heading = group_result(find_for_partial(i)).mat_raw_PSTH.heading_per_trial;
                this_choice = group_result(find_for_partial(i)).mat_raw_PSTH.choice_per_trial;
                
                for tt = 1 : length(FIDora_tCenters)
                    
                    count_win = FIDora_tCenters(tt) - FIDora_binSize/2 <= this_time & this_time <= FIDora_tCenters(tt) + FIDora_binSize/2;
                    
                    for k = 1:3
                        if isempty(find(this_stim_type==k, 1)); continue; end
                        
                        % --- 1. Partial correlation ---
                        X=[];
                        X(:,1) = sum(this_raw_spike_in_bin(this_stim_type==k,count_win),2)...
                            / FIDora_binSize * 1e3; % Average firing rate in Hz
                        X(:,2) = this_heading(this_stim_type==k);
                        X(:,3) = this_choice(this_stim_type==k);
                        
                        r = partialcorr(X);
                        FIDora_partial_corr_all(i,tt,:,k) = r(1,2:3);
                        
                        % --- 2. Partial FI ---
                        % -     2.1 Multivariate linear regression of beta ---
                        normalizedX3 = X(:,3)./rms(X(:,3)).*rms(X(:,2)); % Normalize choice to headings' RMS (to keep sensory heading comparable with the traditional FI)
                        coeff = glmfit([X(:,2) normalizedX3],X(:,1));
                        conditional_sensory_slope = coeff(2) * (180/pi); % Turn to rad
                        conditional_choice_slope = coeff(3) * (180/pi); % Turn to rad, although it's choice 
                        
                        % -     2.2 Compute conditional variance
                        unique_heading = unique(this_heading);
                        conditional_variance_matrix = nan(2,length(unique_heading));
                        
                        for hh = 1:length(unique(this_heading))
                            for cc = LEFT:RIGHT
                                this_rates = X( X(:,2) == unique_heading(hh) & X(:,3) == cc,1);
                                if length(this_rates) >= conditional_variance_min_reps
                                    conditional_variance_matrix(cc,hh) = var(this_rates);
                                end
                            end
                        end
                        
                        conditional_variance_mean = nanmean(conditional_variance_matrix(:));
                        
                        FIDora_conditional_slope2overVar(i,tt,1,k) = conditional_sensory_slope^2/conditional_variance_mean;
                        FIDora_conditional_choiceslope2overVar(i,tt,1,k) = conditional_choice_slope^2/conditional_variance_mean;
                    end
                end
                
                % Use flipped partial corr (Fig.6 of Zaidel 2017)
                if group_result(find_for_partial(i)).PREF_PSTH == LEFT  % Flip according to PREF of this cell
                    FIDora_partial_corr_all_flipped(i,:,:,:) = - FIDora_partial_corr_all(i,:,:,:);
                end
                
                progressbar(i/sum(select_for_partial));
            end
            
            % Use R^2 (Fig.4 of Zaidel 2017)
            FIDora_partial_corr_all_Rsquare = FIDora_partial_corr_all.^2;
        end
        
        % ========  Plot partial correlation over time (Zaidel 2017, Figure 6)  HH20180821
        set(figure(082215),'name',sprintf('Partial correlation over time, j = %g, "any significant" out of N = %g, "%s" cells',...
                    j,sum(select_for_partial),t_criterion_txt)); clf;                 
        set(gcf,'uni','norm','pos',[0.016       0.257       0.965       0.528]);
        
        hhcc_text = {'Heading','Choice'};
        
        for hhcc = 1:2
            
            for k = 1:3
                % Flip
                SeriesComparison(squeeze(FIDora_partial_corr_all_flipped(:,:,:,k)), FIDora_tCenters,...
                    'SEM',1,'Errorbar',2,'Axes',subplot(2,5,k),'Transparent',transparent,'colors',[colors(k,:); 0 0 0]);         legend off;

                % Gaussian vel
                axis tight; plot([0 0],ylim,'k--'); xlim([-100 1600]); plot(xlim,[0 0],'k--'); % ylim([0 max(ylim)*1.1])
                plot([1500 1500], ylim,'k--')                
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);

                % R^2
                SeriesComparison(squeeze(FIDora_partial_corr_all_Rsquare(:,:,:,k)), FIDora_tCenters,...
                    'SEM',1,'Errorbar',2,'Axes',subplot(2,5,5 + k),'Transparent',transparent,'colors',[colors(k,:); 0 0 0]);         legend off;

                % Gaussian vel
                axis tight; plot([0 0],ylim,'k--'); xlim([-100 1600]); plot(xlim,[0 0],'k--'); % ylim([0 max(ylim)*1.1])
                plot([1500 1500], ylim,'k--')                
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
            
            end
            
            SeriesComparison(squeeze(FIDora_partial_corr_all_flipped(:,:,hhcc,:)), FIDora_tCenters,...
                'SEM',1,'Errorbar',2,'Axes',subplot(2,5,3+hhcc),'Transparent',transparent,'colors',colors);     legend off;

            title([hhcc_text{hhcc} ' partial R flip'])

            % Gaussian vel
            xlim([-100 1600]); plot(xlim,[0 0],'k--'); % ylim([0 max(ylim)*1.1])
            plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);

            
            SeriesComparison(squeeze(FIDora_partial_corr_all_Rsquare(:,:,hhcc,:)), FIDora_tCenters,...
                'SEM',1,'Errorbar',2,'Axes',subplot(2,5,5+3+hhcc),'Transparent',transparent,'colors',colors);  legend off;
            title([hhcc_text{hhcc} ' partial R^2'])
            
            %         sumFisherMean = squeeze(mean(sumBoots,1));
            %         sumFisherVestPlusVIs = sumFisherMean(:,1) + sumFisherMean(:,2);
            %         plot(CP_ts{1}(plotRange),sumFisherVestPlusVIs,'m','linew',2);
            
            % Gaussian vel
            axis tight; plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
            xlim([-100 1600]); plot(xlim,[0 0],'k--'); % ylim([0 max(ylim)*1.1])
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
            
        end
        SetFigure(15);
        
        figure(082217); clf
        set(gcf,'uni','norm','pos',[0.053       0.422       0.605       0.387]);
        
        % -------- Plotting FI(t) ---------
        bootN = 2000;
        
        % Get sum of Fisher and std by bootstrap (Gu 2010)
        h = subplot(1,2,1);
        sumBoots = bootstrp(bootN,@(x)sum(x,1),squeeze(FIDora_conditional_slope2overVar));
        sumBoots = reshape(sumBoots,bootN,[],3);
          
        SeriesComparison(sumBoots, FIDora_tCenters,...
            'SEM',0,'Errorbar',2,'axes',h,  'Transparent',transparent,'colors',colors);     legend off;
      
        sumFisherMean = squeeze(nanmean(sumBoots,1));
        % sumFisherVestPlusVIs = sumFisherMean(:,1) + sumFisherMean(:,2);
        sumFisherVestPlusVIs = sumFisherMean(:,1) - nanmean(sumFisherMean(FIDora_tCenters < 200, 1)) ...
                             + sumFisherMean(:,2) - nanmean(sumFisherMean(FIDora_tCenters < 200, 2)) ...
                             + nanmean(sumFisherMean(FIDora_tCenters < 200, 3));

        plot(FIDora_tCenters,sumFisherVestPlusVIs,'m','linew',2);
        
        % Gaussian vel
        xlim([-100 1600]); ylim([0 5000]); % ylim([0 max(ylim)*1.1])
        plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
        plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
        legend off;
        title('Partial FI,Sensory')

        % ----- Partial Choice FI ----
        
        h = subplot(1,2,2);
        sumBoots = bootstrp(bootN,@(x)sum(x,1),squeeze(FIDora_conditional_choiceslope2overVar));
        sumBoots = reshape(sumBoots,bootN,[],3);
          
        SeriesComparison(sumBoots, FIDora_tCenters,...
            'SEM',0,'Errorbar',2,'axes',h, 'Transparent',transparent,'colors',colors);     legend off;
      
        sumFisherMean = squeeze(nanmean(sumBoots,1));
        % sumFisherVestPlusVIs = sumFisherMean(:,1) + sumFisherMean(:,2);
        sumFisherVestPlusVIs = sumFisherMean(:,1) - nanmean(sumFisherMean(FIDora_tCenters < 200, 1)) ...
                             + sumFisherMean(:,2) - nanmean(sumFisherMean(FIDora_tCenters < 200, 2)) ...
                             + nanmean(sumFisherMean(FIDora_tCenters < 200, 3));

        plot(FIDora_tCenters,sumFisherVestPlusVIs,'m','linew',2);
        
        % Gaussian vel
        xlim([-100 1600]); ylim([0 max(ylim)*1.1])
        plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
        plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
        legend off;
        title('Partial FI, Choice')
        SetFigure(15)
    end



    % 2. Decoder of heading (not choice/modality) ====  HH20170810   
    function f6p5p9(debug)
        if debug;  dbstack;  keyboard;  end
 
        % Override
        %         min_reps = 15; % Anne used 20
        %         training_reps = 25; % I artificially added some trials by bootstrapping if reps < training_reps.
        %                             % This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
        
        SVM_training_epoch = [0 1700];    % This make it comparable with PCA results
        
        find_for_SVM = find(select_for_SVM);
        
        % Regroup pseudo-trials pool
        pseudo_trial_pool_training = cell(length(find_for_SVM),4);
        pseudo_trial_pool_testing = cell(length(find_for_SVM),4,length(SVM_testing_tcenters));
        
        for ii = 1:length(find_for_SVM)

            this_cell_id = find_for_SVM(ii);
            
            % ====== Get unfiltered binary spike train data with all (correct and wrong) trials to decode heading!! ======
            % Get data
            spike_aligned_in_bin = group_result(this_cell_id).mat_raw_PSTH.spike_aligned{1,j};
            t_in_bin = group_result(this_cell_id).mat_raw_PSTH.spike_aligned{2,j};
            SVM_training_epoch_ind = min(SVM_training_epoch) <= t_in_bin & t_in_bin <= max(SVM_training_epoch);
            
            raw_this = spike_aligned_in_bin(:,SVM_training_epoch_ind);
            
            % ==== Get trial conditions ====
            stim_type_per_trial = group_result(this_cell_id).mat_raw_PSTH.stim_type_per_trial';
            heading_per_trial = group_result(this_cell_id).mat_raw_PSTH.heading_per_trial';
            unique_heading = unique(heading_per_trial);
            % choice_per_trial = group_result(this_cell_id).mat_raw_PSTH.choice_per_trial; % Not important
            % correct_or_zero_per_trial = (choice_per_trial == ((heading_per_trial>0) + 1))  |(heading_per_trial == 0); % Correct or 0 heading; % Not important
            % pref_null = [group_result(this_cell_id).PREF_PSTH LEFT+RIGHT-group_result(this_cell_id).PREF_PSTH]; % Pref goes first % Not important
            
                        
            % Training pool (one certain time window)
            means_this = cellfun(@(x)mean(x(:,SVM_training_epoch_ind),2),raw_this,'uniformOutput',false);
            [pseudo_trial_pool_training{ii,:}] = means_this{:};
            
            % Testing pool (shifting time windows)
            for ttt = 1:length(SVM_testing_tcenters)
                test_epoch_ind = SVM_testing_tcenters(ttt) - SVM_testing_span/2 <= rate_ts{j_for_SVM} & rate_ts{j_for_SVM} <= SVM_testing_tcenters(ttt) + SVM_testing_span/2;
                means_this = cellfun(@(x)mean(x(:,test_epoch_ind),2),raw_this,'uniformOutput',false);
                
                [pool_testing{ii,:,ttt}] = means_this{:};
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
        if ver_num < 2014
            if matlabpool('size') == 0 ;  matlabpool;  end
        else
            if isempty(gcp('nocreate')); parpool; end
        end
        
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
            svm_training_choice(nn) = svmtrain(firing_for_this_training,answers_choice,'box',1e-5,'tol',1e-7,'autoscale',0);
            svm_training_modality(nn) = svmtrain(firing_for_this_training,answers_modality,'box',1e-5,'tol',1e-7,'autoscale',0);
            
            parfor_progress;
        end
        parfor_progress(0);
        toc; disp('Done!');
        
        %% Averaged weights (bagging)
%         weights_svm_choice_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_choice,'uniform',0));
%         weights_svm_modality_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_modality,'uniform',0));

        % By default, 'autoscale' = true, so here I should rescale back the weight! HH20170613
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
            mean(time_markers{j}(1,3)) - group_result(representative_cell).mat_raw_PSTH.binSize_CP/2   % Pre-sac epoch
            mean(time_markers{j}(1,3)) + group_result(representative_cell).mat_raw_PSTH.binSize_CP/2];  % Post-sac epoch
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,center_t_ind] = min(abs(rate_ts{j} - tuning_centers(1)));
        [~,pre_t_ind] = min(abs(rate_ts{j} - tuning_centers(2)));
        [~,post_t_ind] = min(abs(rate_ts{j} - tuning_centers(3)));
        
        tuning_time_phase = [center_t_ind pre_t_ind post_t_ind];
        tuning_time_phase_title = {'Stimulus', 'Pre-sac', 'Post-sac'};
        
        unique_heading = group_result(representative_cell).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
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
        
        % --- Venn diagram for cell counters -- HH20160915
        % Fig.A: Choice preferences for three modalities
        venn_HD_1 = select_bottom_line_all_monkey & group_ChoicePreference_pvalue(1,:,3)' < 0.05;
        venn_HD_2 = select_bottom_line_all_monkey & group_ChoicePreference_pvalue(2,:,3)' < 0.05;
        venn_HD_3 = select_bottom_line_all_monkey & group_ChoicePreference_pvalue(3,:,3)' < 0.05;
        venn_HD_12 = venn_HD_1 & venn_HD_2;
        venn_HD_13 = venn_HD_1 & venn_HD_3;
        venn_HD_23 = venn_HD_2 & venn_HD_3;
        venn_HD_123 = venn_HD_1 & venn_HD_2 & venn_HD_3;
        
        fprintf('HD_1, HD_2, HD_3: %s\n', num2str([sum(venn_HD_1),sum(venn_HD_2),sum(venn_HD_3),...
            sum(venn_HD_12),sum(venn_HD_13),sum(venn_HD_23),sum(venn_HD_123)]));

        figure(9997);  clf
        [H,S]=venn([sum(venn_HD_1),sum(venn_HD_2),sum(venn_HD_3)],...
            [sum(venn_HD_12),sum(venn_HD_13),sum(venn_HD_23),sum(venn_HD_123)],...
            'FaceColor',{'none','none','none'},'Edgecolor',{colors(1,:),colors(2,:),colors(3,:)});
        for i = 1:length(S.ZoneCentroid)
            text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(S.ZonePop(i)));
        end
        title(sprintf('HD_1, HD_2, HD_3: %s, total = %g',num2str([S.CirclePop S.IntersectPop]),sum(S.ZonePop)));
        
        % Fig.B: Memsac and HD
        venn_HD_any = venn_HD_1 | venn_HD_2 | venn_HD_3;
        venn_HD_all = venn_HD_1 & venn_HD_2 & venn_HD_3;
        venn_Mem = select_bottom_line_all_monkey & (group_MemSac_ps(:,3) < 0.05 | group_TwoPtMemSac_ps(:,3) < 0.05);
        venn_HD_any_all = venn_HD_any & venn_HD_all;
        venn_HD_any_Mem = venn_HD_any & venn_Mem;
        venn_HD_all_Mem = venn_HD_all & venn_Mem; 
        venn_HD_any_all_Mem = venn_HD_any & venn_HD_all & venn_Mem;
        
        fprintf('HD_any, HD_all, Memsac: %s\n', num2str([sum(venn_HD_any),sum(venn_HD_all),sum(venn_Mem),...
            sum(venn_HD_any_all),sum(venn_HD_any_Mem),sum(venn_HD_all_Mem),sum(venn_HD_any_all_Mem)]));
        fprintf('HD_any, Memsac: %s\n', num2str([sum(venn_HD_any),sum(venn_Mem),...
            sum(venn_HD_any_Mem)]));
        
%         figure(9996);  clf 
%         [H,S]=venn([sum(venn_HD_any),sum(venn_HD_all),sum(venn_Mem)],...
%             [sum(venn_HD_any_all),sum(venn_HD_any_Mem),sum(venn_HD_all_Mem),sum(venn_HD_any_all_Mem)],...
%             'FaceColor',{'none','none','none'},'Edgecolor',{colors(1,:),colors(2,:),colors(3,:)});
%         for i = 1:length(S.ZoneCentroid)
%             text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(S.ZonePop(i)));
%         end
%         title(sprintf('%s, total = %g',num2str([S.CirclePop S.IntersectPop]),sum(S.ZonePop)));
        
        % Fig.C: Memsac and HD (simple)
        figure(9995);  clf 
        [H,S]=venn([sum(venn_HD_any),sum(venn_Mem)],...
            [sum(venn_HD_any_Mem)],...
            'FaceColor',{'none','none'},'Edgecolor',{colors(1,:),colors(2,:)}); 
        for i = 1:length(S.ZoneCentroid)
            text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(S.ZonePop(i)));
        end
        title(sprintf('HD_any, Memsac: %s, total = %g',num2str([S.CirclePop S.IntersectPop]),sum(S.ZonePop)));
        
        
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

    function f9p3(debug)    % More analysis of behavioral data
        if debug;  dbstack;  keyboard;  end
        
        % Find unique file names
        file_names = {group_result.fileID};
        same_as_previous = [0 strcmp(file_names(1:end-1),file_names(2:end))]; % Add a 0 to the first element
        unique_file_ind = find((~ same_as_previous) & (Psy_pred_ratio > 0.1)'); % Exclude too crazy data due to bad fitting
               
        monkeys = xls_num{1}(:,header.Monkey);  
        monkey{1} = find(monkeys == 5);
        monkey{2} = find(monkeys == 10);
        
        figure(145957); clf;
        set(gcf,'uni','norm','pos',[0.061       0.096       0.854       0.585]);
        
        for mm = 1:2
            
            this_select = intersect(monkey{mm}, unique_file_ind);
            data_to_average = [Psy_thres(this_select,:) Psy_thres_pred(this_select)];
            
            
            mean_threshold = mean(data_to_average);
            sem_threshold = std(data_to_average)/sqrt(size(data_to_average,1));
            
            % Plotting
            
            color_bar = [colors; 0.1 0.1 0.1];
            
            xlim([0.5 4.5]);
            
            subplot(1,2,mm); hold on;
            for i = 1:4
                bar(i,mean_threshold(i),0.7,'facecol',color_bar(i,:),'edgecol','none');
                h = errorbar(i,mean_threshold(i),sem_threshold(i),'color',color_bar(i,:),'linestyle','none','linewidth',3);
                errorbar_tick(h,13);
            end
            
            % Statistics
            [~,p_vest_vis] = ttest(data_to_average(:,1)-data_to_average(:,2));
            [~,p_comb_vest] = ttest(data_to_average(:,1)-data_to_average(:,3));
            [~,p_comb_vis] = ttest(data_to_average(:,2)-data_to_average(:,3));
            [~,p_comb_pred] = ttest(data_to_average(:,3)-data_to_average(:,4));
            
            % Draw lines
            plot(1:4, data_to_average','-','color',[0.5 0.5 0.5]);
            
            set(gca,'xtick',1:4,'xticklabel',{'Vestibular','Visual','Combined','Optimal'});
            rotateXLabels(gca,45);
            ylabel('Threshold');
            
            text(0,max(ylim)*0.8,sprintf('p vest-vis = %g, p comb-vest = %g\np comb-vis = %g, p comb-pred = %g',p_vest_vis,p_comb_vest,p_comb_vis,p_comb_pred),'fontsize',10.5);
            title(sprintf('n = %g unique recording sessions',length(data_to_average)));
            
        end
        SetFigure();
    end


    function f9p4(debug)      % Correlations 3. Psychophysics v.s. Neural activity
        if debug
            dbstack;
            keyboard;
        end
        
       
        % T-cell is important here because Polo's behavior was getting worse while
        % I was getting better at finding T-cells, so...
        select_psycho = select_bottom_line & select_tcells & (Psy_pred_ratio > 0.1); %& ~[zeros(80,1);ones(138-80,1)];  % Exclude too crazy psychocurve
        
        
       %% 1. With Choice divergence
        j = 1;
        
        % Enhancement of cDiv in combined condition
        % enhance_cdiv = max(ChoiceDiv_ModDiffer{1}(:,0 <= rate_ts{j} & rate_ts{j} <= 1500,2),[],2); % Comb - vis
        
        t_begin = 700; t_end = 800;
        enhance_cdiv = nanmean(ChoiceDiv_ModDiffer{1}(:,700 <= rate_ts{j} & rate_ts{j} <= 800, 3 ),2); % Comb - vis
        
        figure(162824); clf;
        
        h = LinearCorrelation({
            Psy_pred_ratio(select_psycho)
            Psy_pred_ratio(select_psycho)
            Psy_pred_ratio(select_psycho)},...
            {
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 1 ),2);
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 2 ),2);
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 3 ),2);
            },...
            'FaceColors',{colors(1,:),colors(2,:),'k'},'Markers',{'o'},...
            'LineStyles',{'b-','r-','k-'},'MarkerSize',marker_size,...
            'Ylabel',sprintf('\\Delta CDiv (%g - %g ms, 3-1, 3-2, 1-2)',t_begin,t_end),'Xlabel','Psycho prediction ratio',...
            'MarkerSize',12,...
            'figN',162824,'XHist',15,'YHist',15,'logx',1,...
            'XHistStyle','stacked','YHistStyle','stacked','Method','Pearson','FittingMethod',2); 
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        title(sprintf('t cells, n = %g',sum(select_psycho)))
        
        %% 2. With Single cell Fisher info
        % Enhancement of Fisher info in combined condition
        % enhance_cdiv = max(ChoiceDiv_ModDiffer{1}(:,0 <= rate_ts{j} & rate_ts{j} <= 1500,2),[],2); % Comb - vis
        
        
        % Which info? (t-cell, not unique session, two monkeys together)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        findForFisher = find(select_psycho);

        t_begin = 500; t_end = 1500;
        FIScatterTimeRange = t_begin < CP_ts{1} & CP_ts{1} <= t_end;
        FIScatterHeadingRange = 1; % +/-8 degrees
        FIScatterVariance = 2; % Real variance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         if isempty(fisherSimpleGu)
            f6p5p1(0); 
%         end
        
        FI_all_temp = squeeze(mean(fisherSimpleGu(:,FIScatterTimeRange,:,FIScatterHeadingRange,FIScatterVariance),2))';
        FI_temp_comb_minus_vest = FI_all_temp(3,:)-FI_all_temp(1,:);
        FI_temp_comb_minus_vis = FI_all_temp(3,:)-FI_all_temp(2,:);
        FI_temp_vest_plus_vis = FI_all_temp(1,:) + FI_all_temp(2,:);
        
        monkeys = xls_num{1}(:,header.Monkey);  
        monkey1 = monkeys == 5;
        monkey2 = monkeys == 10;
        
        [~ , monkey1_in_findForFisher] = intersect(findForFisher, find(monkey1)); % This is nasty, but works
        [~ , monkey2_in_findForFisher] = intersect(findForFisher, find(monkey2));
        
        optimal_ratio_from_Fisher = @(FI) sqrt((FI(2,:)+FI(1,:))./FI(3,:)); % Optimal ratio predicted from Fisher info: sqrt((I1+I2)/I3)
        
        figure(162844); clf;
        h = LinearCorrelation({
            Psy_pred_ratio(select_psycho & monkey1)
            Psy_pred_ratio(select_psycho & monkey2)
%             Psy_pred_ratio(select_psycho)
%             Psy_pred_ratio(select_psycho)
            },...
            {
%             FI_temp_comb_minus_vest';
%             FI_temp_comb_minus_vis';
%             FI_temp_vest_plus_vis';
              optimal_ratio_from_Fisher(FI_all_temp(:,monkey1_in_findForFisher));
              optimal_ratio_from_Fisher(FI_all_temp(:,monkey2_in_findForFisher));
            },...
            'CombinedIndex',3,'PlotCombinedOnly',0,...
            'FaceColors',{'b','r'},... % {colors(1,:),colors(2,:),'k'}
            'Markers',{'o','^'},...
            'LineStyles',{'b:','r:','k'},... % {'b-','r-','k-'},
            'MarkerSize',marker_size,...
            'Ylabel',sprintf('Sqrt((FI_1+FI_2)/FI_3) (%g - %g ms)',t_begin,t_end),'Xlabel','actual optimal ratio (comb/pred)',...
            'MarkerSize',12,... 
            'figN',162844,'XHist',15,'YHist',15,'logx',1,'logy',1,...
            'XHistStyle','grouped','YHistStyle','grouped','Method','Pearson','FittingMethod',2); figN = figN + 1;
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        title(sprintf('t cells, n = %g, %g',sum(select_psycho & monkey1), sum(select_psycho & monkey2)))
        
    end

    function f9p9(debug)      % Export Associated Memsac Files
        if debug;  dbstack;  keyboard;  end
        
        % Read batch file with all cells
        f = fopen('Z:\Data\Tempo\Batch\20150725_BP_allAreas_m5_m10.m');
        line = fgetl(f);  ii = 0;
        
        while line ~= -1
            if line(1) ~= '%'
                ii = ii + 1;
                all_ori_batch{ii} = line;
                tmp = textscan(line,'%s');
                all_ori_name{ii} = tmp{1}{2};
                all_ori_spikeChan(ii) = str2double(tmp{1}{end});
            end
            line = fgetl(f);
        end
        fclose(f);
        
        % HH20180608
        exportPath = 'Z:\Data\Tempo\Batch\ExportedAssociatedMemsacFiles\';
        new_f = fopen([exportPath 'batch_file.m'],'a');
        
        progressbar;
        for ii = 1:length(group_result)
            
            PATH = 'Z:\Data\MOOG\';
            
            if xls_num{1}(ii,header.Monkey) == 5
                monkeyName = 'Polo';
            else
                monkeyName = 'Messi';
            end
            
            if ~isempty(group_result(ii).mat_raw_MemSac)
                FILE = group_result(ii).mat_raw_MemSac.FILE;
                SpikeChan = group_result(ii).mat_raw_MemSac.SpikeChan;
                %                 htbOK = copyfile([PATH monkeyName '\raw\' FILE '.htb'],[exportPath FILE '.htb']);
                %                 logOK= copyfile([PATH monkeyName '\raw\' FILE '.log'],[exportPath FILE '.log']);
                %                 matOK= copyfile([PATH monkeyName '\Analysis\SortedSpikes2\' FILE '.mat'],[exportPath FILE '.mat']);
                
                % Write new batch file
                where = find(strcmp(all_ori_name,FILE) & all_ori_spikeChan == SpikeChan);
                if length(where) ~= 1 % Should be unique
                    error('No batch info?')
                else
                    fprintf(new_f,'%s\n',all_ori_batch{where});
                end
            end
            
            if ~isempty(group_result(ii).mat_raw_2ptMemSac)
                FILE = group_result(ii).mat_raw_2ptMemSac.FILE;
                SpikeChan = group_result(ii).mat_raw_2ptMemSac.SpikeChan;
                %                 htbOK = copyfile([PATH monkeyName '\raw\' FILE '.htb'],[exportPath FILE '.htb']);
                %                 logOK= copyfile([PATH monkeyName '\raw\' FILE '.log'],[exportPath FILE '.log']);
                %                 matOK= copyfile([PATH monkeyName '\Analysis\SortedSpikes2\' FILE '.mat'],[exportPath FILE '.mat']);
                
                % Write new batch file
                where = find(strcmp(all_ori_name,FILE) & all_ori_spikeChan == SpikeChan);
                if length(where) ~= 1 % Should be unique
                    error('No batch info?')
                else
                    fprintf(new_f,'%s\n',all_ori_batch{where});
                end
            end
            
            
            progressbar(ii/length(group_result));
        end
        
        fclose(new_f);

        edit([exportPath 'batch_file.m']);
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
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
            
            % Annotate tcells
            plot(weights(select_tcells(select),1),weights(select_tcells(select),2),'+','markersize',16,'color',colors(2,:),'linew',2);
            
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
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(weights(select_tcells(select),1),abs(Choice_pref_all(k,select_tcells(select),tt)),'+','markersize',16,'color',colors(2,:),'linew',2);
        
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
%             'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
%         
%         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
%         
%         % Annotate tcells
%         plot(MemSac_indicator(select_tcells),weights(select_tcells(select),1),...
%             '+','markersize',16,'color',colors(2,:),'linew',2);
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
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        delete([h.group(1:2).line]);
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),weights(select_tcells(select),1),...
            '+','markersize',16,'color',colors(2,:),'linew',2);
        
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
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
            
            delete([h.group(1:4).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
            
            % Annotate tcells
            plot((Modality_pref_all(k,select_tcells(select),tt)), weights(select_tcells(select),2),...
                '+','markersize',16,'color',colors(2,:),'linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot((Modality_pref_all(k,:,tt)),weights(:,2),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select});
        end
        
    end


    function h = Weighted_sum_PSTH(weights,txt,selects)     % Plot weighted sum of firing rates
        
        % ------  PSTH correct only, all choices -----------
        
        if ~iscell(weights)  ;  weights = {weights};     end
        if ~iscell(selects) ;  selects = {selects};   end
        
        set(figure(figN),'pos',[34   83  584*length(weights)  504*2]); clf; figN = figN + 1;
        
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
            
            
            % Note here the errorbars should be STD instead of SEM. (Not independent sampling, but bootstrapping)
            h_subplot = subplot(2,length(weights),pp);
            h_series = SeriesComparison({PSTH_projected{1} PSTH_projected{2}},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),colors(3,:),colors(3,:)},'LineStyles',{'-','--','-','--','-','--'},...
                'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Weighted sum of firing','axes',h_subplot);
            hold on;    legend off;     axis tight;
            xlim([min(xlim) + 200 max(xlim)-300]);   % ylim([0.1 .7]);
            
            h_subplot = subplot(2,length(weights),pp+length(weights));
            h_series = SeriesComparison({PSTH_projected{1}(:,:,1:2:end)-PSTH_projected{1}(:,:,2:2:end) PSTH_projected{2}(:,:,1:2:end)-PSTH_projected{2}(:,:,2:2:end)},...
                {rate_ts{1} rate_ts{2} time_markers},...
                'Colors',{colors(1,:),colors(2,:),colors(3,:)},'LineStyles',{'-'},...
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
        h_subplot = tight_subplot(3,2,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
        
        for k = 1:3
            
            for j = 1:2
                
                h.projected_angles_allbootstrap{j} = nan(nbootstrap,size(PSTH_correct_angles_raw{j},2),size(PSTH_correct_angles_raw{j},3),size(PSTH_correct_angles_raw{j},4));
                h.projected_angles_diff{j} = nan(nbootstrap,size(PSTH_correct_angles_raw{j},2),size(PSTH_correct_angles_raw{j},3)/2,size(PSTH_correct_angles_raw{j},4));
                
                for bb = 1:nbootstrap
                    
                    % ------- Different angles -------
                    % Deal with nans of 0 heading for some cells
                    PSTH_correct_angles_raw_this = PSTH_correct_angles_raw{j}(selects{1},:,:,k);
                    yyy = nan(size(PSTH_correct_angles_raw_this,2),size(PSTH_correct_angles_raw_this,3));
                    
                    % HH20170808
                    for hh = 1:size(yyy,2) % For each heading
                        non_nan_this_heading = ~(any(isnan(PSTH_correct_angles_raw_this(:,:,hh)),2));
                        yyy(:,hh) = PSTH_correct_angles_raw_this(non_nan_this_heading,:,hh)' ...
                            * weights{1}(non_nan_this_heading,bb) / sum(weights{1}(non_nan_this_heading,bb));
                    end

                    %{ 
                    % HH20170808
                    PSTH_correct_angles_raw_this(isnan(PSTH_correct_angles_raw_this)) = 0;
                    
                    % Weighted sum with choice weight
                    yyy = reshape(reshape(PSTH_correct_angles_raw_this,size(PSTH_correct_angles_raw_this,1),[])' ...
                        * weights{1}(:,bb) / sum(weights{1}(:,bb)), size(PSTH_correct_angles_raw_this,2),[]);
                    %}
                    
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
                'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Diff','axes',h_subplot(k+3));
            
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
        
        h_marker = guidata(gcbo);
                
        allX = get(h_line,'xData');
        allY = get(h_line,'yData');
        n_single_group = length(allX)/couples;
        
        % ------- Recover the cell number -------
        
        if ismember('control',get(gcf,'currentModifier'))  % Select cell from file name and show in the figure. @HH20150527
            fileN = input('Which cell do you want from the figure?    ','s');
            available_cells = find(select_for_this);
            
            if fileN(1) == '#'  % HH20160905. Direct input the original cell #
                ori_cell_no = str2double(fileN(2:end));
                ind = sum(select_for_this(1:ori_cell_no));
            else
                
                nn = 1; iffind = [];
                while nn <= length(available_cells) && isempty(iffind) % Find the first match
                    iffind = strfind(group_result(available_cells(nn)).cellID{1}{1},fileN);
                    nn = nn + 1;
                end
                
                if ~isempty(iffind)
                    ind = nn - 1;
                else
                    fprintf('Are you kidding...?\n');
                    return
                end
                
                ori_cell_no = find(cumsum(select_for_this)==ind,1);
            end
        else   % Select cell from figure
            pos = get(gca,'currentPoint'); posX = pos(1,1); posY = pos(1,2);
            [min_dis,ind] = min(abs(((posX-allX)/range(xlim)).^2+((posY-allY)/range(ylim)).^2));
            if min_dis > (range(xlim)^2+range(ylim)^2)/100 +inf ; return; end
            
            ind = mod(ind-1,n_single_group) + 1; % Deal with coupled cells
            ori_cell_no = find(cumsum(select_for_this)==ind,1);
        end
        
        % Plotting
        if ~isempty(h_marker) ; try delete(h_marker); catch ; end ;end
        
        all_inds = mod(1:length(allX),n_single_group) == mod(ind,n_single_group); % Deal with coupled cells
        h_marker = plot(allX(all_inds),allY(all_inds),'x','color','m','markersize',15,'linew',3);
        
        % Do plot
        Plot_HD([],[],ori_cell_no);
        
        guidata(gcbo,h_marker);
    end
   
    function Plot_HD(~,~,ori_cell_no,h_1463_axes)    % Plot associated HD traces @HH20150426
        
       
        % ------ PSTH in HD ------
        j_this = 1;
        
        for j = 1:2
            ys_this{j} = group_result(ori_cell_no).mat_raw_PSTH.PSTH{j,1,1}.ys';
            sem_this{j} = group_result(ori_cell_no).mat_raw_PSTH.PSTH{j,1,1}.sem';
            ps_this{j} = group_result(ori_cell_no).mat_raw_PSTH.PSTH{j,1,1}.ps';
        end
        
        % ------ Adapted for batch plot example cells ---
        if nargin < 4
            figure(1463); clf ;
            set(gcf,'uni','norm','pos',[0.005       0.597       0.992       0.304]);
            h_1463_axes = tight_subplot(1,5,[0.1 0.06],[0.2 0.15],[0.05 0.05]);
        end
        
        % --- 1. Raw PSTHs ---
        % h_1463_PSTH = subplot_tight(1,5,1,[0.1 0.06],[0.2 0.15 0.05 0.05]);       
        h_1463_PSTH = h_1463_axes(1);
        axes(h_1463_PSTH)
        
        SeriesComparison({shiftdim(ys_this{1},-1) shiftdim(ys_this{2},-1)},...
            {rate_ts{1} rate_ts{2} time_markers},...
            'OverrideError',{sem_this{1}, sem_this{2}},...
            'OverridePs',{ps_this{1}, ps_this{2}},'ErrorBar',6,...
            'CompareIndex',[1 3 5;2 4 6],'CompareColor',{colors(1,:),colors(2,:),colors(3,:)},...
            'Colors',{colors(1,:),colors(1,:),colors(2,:),colors(2,:),colors(3,:),colors(3,:)},...
            'Transparent',transparent,'LineStyles',{'-','--'},'axes',h_1463_PSTH);
        
        axis tight; legend off;
        xlim([-300 2300]);        
        xticks([0 1000 2050]); xticklabels({'Stim','1000 ms','Sac'});
        
        
        if group_result(ori_cell_no).Waveform_broad
            celltype = 'Pyr';
        else
            celltype = 'Int';
        end
        
        
        id_this = group_result(ori_cell_no).cellID{1}{1};
        name_this = id_this(strfind(id_this,'u'):end);
        pos_this = group_result(ori_cell_no).position;
        
        hemi = 'L' * (pos_this(2) == 1) + 'R' * (pos_this(2) == 2);
                
        title(sprintf('#%g, %s, %s\n%s: AP%gML%.2gD%g',...
            ori_cell_no,celltype, name_this, hemi, pos_this(3:5)));
        
        % xlabel('Time to saccade onset (ms)');      ylabel('Firing rate');
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),0 + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
        
        % --- 2. PSTH_difference ---
        h_1463_PSTH_diff = h_1463_axes(2);
        axes(h_1463_PSTH_diff)

        for j = 1:2
            psth_diff{j} = ys_this{j}(:,[1 3 5]) - ys_this{j}(:,[2 4 6]);
            psth_diff_sem{j} = sqrt(sem_this{j}(:,[1 3 5]).^2 + sem_this{j}(:,[2 4 6]));
        end
        
        SeriesComparison({shiftdim(psth_diff{1},-1) shiftdim(psth_diff{2},-1)},...
            {rate_ts{1} rate_ts{2} time_markers},...
            'OverrideError',{psth_diff_sem{1}, psth_diff_sem{2}},...
            'CompareIndex',[1 2 3;1 2 3],'CompareColor',{colors(1,:),colors(2,:),colors(3,:)},...
            'Colors',{colors(1,:),colors(2,:),colors(3,:)},...
            'Transparent',transparent,'LineStyles',{'-'},'axes',h_1463_PSTH_diff);
        
        axis tight; legend off;
        xlim([-300 2300]);
        xticks([0 1000 2050]); xticklabels({'Stim','1000 ms','Sac'});

        % Gaussian vel
        plot(xlim,[0 0], 'k--')
        % ylim([max(-20,min(ylim)) max(ylim)]);
        
        posi_neg_ratio = max(5, max(ylim)/abs(min(ylim)));
        ylim([-max(ylim)/posi_neg_ratio max(ylim)]);
        
        plot(Gauss_vel(:,1) + time_markers{1}(1), 0 + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);

        title(sprintf('\\DeltaPSTH, p = %.3g,%.3g,%.3g',...
            group_result(ori_cell_no).ChoicePreference_pvalue(3,:)'));
      
        % --- 3. Choice Divergence with newly added p value ---
        h_1463_CD = h_1463_axes(3);
        axes(h_1463_CD)
        
        for j = 1:2
            CD_this{j} = ChoiceDiv_All{j}(ori_cell_no,:,:);
            CD_perm_this{j} = 1.96 * ChoiceDiv_All_perm{j}.std(ori_cell_no,:,:);
            ps_CD_this{j} = squeeze(ChoiceDiv_All_perm{j}.p(ori_cell_no,:,:));
        end
        
        ylim([-0.2 1.1]);
         
        SeriesComparison({CD_this{1}, CD_this{2}},...
            {rate_ts{1} rate_ts{2} time_markers},...
            'Colors',{colors(1,:),colors(2,:),colors(3,:)},...
            'OverridePs',{ps_CD_this{1}, ps_CD_this{2}},'ErrorBar',4,...
            'CompareIndex',[1 2 3 ; 1 2 3],'CompareColor',{colors(1,:),colors(2,:),colors(3,:)},...
            'Transparent',transparent,'LineStyles',{'-'},'axes',h_1463_CD, 'Ylim',[-0.2 1]);
 
        %         hh = SeriesComparison({CD_perm_this{1}, CD_perm_this{2}},...
        %             {rate_ts{1} rate_ts{2} time_markers},...
        %             'Colors',{colors(1,:),colors(2,:),colors(3,:)},...
        %             'LineStyles',{'-'},'figN',1467,'hold',1);
        %         set([hh.h{:}],'linew',0.5)
        legend off; plot(xlim,[0 0],'k--');
        
        title('CD');
        xlim([-300 2300]);
        % ylim([-0.2 1.1]);
        ylim([-max(ylim)/posi_neg_ratio max(ylim)]); % Same as deltaPSTH
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        xticks([0 1000 2050]); xticklabels({'Stim','1000 ms','Sac'});
        yticks([0 0.5 1]);  yticklabels(num2str([0 0.5 1],'%0.1f\n'));
        
        
        % -- PSTH_diff (Hard and Easy) -- HH20160905
%         %{
        figure(1466); clf;
        set(gcf,'uni','norm','pos',[0.013       0.099       0.286       0.352]);
        
        for j = 1:2
            PSTH_diff_hardeasy{j} = PSTH_hard_easy_raw_cellBtB4Bk{j}(ori_cell_no,:,1:2:end,:)...
                - PSTH_hard_easy_raw_cellBtB4Bk{j}(ori_cell_no,:,2:2:end,:);
        end
        
        SeriesComparison({reshape(PSTH_diff_hardeasy{1},1,[],6) reshape(PSTH_diff_hardeasy{2},1,[],6)},...
            {rate_ts{1} rate_ts{2} time_markers},...
            'Colors',hsv2rgb([2/3 .4 1; 2/3 1 1; 0 .4 1; 0 1 1; 1/3 .4 1; 1/3 1 1]),...
            'LineStyles',{'-','-'},'figN',1466);
        title('PSTH\_diff, Hard and Easy');
        axis tight; legend off; plot(xlim,[0 0],'k--');
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        xlim([-300 2300]);
        %}
        
        % -- Dora Tuning -- HH20170327 @ UNIGE
        if ~isempty(dora_tuning_mean_each_cell)
            
            ind_in_partial = sum(select_for_partial(1:ori_cell_no)); % From original_ind to index in select_partial
            unique_heading = group_result(representative_cell).unique_heading;

            set(figure(1465),'name','Dora tuning'); clf;
            set(gcf,'uni','norm','pos',[0.425       0.057       0.567       0.607]);
            to_plot_tt = [1 4 5 6];
            [h_sub,~] = tight_subplot(3,length(to_plot_tt),[0.05 0.02]);
            zero_index = find(unique_heading == 0);
            
            for tt = 1:length(to_plot_tt)
                for k = 1:3
                    axes(h_sub(k+(tt-1)*3));
                    
                    hh_index = {1:zero_index zero_index:length(unique_heading)};
                    cc_marker = {'<','>'};
                    
                    for cc = LEFT:RIGHT
                        for hh = LEFT:RIGHT
                            %                 plot(unique_heading(1:zero_index)-0.2,tuning_mean_dora(1,1:zero_index,tuning_time_phase(pp),k),'v-','markersize',9,'color',colors(k,:),'markerfacecolor',colors(k,:),'LineWid',2);
                            %                 plot(unique_heading(zero_index:end)+0.2,tuning_mean_dora(2,zero_index:end,tuning_time_phase(pp),k),'^-','markersize',9,'color',colors(k,:),'markerfacecolor',colors(k,:),'LineWid',2);
                            %
                            %                 plot(unique_heading(1:zero_index)+0.2,tuning_mean_dora(2,1:zero_index,tuning_time_phase(pp),k),'^-','markersize',9,'color',colors(k,:),'markerfacecolor','none','LineWid',2);
                            %                 plot(unique_heading(zero_index:end)-0.2,tuning_mean_dora(1,zero_index:end,tuning_time_phase(pp),k),'v-','markersize',9,'color',colors(k,:),'markerfacecolor','none','LineWid',2);
                            
                            h = scatter(unique_heading(hh_index{hh}),...
                                dora_tuning_mean_each_cell(ind_in_partial,to_plot_tt(tt),k,cc,hh_index{hh}),...
                                dora_tuning_n_each_cell(ind_in_partial,to_plot_tt(tt),k,cc,hh_index{hh})*20,...
                                colors(k,:),[cc_marker{cc}]); hold on;
                            
                            if cc == hh
                                set(h,'markerfacecolor',colors(k,:));
                            else
                                set(h,'markerfacecolor','none');
                            end
                        end
                        
                        plot(unique_heading,squeeze(dora_tuning_mean_each_cell(ind_in_partial,to_plot_tt(tt),k,cc,:)),'color',colors(k,:),'LineWid',2)
                    end
                    
                    if k == 1
                        title(partial_corr_timewins{to_plot_tt(tt),2});
                    end

                end
            end
        end
        
        % -- Plot visual PSTH grouped by heading and choice (correct only) for Yong Gu -- HH20170327 @ UNIGE
        %{
        set(figure(1466),'name','Yong Gu needs this but I don''t see the point'); clf;
        set(gcf,'uni','norm','pos',[0.541       0.526       0.287       0.353]);
        
        j = 1; k = 2;
        colors_angles_hsv(:,2) = linspace(0.1,1,5);
        colors_angles_hsv(:,[1 3]) = 1;
        colors_angles = hsv2rgb(colors_angles_hsv);
        colors_angles = reshape(repmat(colors_angles,1,2)',3,[])';
        colors_angles = mat2cell(colors_angles,ones(10,1));
        
        % --- Ramping with different angles ---
        for j = 1:2
            yyy{j} = PSTH_correct_angles_raw{j}(ori_cell_no,:,:,k);
            ttt{j} = rate_ts{j};
        end
        
        h = SeriesComparison(yyy,{ttt{1} ttt{2} time_markers},...
            'Colors',colors_angles,'LineStyles',{'-','--'},...
            'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','figN',1466);
        legend off;
        xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  ylim([0.1 max(ylim)]);
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        %}
        
        Plot_memsac([], [], ori_cell_no, h_1463_axes);
    end

    function Plot_memsac(~, ~, ori_cell_no, h_1463_axes)    % Plot associated mem-sac traces @HH20150426
        
        %         set(figure(1462),'pos',[-1243 452 712 394]); clf;
        %         if ishandle(1462)
        %             figure(1462); clf;
        %         else
        %             set(figure(1462), 'unit','norm','pos',[0.00366032210834553 0.0572916666666667 0.418740849194729 0.385416666666667]);
        %         end
        %
        %                 axis_typing = axes('Position',[0.032 0.672 0.459 0.296]);
        %                 axis_polar = axes('Position',[0.039 0.096 0.345 0.586]);
        %                 axis_psth = axes('Position',[0.521 0.17 0.426 0.728]);
        
        h_1463_memPSTH = h_1463_axes(4);
        h_1463_memPolar = h_1463_axes(5);
        
        try
            
            if ~isempty(group_result(ori_cell_no).mat_raw_MemSac)
                
                % plot_time_range = [2 3 4 6]; % 2:end
                plot_time_range = [3]; % Memory
                
                resp_mean = group_result(ori_cell_no).mat_raw_MemSac.resp_mean(plot_time_range);
                resp_err = group_result(ori_cell_no).mat_raw_MemSac.resp_err(plot_time_range);
                p = group_result(ori_cell_no).mat_raw_MemSac.p(plot_time_range);
                vectSum = group_result(ori_cell_no).mat_raw_MemSac.vectSum(plot_time_range);
                temporal_Slice = group_result(ori_cell_no).mat_raw_MemSac.temporal_Slice(plot_time_range,:);
                align_offsets = group_result(ori_cell_no).mat_raw_MemSac.align_offsets;
                align_markers = group_result(ori_cell_no).mat_raw_MemSac.align_markers;  % Desired markers: target onset & target offset & sac onset
                unique_heading = group_result(ori_cell_no).mat_raw_MemSac.unique_heading;
                
                %             PREF_target_location = group_result(ori_cell_no).mat_raw_PSTH.PREF_target_location; % Note this would be [-90,270]! HH20160906
                %             PREF_target_location = mod(PREF_target_location,360); % Changed to [0,360] HH20160906
                                
                % Memsac PSTH: align to VisOn and SacOn. HH20180609
                MemSac_interp_PSTH = group_result(ori_cell_no).mat_raw_MemSac.MemSac_interp_PSTH([1 3]);
                MemSac_interp_PSTH_sem = group_result(ori_cell_no).mat_raw_MemSac.MemSac_interp_PSTH_sem([1 3]);
                MemSac_interp_PSTH_ts = group_result(ori_cell_no).mat_raw_MemSac.t_centers([1 3]);

                mem_sac_align_offset = mean(align_offsets); % This is time of [VisON, VisON, Saccade]
                mem_sac_time_marker = {mem_sac_align_offset - mem_sac_align_offset(1), mem_sac_align_offset - mem_sac_align_offset(3)};
                mem_sac_border = [mem_sac_align_offset(2)-100 -300]; % for SeriesComparison
                mem_sac_lims = [-300 mem_sac_align_offset(2)+700];
                
                MemSac_interp_locations = group_result(ori_cell_no).MemSac_interp_locations;
                
                [~,pref] = min(abs(group_PREF_target_location(ori_cell_no) - MemSac_interp_locations));
                %                     [~,pref] = min(abs(group_PREF_target_location_notebook(ori_cell_no) - MemSac_interp_locations));
                pref = mod(pref-1,length(MemSac_interp_locations)-1)+1;
                null = mod(pref + (length(MemSac_interp_locations)-1)/2 -1, length(MemSac_interp_locations)-1)+1; % The opposite position
                
                % Polar plot
                [~, polarOrder] = sort(max(cell2mat(resp_mean'),[],2),1,'descend'); % The largest goes first in the polar plot
                
                axes(h_1463_memPolar);
                title_text = '';
                
                for sliceN = polarOrder(:)'
                    
                    if strcmp(temporal_Slice{sliceN,5},''); continue; end
                    
                    if p(sliceN) < 0.05
                        sigMarker = '-'; wid = 2;
                    else
                        sigMarker = '-'; wid = 1;
                    end;
                    
                    polarwitherrorbar([unique_heading/180*pi; unique_heading(1)/180*pi], ...
                        [resp_mean{sliceN}'; resp_mean{sliceN}(1)], (p(sliceN) < 0.05 || 1 ) * [resp_err{sliceN}' ; resp_err{sliceN}(1)],...
                        [temporal_Slice{sliceN,5} sigMarker], wid, 1);  % Simple mode = 1
                    
                    %                     set(polar([unique_heading/180*pi; unique_heading(1)/180*pi], ...
                    %                         [resp_mean{sliceN}'; resp_mean{sliceN}(1)],[temporal_Slice{sliceN,5} sigMarker]),'linewidth',2);
                    %
                    hold on;
                    
                    h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,...
                        [0  max(ylim)],[temporal_Slice{sliceN,5} '-']);
                    
                    
                    %                     h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,...
                    %                                 [max(cell2mat(resp_mean))*0.9 max(cell2mat(resp_mean))*1.3],[temporal_Slice{sliceN,5} '-']);
                    
                    if p(sliceN) < 0.05
                        %         h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,[0 vectAmp(sliceN)],[temporal_Slice{sliceN,5} '-']);
                        set(h_p,'LineWidth',3);
                    else
                        set(h_p,'LineWidth',1);
                    end
                    %
                    % axes(axis_typing); axis off;
                    %                     text(0,- sliceN * 0.15+1.1,sprintf('%s:  \\itp_ = %4.2g',temporal_Slice{sliceN,4},p(sliceN)),'color',temporal_Slice{sliceN,5},...
                    %                         'FontSize',14 * (p(sliceN) < 0.05) + 7 * (p(sliceN) >= 0.05 || isnan(p(sliceN))));
                    drawnow;
                    title_text = [title_text sprintf('%s:  \\itp_ = %4.2g\n',temporal_Slice{sliceN,4},p(sliceN))];
                end
                
                % Target position
                axes(h_1463_memPolar);
                radi = max(xlim);
                pref_theta = MemSac_interp_locations(pref)/180*pi;
                plot(radi*cos(pref_theta),radi*sin(pref_theta),'ok','markersize',5,'linew',2,'markerfacecol','k');
                
                null_theta = MemSac_interp_locations(null)/180*pi;
                plot(radi*cos(null_theta),radi*sin(null_theta),'ok','markersize',5,'linew',2);
                
                text(min(xlim),min(ylim),title_text);
                
                
                % ----- PSTH plot -----
                
                % (1) The one which is nearest to the vect_sum of Mem-period
                
                %{
                which_is_memory = (plot_time_range == 3);
                
                [~,pref] = min(abs(vectSum(which_is_memory) - MemSac_interp_locations));
                pref = mod(pref-1,length(MemSac_interp_locations)-1)+1;
                null = mod(pref + (length(MemSac_interp_locations)-1)/2 -1, length(MemSac_interp_locations)-1)+1; % The opposite position
                
                for jjj = 1:2
                    ys{jjj} = shiftdim(MemSac_interp_PSTH{jjj}(:,[pref null]),-1);
                    sems{jjj} = MemSac_interp_PSTH_sem{jjj}(:,[pref null]);
                end
                
                axes(h_1463_memPSTH); hold on;
                
                SeriesComparison({ys{1},ys{2}},...
                                 {MemSac_interp_PSTH_ts{1},MemSac_interp_PSTH_ts{2},mem_sac_time_marker},...
                                'OverrideError',{sems{1},sems{2}},...
                                'ErrorBar', 2, 'CompareIndex',[],...
                                'Colors',{colors(3,:),colors(3,:)},...
                                'Border', mem_sac_border,  'LineStyles',{'-','--'},'axes',h_1463_memPSTH);
                
                text(-1000,max(ylim)*.95,'MemsacPref\_interp','color',colors(3,:));
                %}
                
                % (2) The two which are nearest to actual PREF location. @HH20150524
                if isempty(group_result(ori_cell_no).TwoPtMemSac_PSTH) % Only if there is not TwoPtMemSac

                    
                    %                 set(polar([MemSac_interp_locations(pref) MemSac_interp_locations(pref)]/180*pi,[0 max(cell2mat(resp_mean))*1.3],'k-'),'linew',3);
                    %                 set(polar([MemSac_interp_locations(null) MemSac_interp_locations(null)]/180*pi,[0 max(cell2mat(resp_mean))*1.3],'k--'),'linew',3);
                    
                    
                    for jjj = 1:2
                        ys{jjj} = shiftdim(MemSac_interp_PSTH{jjj}(:,[pref null]),-1);
                        sems{jjj} = MemSac_interp_PSTH_sem{jjj}(:,[pref null]);
                    end
                    
                    axes(h_1463_memPSTH); hold on;
                    
                    SeriesComparison({ys{1},ys{2}},...
                        {MemSac_interp_PSTH_ts{1},MemSac_interp_PSTH_ts{2},mem_sac_time_marker},...
                        'OverrideError',{sems{1},sems{2}},...
                        'ErrorBar', 2, 'CompareIndex',[],...
                        'Colors',{'k'},'hold', 1,'Ylim',[0 max(max(max((ys{:}))))], ...
                        'Transparent',transparent,'Border', mem_sac_border, 'LineStyles',{'-','--'},'axes',h_1463_memPSTH);
                    
                    text(0,max(ylim)*.88,'ActualFixP\_interp','color','k');
                end
            end
            
            % (3) The two in 2pt-memsac. @HH20150524
            TwoPtMemSac_PSTH = group_result(ori_cell_no).TwoPtMemSac_PSTH;
            
            if ~isempty(TwoPtMemSac_PSTH) 
                
                % Memsac PSTH 2pt: align to VisOn and SacOn. HH20180609
                MemSac_interp_PSTH_ts = group_result(ori_cell_no).mat_raw_2ptMemSac.t_centers([1 3]);
                
                align_offsets = group_result(ori_cell_no).mat_raw_MemSac.align_offsets;
                mem_sac_align_offset = mean(align_offsets); % This is time of [VisON, VisON, Saccade]
                mem_sac_time_marker = {mem_sac_align_offset - mem_sac_align_offset(1), mem_sac_align_offset - mem_sac_align_offset(3)};
                mem_sac_border = [mem_sac_align_offset(2)-100 -300]; % for SeriesComparison
                mem_sac_lims = [-300 mem_sac_align_offset(2)+700];
               
                plot_aligns = [1 3];
                
                for jjj = 1:2
                    MemSac_interp_PSTH{jjj} = reshape...
                        ([group_result(ori_cell_no).mat_raw_2ptMemSac.result_PSTH_anne_mean{plot_aligns(jjj),:}],...
                        length(MemSac_interp_PSTH_ts{jjj}),[]);
                    MemSac_interp_PSTH_sem{jjj} = reshape...
                        ([group_result(ori_cell_no).mat_raw_2ptMemSac.result_PSTH_anne_sem{plot_aligns(jjj),:}],...
                        length(MemSac_interp_PSTH_ts{jjj}),[]);
                end
                
                % If HD_pref is on the left, flip 2pt-MemSac to let PREF go first
                if 90 <= group_PREF_target_location(ori_cell_no) && group_PREF_target_location(ori_cell_no) <=270
                    for jjj = 1:2
                        MemSac_interp_PSTH{jjj} = fliplr(MemSac_interp_PSTH{jjj});
                        MemSac_interp_PSTH_sem{jjj} = fliplr(MemSac_interp_PSTH_sem{jjj});
                    end
                end
                
                for jjj = 1:2
                    ys{jjj} = shiftdim(MemSac_interp_PSTH{jjj},-1);
                    sems{jjj} = MemSac_interp_PSTH_sem{jjj};
                end
                
                axes(h_1463_memPSTH); hold on;
                
                SeriesComparison({ys{1},ys{2}},...
                    {MemSac_interp_PSTH_ts{1},MemSac_interp_PSTH_ts{2},mem_sac_time_marker},...
                    'OverrideError',{sems{1},sems{2}},...
                    'ErrorBar', 2, 'CompareIndex',[],...
                    'Colors',{'k'}, 'hold', 1, 'Ylim',[0 max(max(max((ys{:}))))],...
                    'Transparent',transparent,'Border', mem_sac_border, 'LineStyles',{'-','--'},'axes',h_1463_memPSTH);
                
                text(0,max(ylim)*.81,'2pt','color',hsv2rgb([0.8 1 0.6]));
                
                %                 temporal_Slice = group_result(ori_cell_no).mat_raw_2ptMemSac.temporal_Slice(2:end,:);
                %                 align_offsets = group_result(ori_cell_no).mat_raw_2ptMemSac.align_offsets;
                %                 align_markers = group_result(ori_cell_no).mat_raw_2ptMemSac.align_markers;  % Desired markers: target onset & target offset & sac onset

            end
            
            % Annotating
            % xlabel('Time to saccade onset (ms)');   
            ylabel('Firing rate (Hz)');
            xlim(mem_sac_lims)
            legend off;
            
            % title(sprintf('PREF: %g, %s',group_PREF_target_location(ori_cell_no),group_result(ori_cell_no).cellID{2}{1}(30:end)));
            title(sprintf('%s',group_result(ori_cell_no).cellID{2}{1}(30:end)));
            
            for sliceN = 1:size(temporal_Slice,1)
                if temporal_Slice{sliceN,3} == 4 % Precise windows
                    plot([temporal_Slice{sliceN,1} temporal_Slice{sliceN,2}],[max(ylim) max(ylim)],temporal_Slice{sliceN,5},'linew',5);
                else  % Windows that is not so precise
                    % Mean shift
                    meanShift = mean(align_offsets(:,align_markers == temporal_Slice{sliceN,3})-align_offsets(:,align_markers==4),1);
                    plot( [meanShift meanShift],ylim,'k--','LineWidth',1);
                    plot(meanShift + [temporal_Slice{sliceN,1} temporal_Slice{sliceN,2}],[max(ylim) max(ylim)],temporal_Slice{sliceN,5},'linew',5);
                end
            end
            
            % Plot target manually
            plot([0 mem_sac_time_marker{1}(2)],[max(ylim) max(ylim)],'r','linew',5);

            xticks([0 1000 sum(abs(mem_sac_border))+100]); xticklabels({'TargOn','1000 ms','Sac'});
            ylim([0 max(ylim)])
            
            SetFigure(13)
            
        catch err
            err
            % keyboard
        end
    end


end

