%%   LIP HD Pooled Data
% Begin @ HH 201406

function Group_HD(XlsData)

try
    
    %% Get data
    num = XlsData.num;
    txt = XlsData.txt;
    raw = XlsData.raw;
    header = XlsData.header;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Batch address
    mat_address = {
        % Major protocol goes here (Address, Suffix)
       'Z:\Data\Tempo\Batch\20150418_LIP_HD_m5_modifiedBatch','PSTH';
        %         'Z:\Data\Tempo\Batch\20150411_LIP_HD_m5','PSTH';
        %     'Z:\Data\Tempo\Batch\20141118_LIP_Decision_WithAUC_ResultPSTHOrderChanged','PSTH';
        % Associative protocols
       'Z:\Data\Tempo\Batch\20150411_LIP_memSac_m5','MemSac';
        };
    
    % Global Mask for the each protocol (for XLS)
    mask_all = {
        strcmp(txt(:,header.Protocol),'HD') & strcmp(txt(:,header.Area),'LIP') & (num(:,header.Monkey) == 5);
        % &(num(:,header.HD_rep) >= 8); % This has been moved to sections below (I will find repN from .mat files)
        strcmp(txt(:,header.Protocol),'MemSac') & strcmp(txt(:,header.Area),'LIP') & (num(:,header.Monkey) == 5);
        };
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
    
    %
    % Load .mat files and put the data into a large structure array "group_result(i).mat_raw"
    group_result(size(xls_txt{1},1)).cellID = [];  % The unique cellID
    
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
            
            
            try
                % Load .mat for major and associative protocols
                if ~isnan(file_i)
                    group_result(major_i).cellID{pp} = cell_info{pp}{file_i};
                    
                    mat_file_name = sprintf('%s_%g',xls_txt{pp}{file_i,header.FileNo},xls_num{pp}(file_i,header.Chan1));
                    mat_file_fullname = [mat_address{pp,1} '\' mat_file_name '_' mat_address{pp,2}];
                    
                    raw = load(mat_file_fullname);
                    group_result(major_i).(['mat_raw_' mat_address{pp,2}]) = raw.result;  % Dynamic structure
                end
            catch
                fprintf('Error Loading %s\n',[mat_file_name '_' mat_address{pp,2}]);
                keyboard;
            end
        end
        
        progressbar(major_i/length(group_result));
        
    end
    fprintf('Loaded %g files (%g of them are incomplete).\n',length(group_result),not_match);
    toc
    
    %% Get some values from group_result(i).mat_raw_xxx our to group_result(i) for easier access
    %  Preprocess each .mat_raw for different purposes (which have not been done in Batch Processing)
    
    % HH20150414: Modality divergence added.
    % HH20150419: Choice Divergence/Preference and Modality Divergence/Preference have been moved to Batch processing
    % HH20150419: Different temporal alignment have been included
    
    %%%%%%% Order corresponding to "Sort_Id" in TEMPO_GUI processing %%%%
    ALL_CHOICE = 1; CHOICE_ANGLE = 2; CHOICE_DIFFICULT = 3; OUTCOME = 4;
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
        if_contralateral = xls_num{1}(i,header.Hemisphere) ~= group_result(i).PREF_PSTH;
        group_result(i).ChoicePreference = group_result(i).ChoicePreference * sign(if_contralateral - 0.5);
        
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
    
    % Load Gaussian velocity (measured by accelerometer at 109. Should retest it on 103. HH20150422)
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
    memsac_ts = group_result(end).MemSac_ts;
    
    % Initialization
    
    % j-sensitive
    for j = 1:2
        PSTH_all_Norm{j} = NaN(N,length(rate_ts{j}),6);
        PSTH_angles_Norm{j} = NaN(N,length(rate_ts{j}),1+length(group_result(end).unique_heading),3);
        PSTH_outcomes_Norm{j} = NaN(N,length(rate_ts{j}),4,3);
        
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
    MemSac_PREFmNULL_PSTH = NaN(N,length(group_result(2).MemSac_ts));
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
            
            offset = min(min(PSTH_all_Norm_this(:,norm_range)));
            gain = max(max(PSTH_all_Norm_this(:,norm_range))) - offset;
            
            PSTH_all_Norm_this = PSTH_all_Norm_this - offset;
            PSTH_all_Norm_this = PSTH_all_Norm_this / gain;
            
            for stim_type = 1:3 % Stim_type check
                
                % I have done "always output three conditions" in TEMPO_GUI for all data except CP and PSTH. @HH20150419
                k = find(stim_type == group_result(i).mat_raw_PSTH.unique_stim_type);
                if ~isempty(k)   % We have this condition
                    
                    % Pack PSTH_all_Norm
                    PSTH_all_Norm{j}(i,:,k*2-1) = PSTH_all_Norm_this(k*2-1,:);
                    PSTH_all_Norm{j}(i,:,k*2) = PSTH_all_Norm_this(k*2,:);
                    
                    % Normalize and pack PSTH_angles_Norm
                    PSTH_angles_norm_this = group_result(i).mat_raw_PSTH.PSTH{j,CHOICE_ANGLE,k}.ys;
                    PSTH_angles_norm_this = PSTH_angles_norm_this - offset;
                    PSTH_angles_norm_this = PSTH_angles_norm_this / gain;
                    
                    if size(PSTH_angles_norm_this,1) == size(PSTH_angles_Norm{j},3)
                        PSTH_angles_Norm{j}(i,:,:,k) = PSTH_angles_norm_this';
                    elseif size(PSTH_angles_norm_this,1) == size(PSTH_angles_Norm{j},3) - 2 % Without zero heading
                        PSTH_angles_Norm{j}(i,:,3:end,k) = PSTH_angles_norm_this';
                    else
                        disp('No match PSTH_angles_norm...');
                    end
                    
                    % Normalize and pack PSTH_outcome_Norm
                    PSTH_outcome_norm_this = group_result(i).mat_raw_PSTH.PSTH{j,OUTCOME,k}.ys;
                    PSTH_outcome_norm_this = PSTH_outcome_norm_this - offset;
                    PSTH_outcome_norm_this = PSTH_outcome_norm_this / gain;
                    
                    PSTH_outcomes_Norm{j}(i,:,:,k) = PSTH_outcome_norm_this';
                    
                end
                
                % Stim type check already done
                ChoiceDiv_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_ALL{j}(stim_type,:),smooth_factor_for_divergence);
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
        
        % Mem-sac stuff (Note some are NaNs)
        if ~isempty(group_result(i).MemSac_DDI)
            % Global memsac indicator
            MemSac_DDI(i,:) = group_result(i).MemSac_DDI;
            
            % Local memsac indicator (Pref/null of HD task)
            MemSac_PREFmNULL_PSTH(i,:) = (group_result(i).MemSac_PSTH(:,1) - group_result(i).MemSac_PSTH(:,5)) *  sign((group_result(i).PREF_PSTH == 2)-0.5);  % (Right - Left)* Right is pref
            MemSac_PSTH_AngDiff(i,:) = mod(group_result(i).MemSac_vectSum - ((group_result(i).PREF_PSTH == 2)*0 + (group_result(i).PREF_PSTH == 1)*180),360);
        end
        
        MemSac_PSTH_AngDiff(MemSac_PSTH_AngDiff > 180) = 360 - MemSac_PSTH_AngDiff(MemSac_PSTH_AngDiff > 180); % Ang diffs are less than 180
        
    end
    
    % Psychometric parameters
    Psy_all = [group_result.Psy_para];
    Psy_thres = reshape(Psy_all(1,:),3,[])';
    
    bayes_pred = @(x,y) sqrt(1./(1./x.^2 + 1./y.^2));
    Psy_thres_pred = bayes_pred(Psy_thres(:,1),Psy_thres(:,2));
    Psy_pred_ratio = Psy_thres(:,3)./Psy_thres_pred;
        
    %% Cell Selection and Cell Counter
    % @ HH20150413
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Limitations on Repetition number & Target first
    select_all = ([group_result.repN]' >= 8) & ([group_result.length_unique_stim_type]' == 3 ) & (xls_num{1}(:,header.HD_TargFirst)~=0);
    
    % + SUs
    select_sus = select_all & (xls_num{1}(:,header.Chan1) >= 5) & (xls_num{1}(:,header.Chan1) < 20);
    
    % + T(ypical) Cells
    select_tcells = select_sus & (xls_num{1}(:,header.HD_MemSac) == 1);
    
    % Non-T(ypical) Cells
    select_no_tcells = select_sus & (xls_num{1}(:,header.HD_MemSac) < 1);
    
    % Bottom line for most figures
    select_bottom_line = select_sus;
    
    % Psychometric good (prediction ratio <= 1.3)
    select_psy_good = Psy_pred_ratio <= 1.3;
    select_psy_bad = Psy_pred_ratio > 1.3;
    
    % select_memsac_good = xls_num{1}(:,header.HD_MemSac) == 1; % MemSac_DDI(:,4) >= 0.55;
    % select_memsac_bad = xls_num{1}(:,header.HD_MemSac) < 1; % MemSac_DDI(:,4) < 0.55;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update cell counter
    h_all = findall(gcbf,'tag','num_all_units');
    set(h_all,'string',num2str(sum(select_all)));
    h_su = findall(gcbf,'tag','num_su');
    set(h_su,'string',num2str(sum(select_sus)));
    h_tcell = findall(gcbf,'tag','num_t_cells');
    set(h_tcell,'string',num2str(sum(select_tcells)));
    
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
    
    datetick('x',26);
    
    rotateXLabels(gca,30);
    axis tight;
    SetFigure();
    
    % After running this part of code, I am literally crying...
    % Come on man!!!!!!!
    
    % Plot prediction ratio
    set(figure(9998),'name','Prediction ratio','position',[686 539 900 412]); clf;
    h_line = plot(Psy_pred_ratio,'o','markersize',10); hold on;
    plot(find(select_tcells),Psy_pred_ratio(select_tcells),'+r','linew',2,'markersize',10);
    plot(xlim,[1 1],'k--'); ylabel('Prediction ratio'); set(gca,'yscale','log','ytick',[0.6 1 2]);
    SetFigure(15);
    
    % Show individual cell selected from the figure. HH20150424
    set([gca h_line],'ButtonDownFcn',{@Show_individual_cell,h_line,true(size(Psy_pred_ratio))});
    
    % ====================================== Final Preparation ============================================
    
    set(0,'defaultAxesColorOrder',[0 0 1; 1 0 0; 0 0.8 0;]);
    colors = [0 0 1; 1 0 0; 0 0.8 0.4];
    transparent = 0.7;
    figN = 2499;

    
    find_bottom_line = find(select_bottom_line);
    
    % Overall time marker (Stim on, Stim off, Sac on)
    tt = reshape(cell2mat([group_result.time_marker]),6,[])';
    
    for j = 1:2
        time_markers{j} = [mean(tt(:,(j-1)*3+1:(j-1)*3+3)); std(tt(:,(j-1)*3+1:(j-1)*3+3))];
    end
    marker_for_time_markers{1} = {'-','-','--'};
    marker_for_time_markers{2} = {'--','--','-'};

%% ====================================== Data Prepration Done =============================================%

    %% PSTH 1. Correct only, all choices
    
    set(figure(999),'name','Average PSTH (Correct only, all choices)','pos',[27 57 919 898]); clf
    h_subplot = tight_subplot(3,2,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
    
    methods_of_select = {
        select_bottom_line, 'All cells';
        select_tcells, 'Typical cells';
        select_no_tcells, 'Non-typical cells'};
    
    for j = 1:2
        for ms = 1:3
            SeriesComparison(PSTH_all_Norm{j}(methods_of_select{ms,1},:,:),rate_ts{j},'Colors',{'b','b','r','r','g','g'},'LineStyles',{'-','--','-','--','-','--'},...
                'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(ms+(j-1)*3));
            if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            legend off;
            title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
            
            for tt = 1:3
                plot([1 1] * time_markers{j}(1,tt),[0 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
            end
            
            xlim([rate_ts{j}(10) rate_ts{j}(end-10)]); ylim([0 .8]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*max(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
        end
    end
    SetFigure(15);
    
    %% PSTH 2. Different headings
    
    
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
        h_subplot = tight_subplot(3,2,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
        
        for k = 1:3
            for j = 1
                
                yyy = PSTH_angles_Norm{j}(methods_of_select{ms,1},:,:,k);
                ttt = rate_ts{j};
                
                h = SeriesComparison(yyy,ttt,...
                    'Colors',colors_angles,'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k+(j-1)*3));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                if j == 1 && k == 1
                    title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                for tt = 1:3
                    plot([1 1] * time_markers{j}(1,tt),[0 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
                end
                
                xlim([rate_ts{j}(10) rate_ts{j}(end-10)]); ylim([0.1 .7]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
            end
            
            % Difference
            yyy_mean = squeeze(nanmean(yyy,1));
            yyy_diff =  yyy_mean(:,1:2:end) - yyy_mean(:,2:2:end);
            
            h = SeriesComparison(shiftdim(yyy_diff,-1),ttt,...
                'Colors',colors_angles(1:2:end),'LineStyles',{'-'},...
                'ErrorBar',0,'Xlabel',[],'Ylabel','Diff','axes',h_subplot(k+(2-1)*3));
            
            if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            legend off;
            
            xlim([rate_ts{j}(10) rate_ts{j}(end-10)]); ylim([-.2 .3]);
            
            for tt = 1:3
                plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
            end
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
            
            
        end
        SetFigure(15);
    end
    
    %% PSTH 3. Correct / Wrong Trials ------------
    
    methods_of_select = {
        select_bottom_line, 'All cells';
        select_tcells, 'Typical cells';
        select_no_tcells, 'Non-typical cells'};
    
    for ms = 1:3
        set(figure(999-ms),'name',['Average PSTH (Correct vs Wrong), ' methods_of_select{ms,2}],'pos',[27 57 919 898]); clf
        h_subplot = tight_subplot(3,2,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
        
        for j = 1:2
            for k = 1:3
                
                h = SeriesComparison(PSTH_outcomes_Norm{j}(methods_of_select{ms,1},:,:,k),rate_ts{j},...
                    'Colors',{'k','k','m','m'},'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k+(j-1)*3));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                if j == 1 && k == 1
                    title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                for tt = 1:3
                    plot([1 1] * time_markers{j}(1,tt),[0 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
                end
                
                xlim([rate_ts{j}(10) rate_ts{j}(end-10)]); ylim([0.1 .7]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
            end
        end
        
        SetFigure(15);
    end
    
    %% Correlations 1. Mem-sac v.s. Choice Divergence
    
    
    % mean(A(:,2:3),2)+0.4,mean(ChoiceDiv_All(selectCells,rate_ts > -500 & rate_ts < 00,1),2),'o'
    
    % Maybe I should calculate a corresponding ChoiceDiv for MemSac at PREF/NULL target locations??  @HH20150419
    
    j = 2;  % Not much time-sensitive, so I use j = 2 here.
    
    nHist = 15;
    
    LinearCorrelation(mean(MemSac_DDI(select_bottom_line,4),2),...
        [  mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),...
        mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),...
        mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,3),2)],...
        'Xlabel','MemSac DDI (pre)','Ylabel','Pre-sac CDiv','FaceColors',{'b','r','g'},'Markers',{'o'},...
        'LineStyles',{'b-','r-','g-'},'MarkerSize',12,'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
    set(gcf,'name',['j = ' num2str(j)]);
    
    % Annotate tcells
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),...
        '+','markersize',16,'color','k','linew',2);
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),...
        '+','markersize',16,'color','k','linew',2);
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,3),2),...
        '+','markersize',16,'color','k','linew',2);
    
    LinearCorrelation(mean(MemSac_DDI(select_bottom_line,4),2),...
        [  mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > 100 & rate_ts{j} < 500,1),2),...
        mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > 100 & rate_ts{j} < 500,2),2),...
        mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > 100 & rate_ts{j} < 500,3),2)],...
        'Xlabel','MemSac DDI (pre)','Ylabel','Post-sac CDiv','FaceColors',{'b','r','g'},'Markers',{'o'},...
        'LineStyles',{'b-','r-','g-'},'MarkerSize',12,'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
    set(gcf,'name',['j = ' num2str(j)]);
    % Annotate tcells
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > 100 & rate_ts{j} < 500,1),2),...
        '+','markersize',16,'color','k','linew',2);
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > 100 & rate_ts{j} < 500,2),2),...
        '+','markersize',16,'color','k','linew',2);
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > 100 & rate_ts{j} < 500,3),2),...
        '+','markersize',16,'color','k','linew',2);
    
    
    % LinearCorrelation(mean(MemSac_PREFmNULL_PSTH(select_bottom_line,memsac_ts > -300 & memsac_ts < -100),2),...
    %     [  mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),...
    %        mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),...
    %        mean(ChoiceDiv_All{j}(select_bottom_line,rate_ts{j} > -400 & rate_ts{j} < -100,3),2)],...
    %         'Xlabel','MemSac Rate Diff (pre)','Ylabel','Pre-sac Cdiv','FaceColors',{'b','r','g'},'Markers',{'o'},...
    %         'LineStyles',{'b-','r-','g-'},'MarkerSize',12,'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
    
    %% Correlations 2. Mem-sac v.s. CP
    
    j = 2; % Not much time-sensitive, so I use j = 2 here.
    CP_interest = -500;
    
    CP_interval = CP_ts{j} == CP_interest ;
    
    for k = 1:3
        CP_NS{k} = find_bottom_line(CP_p{j}(select_bottom_line,CP_interval,k) > 0.05);
        CP_S{k} = find_bottom_line(CP_p{j}(select_bottom_line,CP_interval,k) <= 0.05);
    end
    
    
    h = LinearCorrelation(...
        {   mean(MemSac_DDI(CP_NS{1},4),2),...
        mean(MemSac_DDI(CP_S{1},4),2),...
        mean(MemSac_DDI(CP_NS{2},4),2),...
        mean(MemSac_DDI(CP_S{2},4),2),...
        mean(MemSac_DDI(CP_NS{3},4),2),...
        mean(MemSac_DDI(CP_S{3},4),2),...
        },...
        {  mean(CP{j}(CP_NS{1},CP_interval,1),2),...
        mean(CP{j}(CP_S{1},CP_interval,1),2),...
        mean(CP{j}(CP_NS{2},CP_interval,2),2),...
        mean(CP{j}(CP_S{2},CP_interval,2),2),...
        mean(CP{j}(CP_NS{3},CP_interval,3),2),...
        mean(CP{j}(CP_S{3},CP_interval,3),2)},...
        'CombinedIndex',[3 12 48],...
        'Xlabel','MemSac DDI (pre)','Ylabel',['CP at Sac ' num2str(CP_interest)],...
        'FaceColors',{'none','b','none','r','none','g'},'Markers',{'o'},...
        'LineStyles',{'b:','b:','r:','r:','g:','g:','b-','r-','g-'},'MarkerSize',12,...
        'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
    set(gcf,'name',['j = ' num2str(j)]);
    delete([h.group(1:6).line]);
    plot(xlim,[0.5 0.5],'k--');
    
    % Annotate tcells
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(CP{j}(select_tcells,CP_interval,1),2),...
        '+','markersize',16,'color','k','linew',2);
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(CP{j}(select_tcells,CP_interval,2),2),...
        '+','markersize',16,'color','k','linew',2);
    plot(mean(MemSac_DDI(select_tcells,4),2),mean(CP{j}(select_tcells,CP_interval,3),2),...
        '+','markersize',16,'color','k','linew',2);
    
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
        'LineStyles',{'b:','b:','r:','r:','g:','g:','b-','r-','g-'},'MarkerSize',12,...
        'figN',figN,'XHist',nHist,'YHist',nHist,'SameScale',1); figN = figN + 1;
    delete([h.group(1:6).line]);
    plot(xlim,[0.5 0.5],'k--'); plot([0.5 0.5],ylim,'k--');
    
    
    set(gcf,'name',['j = ' num2str(j)]);
    % Annotate tcells
    plot(.5+mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,1),2),mean(CP{j}(select_tcells,CP_interval,1),2),...
        '+','markersize',16,'color','k','linew',2);
    plot(.5+mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,2),2),mean(CP{j}(select_tcells,CP_interval,2),2),...
        '+','markersize',16,'color','k','linew',2);
    plot(.5+mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > -400 & rate_ts{j} < -100,3),2),mean(CP{j}(select_tcells,CP_interval,3),2),...
        '+','markersize',16,'color','k','linew',2);
   %} 
   
    %% Correlations 3. Choice preference v.s. modality preference (Anne Fig.2f-h; Fig.3a)
    
    select_cpref_mpref = select_bottom_line;
    
    Choice_pref_all = reshape([group_result(select_cpref_mpref).ChoicePreference]',3,[],2);  % Stim, Cell No, Pre/Post
    Choice_pref_p_value_all = reshape([group_result(select_cpref_mpref).ChoicePreference_pvalue]',3,[],2);
    Modality_pref_all = reshape([group_result(select_cpref_mpref).ModalityPreference]',3,[],2);
    Modality_pref_p_value_all = reshape([group_result(select_cpref_mpref).ModalityPreference_pvalue]',3,[],2);
    
    % ---------------- 3.1  Choice pref v.s. modality pref (Anne Fig.3a)-------------------
    tt = 1;
    
    for k = 1:3
        set(figure(figN),'name','CDiv vs. MDiv','pos',[17 514 1151 449]); 
        cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
        mpref_sig = Modality_pref_p_value_all(1,:,tt) < 0.05;
        
        h = LinearCorrelation({
            (Modality_pref_all(1, ~cpref_sig & ~mpref_sig,tt)) ;
            (Modality_pref_all(1, cpref_sig | mpref_sig,tt));
            (Modality_pref_all(1, cpref_sig & mpref_sig,tt))},...
            {
            (Choice_pref_all(k,~cpref_sig & ~mpref_sig,tt)) ;
            (Choice_pref_all(k,cpref_sig | mpref_sig,tt)) ;
            (Choice_pref_all(k,cpref_sig & mpref_sig,tt)) },...
            'CombinedIndex',[7],...
            'Ylabel','Choice preference (pre)','Xlabel','Modality preference (pre)',...
            'FaceColors',{'none',colors(k,:)*0.2 + [0.8 0.8 0.8],colors(k,:)},'Markers',{'o'},...
            'LineStyles',{'k:','k:','k:','k-'},'MarkerSize',12,...
            'figN',figN,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,'Method','Spearman'); figN = figN + 1;
        delete([h.group(1:3).line h.diag]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
        
        % Annotate tcells
        plot(Modality_pref_all(1,select_tcells(select_bottom_line),tt),Choice_pref_all(k,select_tcells(select_bottom_line),tt),'+','markersize',16,'color','k','linew',2);
        
        %     for i = 1:3
        %         set(h.group(i).dots,'color',colors(k,:));
        %     end
    end
    
    % ---------------- 3.2 Choice pref (vest and visual) -------------------
    
    
    tt = 1;
    
    set(figure(figN),'name','CDiv (visual) vs. CDiv (vest)','pos',[17 514 1151 449]);
    
    cpref_sig_1 = Choice_pref_p_value_all(1,:,tt) < 0.05;
    cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
    
    h = LinearCorrelation({
        (Choice_pref_all(2,~cpref_sig_1 & ~cpref_sig_2,tt)) ;
        (Choice_pref_all(2,cpref_sig_1 | cpref_sig_2,tt)) ;
        (Choice_pref_all(2,cpref_sig_1 & cpref_sig_2,tt))
        },...
        {
        (Choice_pref_all(1,~cpref_sig_1 & ~cpref_sig_2,tt)) ;
        (Choice_pref_all(1,cpref_sig_1 | cpref_sig_2,tt)) ;
        (Choice_pref_all(1,cpref_sig_1 & cpref_sig_2,tt))
        },...
        'CombinedIndex',[7],...
        'Ylabel','Vestibular choice preference','Xlabel','Visual choice preference',...
        'FaceColors',{'none',[0.8 0.8 0.8],'k'},'Markers',{'o'},...
        'LineStyles',{'k:','k:','k:','k-'},'MarkerSize',12,...
        'figN',figN,'XHist',20,'YHist',20,...
        'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,'Method','Spearman'); figN = figN + 1;
    delete([h.group(1:3).line]);
    plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
    
    % Annotate tcells
    plot(Choice_pref_all(2,select_tcells(select_bottom_line),tt),Choice_pref_all(1,select_tcells(select_bottom_line),tt),...
        '+','markersize',16,'color','m','linew',2);
    
    %     for i = 1:3
    %         set(h.group(i).dots,'color',colors(k,:));
    %     end
    
    % ---------------- 3.3 Choice pref (single vs comb) -------------------
    
    tt = 1;
    
    set(figure(figN),'name','CDiv (single) vs. CDiv (comb)','pos',[17 514 1151 449]); 
    
    cpref_sig_1 = Choice_pref_p_value_all(1,:,tt) < 0.05;
    cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
    cpref_sig_3 = Choice_pref_p_value_all(3,:,tt) < 0.05;
    
    h = LinearCorrelation({
        (Choice_pref_all(1,~cpref_sig_1 & ~cpref_sig_3,tt)) ;
        (Choice_pref_all(1,cpref_sig_1 | cpref_sig_3,tt)) ;
        (Choice_pref_all(1,cpref_sig_1 & cpref_sig_3,tt)) ;
        (Choice_pref_all(2,~cpref_sig_2 & ~cpref_sig_3,tt)) ;
        (Choice_pref_all(2,cpref_sig_2 | cpref_sig_3,tt)) ;
        (Choice_pref_all(2,cpref_sig_2 & cpref_sig_3,tt))
        
        },...
        {
        (Choice_pref_all(3,~cpref_sig_1 & ~cpref_sig_3,tt)) ;
        (Choice_pref_all(3,cpref_sig_1 | cpref_sig_3,tt)) ;
        (Choice_pref_all(3,cpref_sig_1 & cpref_sig_3,tt)) ;
        (Choice_pref_all(3,~cpref_sig_2 & ~cpref_sig_3,tt)) ;
        (Choice_pref_all(3,cpref_sig_2 | cpref_sig_3,tt)) ;
        (Choice_pref_all(3,cpref_sig_2 & cpref_sig_3,tt))
        
        },...
        'CombinedIndex',[7 56],...
        'Xlabel','Single modality choice preference','Ylabel','Combined choice preference',...
        'FaceColors',{'none',[0.8 0.8 1],'b','none',[1 0.8 0.8],'r'},'Markers',{'o'},...
        'LineStyles',{'k:','k:','k:','k:','k:','k:','b-','r-'},'MarkerSize',12,...
        'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
        'SameScale',1,'Method','Spearman'); figN = figN + 1;
    
    delete([h.group(1:6).line]);
    plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
    
    % Annotate tcells
    % plot(Choice_pref_all(2,select_tcells(select_bottom_line),tt),Choice_pref_all(1,select_tcells(select_bottom_line),tt),...
    %     '+','markersize',16,'color','m','linew',2);
    
    %     for i = 1:3
    %         set(h.group(i).dots,'color',colors(k,:));
    %     end
    
    % ---------------- 3.4 Choice pref (pre vs post) -------------------
    
    set(figure(figN),'name','CDiv (pre-sac) vs. CDiv (post-sac)','pos',[17 514 1151 449]); 
    
    k = 3;
    
    cpref_sig_pre = Choice_pref_p_value_all(k,:,1) < 0.05;
    cpref_sig_post = Choice_pref_p_value_all(k,:,2) < 0.05;
    
    h = LinearCorrelation({
        (Choice_pref_all(k,~cpref_sig_pre & ~cpref_sig_post,1)) ;
        (Choice_pref_all(k,cpref_sig_pre | cpref_sig_post,1)) ;
        (Choice_pref_all(k,cpref_sig_pre & cpref_sig_post,1)) ;
        },...
        {
        (Choice_pref_all(k,~cpref_sig_pre & ~cpref_sig_post,2)) ;
        (Choice_pref_all(k,cpref_sig_pre | cpref_sig_post,2)) ;
        (Choice_pref_all(k,cpref_sig_pre & cpref_sig_post,2)) ;
        
        },...
        'CombinedIndex',[7],...
        'Xlabel','Choice preference: decision formation','Ylabel','Choice preference: Movement',...
        'FaceColors',{'none',colors(k,:)*0.2 + [0.8 0.8 0.8],colors(k,:)},'Markers',{'o'},...
        'LineStyles',{'k:','k:','k:','k-'},'MarkerSize',12,...
        'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
        'SameScale',1,'Method','Spearman'); figN = figN + 1;
    
    delete([h.group(1:3).line]);
    plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
    
    % Annotate tcells
    plot(Choice_pref_all(k,select_tcells(select_bottom_line),1),Choice_pref_all(k,select_tcells(select_bottom_line),2),...
        '+','markersize',16,'color','k','linew',2);
    
    %     for i = 1:3
    %         set(h.group(i).dots,'color',colors(k,:));
    %     end
    
 %%{   
 
    %% Correlations 4. Psychophysics v.s. CD
    
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
        'LineStyles',{'b-','r-','k-'},'MarkerSize',12,...
        'Ylabel',sprintf('\\Delta CDiv (%g - %g ms, 3-1, 3-2, 1-2)',t_begin,t_end),'Xlabel','Psycho prediction ratio',...
        'MarkerSize',12,...
        'figN',figN,'XHist',15,'YHist',15,'logx',1,...
        'XHistStyle','stacked','YHistStyle','stacked','Method','Pearson'); figN = figN + 1;
    plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
    
    %% State Space 1. Hot-gram
    
    % Pack all data into matrix A: [memsacDDI, vest_div, vis_div, comb_div]
    
    % Now it becomes A: [memsacDDI, vest_div, vis_div, comb_div, 2-1 mod_div, 3-1 mod_div, 3-2 mod_div]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    j = 1;
    enlarge_factor = 30; % Enlarge memsac DDIs
    
    sort_time_interval1 = [find(rate_ts{j} > time_markers{j}(1,1), 1)  find(rate_ts{j} >= time_markers{j}(1,3),1)-1]; % Stim on to saccade on
    sort_time_interval2 = [find(rate_ts{j} > time_markers{j}(1,3),1)  length(rate_ts{j})]; % Post sac
    sort_time_interval3 = [1  length(CP_ts{j})];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Exclude any cells that have NaNs
    % for k = 1:3
    %     selectCells = selectCells & ~isnan(ChoiceDiv_All(:,1,k));
    % end
    
    
    % --- Sort A according to different methods for visualizing ---
    
    A_CP = reshape(CP{j}(select_bottom_line,:,:),sum(select_bottom_line),[])-0.5;
    A_choicediv= reshape(ChoiceDiv_All{j}(select_bottom_line,:,:),sum(select_bottom_line),[]);
    
    % A_memSac = xls_num{1}(selectCells,header.DDI_LS:header.DDI_post) - 0.4; % All memSac data
    % LS, M, Pre, Co, Post
    A_memSac = MemSac_DDI(select_bottom_line,2:end) - 0.4;  % Read memsac from .mat file. HH20150413
    
    
    % Now I add mod_div. HH20150415
    A_moddiv = reshape(ModDiv_All{j}(select_bottom_line,:,:),sum(select_bottom_line),[]);
    
    A = [A_memSac A_choicediv A_moddiv]; % A_CP];
    
    sort_method = {
        [1 4];
        %     [3 4];
        
        %     5 + 0*length(rate_ts{j}) + sort_time_interval1;
        %     5 + 1*length(rate_ts{j}) + sort_time_interval1;
        5 + 2*length(rate_ts{j}) + sort_time_interval1;
        
        %     5 + 0*length(rate_ts{j}) + sort_time_interval2;
        %     5 + 1*length(rate_ts{j}) + sort_time_interval2;
        5 + 2*length(rate_ts{j}) + sort_time_interval2;
        
        5 + 3* length(rate_ts{j}) + sort_time_interval1;
        %     5 + 4* length(rate_ts{j}) + sort_time_interval1;
        %     5 + 5* length(rate_ts{j}) + sort_time_interval1;
        
        
        %     5 + 6* length(rate_ts{j}) + 0*length(CP_ts{j}) + sort_time_interval3;
        %     5 + 6* length(rate_ts{j}) + 1*length(CP_ts{j}) + sort_time_interval3;
        %     5 + 6* length(rate_ts{j}) + 2*length(CP_ts{j}) + sort_time_interval3;
        };
    
    for sort_ind = 1:length(sort_method)
        
        sort_begin = sort_method{sort_ind}(1);
        sort_end = sort_method{sort_ind}(end);
        [~,sort_order] = sort(mean(A(:,sort_begin:sort_end),2));
        
        % Enlarge mem-sac part for clarity
        A_forplot = [reshape(repmat(A_memSac,enlarge_factor,1),size(A_memSac,1),[]) A_choicediv A_moddiv A_moddiv(:,end)];% A_CP A_CP(:,end)]; % To ensure the last value to be shown
        A_forplot = A_forplot(sort_order,:);
        
        A_forplot = [A_forplot; A_forplot(end,:)]; % To ensure the last cell to be shown
        
        
        [nn, tt] = meshgrid(1:size(A_forplot,1), 1:size(A_forplot,2));
        
        set(figure(2099+figN),'name',['CDiv and MDiv, hot-gram, j = ' num2str(j)],'pos',[409 182 913 775]); clf; figN = figN+1;
        surf(tt,nn,A_forplot','Edgecolor','none');
        view(0,90);
        
        % Annotate time separations
        hold on;
        CD_begin_at = enlarge_factor * 5 + 1;
        
        
        % Stim on / stim off / sac on
        for ttt = 1:3
            for tt = 0:5
                plot3(CD_begin_at + tt*length(rate_ts{j}) + find(rate_ts{j} >= time_markers{j}(1,ttt),1) * [1 1],...
                    [1 size(A_forplot,1)],[1 1],'k','linesty',marker_for_time_markers{j}{ttt},'linewid',1.5);
            end
        end
        
        % End
        plot3([CD_begin_at CD_begin_at],[1 size(A_forplot,1)],[1 1],'k','linewid',3);
        for tt = 1:6
            plot3(CD_begin_at + tt*[length(rate_ts{j})  length(rate_ts{j})],[1 size(A_forplot,1)],[1 1],'k','linewid',3);
        end
        
        % Annotate tcells
        cell_loc = select_tcells(select_bottom_line);
        plot3(0,find(cell_loc(sort_order))+.5,1,'+','markersize',5,'color','k','linew',1.5);
        
        % Annotate sort methods
        temp_beg = @(x)((x<=5)*((x-1)*enlarge_factor) + 1 + (x>5)*(CD_begin_at+(x-5)-1));
        temp_end = @(x)((x<=5)*((x)*enlarge_factor) + 1 + (x>5)*(CD_begin_at+(x-5)));
        
        sort_begin_forplot = temp_beg(sort_begin);
        sort_end_forplot = temp_end(sort_end);
        
        plot3([sort_begin_forplot sort_end_forplot],[size(A_forplot,1)+1 size(A_forplot,1)+1],[1 1],'r','linewid',5);
        xlabel('Temporal features');
        ylabel('Cell Number');
        
        set(gca,{'xtick','ytick','ztick'},{[],10:10:size(A_forplot,1),[]},'linewidth',0.00000001);
        axis tight; colorbar; ylim([1 size(A,1)+2.5]); SetFigure();
        
        
    end
    
    %% State Space 2. PCA along feature dimensions (Eigenfeature)
    
    % Use the original A instead of A_forplot
    
    dataPCA = [A_memSac A_choicediv];% A_moddiv]; % A_CP ];
    dataPCA(sum(isnan(dataPCA),2)>0,:) = [];
    
    % [~,sort_order] = sort(mean(dataPCA(:,3:4),2)); % Sort according to mem-sac
    [~,sort_order] = sort(mean(dataPCA(:, 5 + 0*length(rate_ts{j}) + sort_time_interval1),2));
    dataPCA = dataPCA(sort_order,:);
    
    % Do PCA
    [PC, ~] = pca(dataPCA);
    
    % The first three eigenvectors that have the first three largest
    % eigenvalues.
    PC1 = PC(:,1);
    PC2 = PC(:,2);
    % PC3 = sortedEigVectors(:,3);
    
    % Projecting the raw data onto the first three eigenvectors
    projPC1 = dataPCA * PC1;
    projPC2 = dataPCA * PC2;
    % projPC3 = dataPCA(:,2:end) * PC3;
    
    set(figure(2099+figN),'name',['Scatter of Eigenfeature, j = ' num2str(j)],'pos',[41 445 543 451]); clf; figN = figN+1;
    
    colorOrder = jet(length(projPC1));
    
    scatter(projPC1, projPC2, 150, colorOrder,'fill');
    grid on; hold on;
    xlabel('PC1'); ylabel('PC2');
    SetFigure();
    
    set(figure(2099+figN),'pos',[602 446 978 450],'name',['Eigenfeature, j = ' num2str(j)]); clf; figN = figN+1;
    
    for ss = 1:2
        subplot(2,1,ss);
        
        eval(['plot(PC' num2str(ss) ',''lineW'',2)']);
        
        hold on;
        
        plot([6 6],ylim,'k','linewid',1);
        plot(xlim,[0 0],'k--');
        
        % Stim on / stim off / sac on
        for ttt = 3
            for tt = 0:5
                plot(6 + tt*length(rate_ts{j}) + find(rate_ts{j} >= time_markers{j}(1,ttt),1) * [1 1],...
                    ylim,'k','linesty',marker_for_time_markers{j}{ttt},'linewid',1.5);
            end
        end
        
        % End
        plot([6 6],ylim,'k','linewid',3);
        for tt = 1:6
            plot(6 + tt*[length(rate_ts{j})  length(rate_ts{j})],ylim,'k','linewid',3);
        end
        
        axis tight; axis off
    end
    
    
    SetFigure();
  %}  
    
    %% State Space 3. PCA along neuron dimensions (EigenNeuron)  % HH20150413
    
    % Pack all data (PSTH) into matrix B: [vest_pref, vest_null, vis_pref, vis_null, comb_pref, comb_null]
    % Only decicion period is included
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = 1;
    select_for_PCA_B = select_bottom_line ;
    PCA_B_time_range = min(rate_ts{j})+100 <= rate_ts{j} & rate_ts{j} <= time_markers{j}(1,3);  % Before saccade
    PCA_B_times = rate_ts{j}(PCA_B_time_range);
    denoised_dim = 12;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    find_PCA_B = find(select_for_PCA_B);
    
    B = nan(sum(select_for_PCA_B),6 * sum(PCA_B_time_range));
    
    for i = 1:sum(select_for_PCA_B)
        raw_PSTH = group_result(find_PCA_B(i)).mat_raw_PSTH.PSTH{j,ALL_CHOICE,1}.ys;
        if size(raw_PSTH,1)==6
            B(i,:) = reshape(raw_PSTH(:,PCA_B_time_range)',[],1)';
        end
    end
    
    B(sum(isnan(B),2) > 0,:) = [];
    
    % Do PCA
    [PC, score, latent, ~, explained] = pca(B');
    
    % Manually rotate PC1 and PC2 (they're also orthogonal, so there's no
    % change in PCA results). @HH20150424
    theta = 55 / 180 * pi; % Rotate counterclockwise the two PC bases
    THETA = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
    PC(:,1:2) = PC(:,1:2) * THETA;
    PC(:,2) = -PC(:,2);
    
    % Projecting the raw data onto the first three eigenvectors
    for dim = 1:denoised_dim
        projPC{dim} = (reshape(PC(:,dim)' * B,[],6))';
        %     projPC{dim} = reshape(score(:,dim),[],6)';
    end
    
        %% PCA 1. Variance explained
    set(figure(2099+figN),'pos',[936 491 733 463],'name',['Variance explained, j = ' num2str(j)]); clf; figN = figN+1; hold on;
    plot((1:length(explained))', cumsum(explained),'o-','markersize',8,'linew',1.5);
    plot((1:denoised_dim)',cumsum(explained(1:denoised_dim)),'ro-','markersize',8,'linew',1.5,'markerfacecol','r');
    plot([0 1],[0 explained(1)],'r-','linew',1.5);
    plot(xlim,[1 1]*sum(explained(1:denoised_dim)),'r--');
    text(denoised_dim,sum(explained(1:denoised_dim))*0.9,[num2str(sum(explained(1:denoised_dim))) '%'],'color','r');
    SetFigure(); xlabel('Num of principal components'); ylabel('Explained variability (%)'); ylim([0 100]);
    
        %% PCA 2. PCA weights for naive 2-D space
    set(figure(2099+figN),'pos',[13 53 836 702],'name',['Variance explained, j = ' num2str(j)]); clf; figN = figN+1; hold on;

    tcell_in_bottom_line = select_tcells(select_for_PCA_B);
    plot3(PC(tcell_in_bottom_line,1),PC(tcell_in_bottom_line,2),PC(tcell_in_bottom_line,3),'r+','markersize',10,'linew',1.5);

    h_line = plot3(PC(:,1),PC(:,2),PC(:,3),'o','markersize',10,'linew',1.5);
    SetFigure(); xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); grid on; axis tight;
    
    % Show individual cell selected from the figure. HH20150424
    set([gca h_line],'ButtonDownFcn',{@Show_individual_cell,h_line,select_for_PCA_B});
    
    % ======= Correlation between PCA weight 1 and Choice preference =========
    tt = 1;
    k = 1;
    
    set(figure(figN),'name','PCA weight 1 vs. Choice div','pos',[17 514 1151 449]);
    
    cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
    
    h = LinearCorrelation({
        (PC(~cpref_sig,1));
        (PC(cpref_sig,1));
        },...
        {
        abs(Choice_pref_all(k,~cpref_sig ,tt)) ;
        abs(Choice_pref_all(k,cpref_sig,tt)) ;
        },...
        'CombinedIndex',[3],...
        'Ylabel','Abs(choice preference)','Xlabel','Contribution of Eigen-neuron 1',...
        'FaceColors',{'none',(colors(k,:)-1)*0.5 + 1},'Markers',{'o'},...
        'LineStyles',{'k-'},'MarkerSize',12,...
        'figN',figN,'XHist',20,'YHist',20,...
        'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
    delete([h.group(1:2).line]);
    plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
    
    % Annotate tcells
    plot(PC(select_tcells(select_for_PCA_B),1),abs(Choice_pref_all(k,select_tcells(select_for_PCA_B),tt)),'+','markersize',16,'color','k','linew',2);

    % ======= Correlation between Memsac DDI and PCA weight 1  =========
    
        MemSac_DDI_phase = [2:4];
    
    h = LinearCorrelation({
        mean(MemSac_DDI(select_for_PCA_B,MemSac_DDI_phase),2); 
        },...
        {
        (PC(:,1));
        },...
        'CombinedIndex',[],...
        'Xlabel',['Memsac DDI (' num2str(MemSac_DDI_phase) ')'],'Ylabel','Contribution of Eigen-neuron 1',...
        'FaceColors',{'none',(colors(k,:)-1)*0.5 + 1},'Markers',{'o'},'Markers',{'o'},...
        'LineStyles',{'k-'},'MarkerSize',12,...
        'figN',figN,'XHist',20,'YHist',20,...
        'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
    
    plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
    
    % Annotate tcells
    plot(mean(MemSac_DDI(select_tcells,MemSac_DDI_phase),2),PC(select_tcells(select_for_PCA_B),1),...
        '+','markersize',16,'color','k','linew',2);
    
    % ======= Correlation between Memsac DDI and Choice Preference  =========
        
    k = 1;
    tt = 1;
    MemSac_DDI_phase = [2:4];
    cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.05;
    
    MemSac_DDI_selected = MemSac_DDI(select_bottom_line,:);
    h = LinearCorrelation({
        mean(MemSac_DDI_selected((~cpref_sig),MemSac_DDI_phase),2); 
        mean(MemSac_DDI_selected((cpref_sig),MemSac_DDI_phase),2); 
        },...
        {
        abs(Choice_pref_all(k,~cpref_sig,tt)) ;
        abs(Choice_pref_all(k,cpref_sig,tt)) ;
        },...
        'CombinedIndex',[3],...
        'Xlabel',['Memsac DDI (' num2str(MemSac_DDI_phase) ')'],'Ylabel','Abs(choice preference)',...
        'FaceColors',{'none',(colors(k,:)-1)*0.5 + 1},'Markers',{'o'},...
        'LineStyles',{'k-'},'MarkerSize',12,...
        'figN',figN,'XHist',20,'YHist',20,...
        'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
    
    plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
    delete([h.group(1:2).line]);
    % Annotate tcells
    plot(mean(MemSac_DDI_selected(select_tcells(select_for_PCA_B),MemSac_DDI_phase),2),abs(Choice_pref_all(k,select_tcells(select_for_PCA_B),tt)),...
        '+','markersize',16,'color','k','linew',2);
    


    % ======= Correlation between PCA weight 2 and Modality preference =========
    tt = 1;
    k = 1; % Vis-vest
    
    set(figure(figN),'name','PCA weight 1 vs. Choice div','pos',[17 514 1151 449]);
    
    mpref_sig = Modality_pref_p_value_all(1,:,tt) < 0.05;
    
    h = LinearCorrelation({
        (PC(~mpref_sig,2));
        (PC(mpref_sig,2));
        },...
        {
        (Modality_pref_all(k,~mpref_sig,tt)) ;
        (Modality_pref_all(k,mpref_sig,tt)) ;
        },...
        'CombinedIndex',[3],...
        'Ylabel','Modality preference (vis - vest)','Xlabel','Contribution of Eigen-neuron 2',...
        'FaceColors',{'none',(colors(k,:)-1)*0.5 + 1},'Markers',{'o'},...
        'LineStyles',{'k-'},'MarkerSize',12,...
        'figN',figN,'XHist',20,'YHist',20,...
        'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman'); figN = figN + 1;
    delete([h.group(1:2).line]);
    plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
    
    % Annotate tcells
    plot(PC(select_tcells(select_bottom_line),2),(Modality_pref_all(k,select_tcells(select_bottom_line),tt)),'+','markersize',16,'color','k','linew',2);
    
        %% PCA 3. Show eigenneurons' trajectory
        
    % ======== 1D ========= %
    set(figure(2099+figN),'pos',[462 67 1207 889],'name',['Population Dynamics, j = ' num2str(j)]); clf; figN = figN+1; hold on;
    ds = [1:9];
    
    for d = 1:length(ds)
        ranges(d) = range(projPC{ds(d)}(:));
    end
    
    % Plotting
    for d = 1:length(ds)
        % Normalize to [-1,1]
        gain = max(ranges) / (2 * 0.9) ;
        offset = mean(projPC{ds(d)}(:));
        norm_proj_PC_this = (projPC{ds(d)}-offset)/gain;

        h = subplot(fix(sqrt(length(ds))),ceil(length(ds)/fix(sqrt(length(ds)))),d);
        SeriesComparison(shiftdim(norm_proj_PC_this',-1),PCA_B_times,'Colors',{'b','b','r','r',[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'axes',h); 
        axis tight; ylim([-1 1]); 
        for tt = 1:3
            plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
        end
        set(gca,'ytick',[-1 0 1]);
        xlabel([]); ylabel('Amplitude (a.u.)'); title(['Eigen-neuron ' num2str(ds(d))]); legend off;
    end
    % set(get(gcf,'children'),'ylim',[min_ylim max_ylim]);
    SetFigure(15);
        
    % ======== 2D ========= %
    set(figure(2099+figN),'pos',[18 170 898 786],'name',['Population Dynamics, j = ' num2str(j)]); clf; figN = figN+1; hold on;
    
    which_two_dimension = [1,2];
    
    for k = 1:3
        %     % Pref
        %     plot3(projPC1((k-1)*2+1,:),projPC2((k-1)*2+1,:),projPC3((k-1)*2+1,:),'-','color',colors(k,:),'linew',3);
        %     % Null
        %     plot3(projPC1(k*2,:),projPC2(k*2,:),projPC3(k*2,:),':','color',colors(k,:),'linew',3);
        %
        %     % Time markers
        %     tt = find(rate_ts{j} > -500 , 1);
        %     % Pref
        %     plot3(projPC1((k-1)*2+1,tt),projPC2((k-1)*2+1,tt),projPC3((k-1)*2+1,tt),'o','color',colors(k,:),'markersize',20,'markerfacecol',colors(k,:));
        %     % Null
        %     plot3(projPC1(k*2,tt),projPC2(k*2,tt),projPC3(k*2,tt),'o','color',colors(k,:),'markersize',20,'linew',2);
        
        % Time markers
        start_time = 1; % Start point
        % Pref
        h_pref(k) = plot(projPC{which_two_dimension(1)}((k-1)*2+1,start_time),projPC{which_two_dimension(2)}((k-1)*2+1,start_time),'o','color',colors(k,:),'markersize',20,'markerfacecol',colors(k,:));
        % Null
        h_null(k) = plot(projPC{which_two_dimension(1)}(k*2,start_time),projPC{which_two_dimension(2)}(k*2,start_time),'o','color',colors(k,:),'markersize',20,'linew',3);
        
        % Pref
        plot(projPC{which_two_dimension(1)}((k-1)*2+1,:),projPC{which_two_dimension(2)}((k-1)*2+1,:),'-','color',colors(k,:),'linew',3);
        % Null
        plot(projPC{which_two_dimension(1)}(k*2,:),projPC{which_two_dimension(2)}(k*2,:),'--','color',colors(k,:),'linew',3);
        
        
    end
    
    axis tight;  grid off;
    axis off
    % xlabel('PC1'); ylabel('PC2');
    
    % xlims = xlim; ylims = ylim;
    % h_text = text(xlims(1),ylims(2),'');
    
    
    % ======= Euler distance =========
    h_timeaxis = axes('pos', [0.026 0.716 0.344 0.263] ,'color','none');  hold on

    prefs = [1 3 5];
    nulls = [2 4 6];
    
    distance = sqrt((projPC{1}(prefs,:) - projPC{1}(nulls,:)).^2 + (projPC{2}(prefs,:) - projPC{2}(nulls,:)).^2 + (projPC{3}(prefs,:) - projPC{3}(nulls,:)).^2);
    
    % Normalize distance
    distance = distance';
    % distance = (distance - repmat(min(distance,[],1),size(distance,1),1))./repmat(range(distance),size(distance,1),1);
    % set(figure(2099+figN),'pos',[18 355 766 601],'name',['Euler distance, j = ' num2str(j)]); clf; figN = figN+1; hold on;
    
    for k = 1:3
        plot(repmat(PCA_B_times',1,3),distance(:,k),'color',colors(k,:),'linew',3);
    end
    
    % Gaussian vel
    plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*max(ylim)/4,'--','linew',3,'color',[0.6 0.6 0.6]);
    axis tight
    
    for tt = 1:3
        plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
    end
    
    text(mean(xlim)/3*2,max(ylim)*.8,['N = ' num2str(sum(select_for_PCA_B))]);
    
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
    uicontrol('Style','slider','Unit','norm','Position',[0.026 0.63 0.344 0.03],...
        'Min',1,'Max',length(PCA_B_times),'Value',1,...
        'Callback',{@Population_dynamics,999});
 
    %%{
    
    %% CP: Targ first, all units
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    select_for_CP = select_bottom_line;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    set(figure(2099+figN),'name','CP','pos',[287 483 1161 475]); figN = figN+1; clf;
    
    for j = 1:2
        subplot(1,2,j);
        
        for k = 1:3
            
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
    
    %% CDiv: Time-course (Modality divergence)
    
    % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) >1) & (xls_num{1}(:,header.HD_TargFirst)~=0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    select_for_div = select_bottom_line;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        for k = 1:3
            selectCells_notNaN =  select_for_div & (~isnan(ChoiceDiv_All{j}(:,1,k)));
            
            ys_CD(k,:) = mean(ChoiceDiv_All{j}(selectCells_notNaN,:,k));
            errors = std(ChoiceDiv_All{j}(selectCells_notNaN,:,k))/sqrt(sum(selectCells_notNaN));
            h = shadedErrorBar(rate_ts{j},ys_CD(k,:),errors,{'Color',colors(k,:),'linestyle','-'},transparent);
            set(h.mainLine,'LineWidth',2);
        end
        
        axis tight;  ylim([-0.1 0.25]);
        plot(xlim,[0 0],'k--');
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',3,'color',[0.6 0.6 0.6]);
        
        title(sprintf('Choice divergence (SU, N = %g)',sum(select_for_div)));
        
        for tt = 1:3
            plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
        end
        
    end
    
    SetFigure(17);
    
    %% CDiv: Multisensory Enhancement
    
    % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) > 0) & (xls_num{1}(:,header.HD_TargFirst)~=0);% ...
    % & (~isnan(ChoiceDiv_All(:,1,k))) ;% & MemSac_DDI(:,4)<0.55 ;
    
    set(figure(2099+figN),'name','Multisensory Enhancement of CDiv','pos',[12 276 1272 684]); clf; figN = figN+1;
    
    modality_diff_colors = colors;
    modality_diff_colors(3,:) = [0 0 0];
    
    p_critical = 0.01;
    
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
            
            if sum(ps < p_critical) > 0
                plot(rate_ts{j}(ps<p_critical),max(ylim)*.9,'.','color',modality_diff_colors(k,:));
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
            
            if sum(ps < p_critical)>0
                plot(rate_ts{j}(ps < p_critical),max(ylim)*.9,'.','color',modality_diff_colors(k,:));
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
    select_for_tuning = select_bottom_line;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    find_for_tuning = find(select_for_tuning);
    
    % Time centers (of CP time windows, width = 500 ms)
    t_stim_center = mean(time_markers{j}(1,1:2)) + 150; % Stim center + sensory delay
    t_pre_center = mean(time_markers{j}(1,3)) - group_result(end).mat_raw_PSTH.binSize_CP/2;   % Pre-sac epoch
    t_post_center = mean(time_markers{j}(1,3)) + group_result(end).mat_raw_PSTH.binSize_CP/2;  % Post-sac epoch
    
    [~,center_t_ind] = min(abs(CP_ts{j} - t_stim_center));
    [~,pre_t_ind] = min(abs(CP_ts{j} - t_pre_center));
    [~,post_t_ind] = min(abs(CP_ts{j} - t_post_center));
    
    tuning_time_phase = [center_t_ind pre_t_ind post_t_ind];
    tuning_time_phase_title = {'Stimulus', 'Pre-sac', 'Post-sac'};
    
    unique_heading = group_result(find_for_tuning(end)).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters for tuning curve normalization
    
    time_for_dynamic_range = find(CP_ts{j} < time_markers{j}(1,2)); % all pre-sac times
    % time_for_dynamic_range = [center_t_ind - 5 center_t_ind + 5]; % Around stimulus center
    modalities_share_same_dynamic_range = 1;
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
    tuning_mean_correctonly = tuning_mean_all;
    tuning_sem_correctonly = tuning_mean_all;
    
    for pp = 1:length(CP_ts{j})
        
        for k = 1:3
            % Mean and sem
            this_tuning_all = tuning_pack{1}{k,pp}(:,~isnan(tuning_pack{1}{k,pp}(1,:)));
            tuning_mean_all(:,pp,k) = mean(this_tuning_all,2);
            tuning_sem_all(:,pp,k) = std(this_tuning_all,[],2)/sqrt(size(this_tuning_all,2));
            
            this_tuning_correctonly =  tuning_pack{2}{k,pp}(:,~isnan(tuning_pack{2}{k,pp}(1,:)));
            tuning_mean_correctonly(:,pp,k) = mean(this_tuning_correctonly,2);
            tuning_sem_correctonly(:,pp,k) = std(this_tuning_correctonly,[],2)/sqrt(size(this_tuning_correctonly,2));
        end
        
    end
    
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
            plot(unique_heading,tuning_mean_correctonly(:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
            h = errorbar(unique_heading,tuning_mean_correctonly(:,tuning_time_phase(pp),k),tuning_sem_correctonly(:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2);
            errorbar_tick(h,10000);
            title([tuning_time_phase_title{pp} ', n = ' num2str(size(this_tuning_correctonly,2))]);
            axis tight; xlim(xlim*1.1);
            
            SetFigure(15);
        end
    end
    
    %% Tuning: Animation
    
    % Plotting tuning curves at three time points
    set(figure(2099+figN),'pos',[9 309 788 650],'name',['Tuning curve, j = ' num2str(j)]); clf; figN = figN+1;
    ylabel('Correct only');
    SetFigure(15); drawnow;
    
    % Preallocate movie structure.
    mov_tuning(1:length(CP_ts{j})) = struct('cdata', [], 'colormap', []);
    
    for pp = 1:length(CP_ts{j})
        hold off;
        for k = 1:3  % For each stim type
            % Plotting
            plot(unique_heading,tuning_mean_correctonly(:,pp,k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
            hold on;
            h = errorbar(unique_heading,tuning_mean_correctonly(:,pp,k),tuning_sem_correctonly(:,pp,k),'color',colors(k,:),'LineWid',2);
            errorbar_tick(h,10000);
        end
        axis([-7 7 0.25 0.8]);
        title(CP_ts{j}(pp));
        axis tight;
        
        %     mov_tuning(pp) = getframe(gcf);
        
        drawnow;
        pause(0.02);
    end
    
    % movie2avi(mov_tuning,'LIP_tuning_evolve_j_1_3.avi','compression','none');
    
    %% Tuning: Hotgram
    to_evolve = {tuning_mean_all tuning_mean_correctonly};
    to_evolve_title = {'All' 'Correct only'};
    
    for co = 1:2
        set(figure(2099+figN),'pos',[82 99 927 855],'name','Evolve of tuning curves'); clf; figN = figN+1;
        
        [X,Y] = meshgrid(CP_ts{j},unique_heading);
        [Xq,Yq] = meshgrid(linspace(min(CP_ts{j}),max(CP_ts{j}),1000),linspace(min(unique_heading),max(unique_heading),100));
        
        
        for k = 1:3  % For each stim type
            subplot_tight(1,3,k,0.02,0.1);
            surf(Xq,Yq,interp2(X,Y,to_evolve{co}(:,:,k),Xq,Yq),'edgecolor','none'); view(90,-90); hold on;
            
            if k>1 ;set(gca,'xtick',[]); end
            if k==2 ; title(to_evolve_title{co}); end
            
            caxis([min(to_evolve{co}(:)) max(to_evolve{co}(:))]); % Use the same color range
            axis tight;
            
            for tt = 1:3
                plot3([1 1] * time_markers{j}(1,tt),ylim,-3*[1 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
            end
            
            axis tight;
            colormap(hsv);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
            
        end
        
        SetFigure();
    end
    
    %% Targ first vs Target Last (All SU + MU)
    
    j = 1;
    selectCells = (xls_num{1}(:,header.Units_RealSU) >=0 ) & (xls_num{1}(:,header.HD_TargFirst)~=0) & (~isnan(ChoiceDiv_All{j}(:,1,k)));
    
    figure(2099+figN); clf; figN = figN+1;
    
    for k = 1:3
        ys = nanmean(ChoiceDiv_All{j}(selectCells,:,k));
        errors = nanstd(ChoiceDiv_All{j}(selectCells,:,k))/sqrt(sum(selectCells));
        h1 = shadedErrorBar(rate_ts{j},ys,errors * 0,{'Color',colors(k,:)},transparent);
        set(h1.mainLine,'LineWidth',2);
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
        h2 = shadedErrorBar(rate_ts{j},ys,errors * 0,{'Linestyle','-.','Color',colors(k,:) * 0.5 + [.5 .5 .5] },transparent);
        set(h2.mainLine,'LineWidth',2)
        hold on;
    end
    
    str2 = sprintf('Choice divergence (Target Last, N = %g)',sum(selectCells));
    
    for tt = 1:3
        plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
    end
    
    SetFigure(); axis tight;
    
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
    
    keyboard;
    
catch err
    err
    opentoline(err.stack(1).file,err.stack(1).line);
    keyboard
end

    % Show individual cell selected from the figure. HH20150424
    function Show_individual_cell(~,~,h_line, select_for_this)  % Note I use a nested function here
        persistent h_marker;
        
        % Recover the cell number
        pos = get(gca,'currentPoint'); posX = pos(1,1); posY = pos(1,2);
        allX = get(h_line,'xData');
        allY = get(h_line,'yData');
        [min_dis,ind] = min(abs(((posX-allX)/range(xlim)).^2+((posY-allY)/range(ylim)).^2));
        if min_dis > (range(xlim)^2+range(ylim)^2)/100 ; return; end
        
        real_cell_no = find(cumsum(select_for_this)==ind,1);
        
        % Plotting
        if ~isempty(h_marker) ; try delete(h_marker); catch ; end ;end
        h_marker = plot(allX(ind),allY(ind),'kx','markersize',15,'linew',2);
        
        j_this = 1;
        
        ys_this = group_result(real_cell_no).mat_raw_PSTH.PSTH{j_this,1,1}.ys';
        SeriesComparison(shiftdim(ys_this,-1),rate_ts{j_this},'Colors',{'b','b','r','r',[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'figN',1463);
        axis tight; legend off;
        set(gcf,'pos',[-720 321 714 525]);
        
        for t = 1:3
            plot([1 1] * time_markers{j_this}(1,t),ylim,'k','linestyle',marker_for_time_markers{j_this}{t},'linew',1.5);
        end
        
        text(min(xlim),max(ylim),sprintf('X = %g, Y = %g, real cell# = %g, %s',...
            posX,posY,real_cell_no,group_result(real_cell_no).cellID{1}(33:end)),'fontsize',10);
        xlabel('t'); ylabel('Firing rate');
        xlim(xlim * 0.9);
        
    end

    % Play population dynamics. HH20150424
    function Population_dynamics(~,~,flag)
        switch flag
            case {0,10}   % Play & record
                % Preallocate movie structure.
                mov(1:size(projPC{1},2)) = struct('cdata', [],...
                    'colormap', []);
                
                for t_ind = 1 : size(projPC{1},2)
                    for kk = 1:3
                        % Update
                        set(h_pref(kk),'xdata',projPC{which_two_dimension(1)}((kk-1)*2+1,t_ind),'ydata',projPC{which_two_dimension(2)}((kk-1)*2+1,t_ind));
                        set(h_null(kk),'xdata',projPC{which_two_dimension(1)}(kk*2,t_ind),'ydata',projPC{which_two_dimension(2)}(kk*2,t_ind));
                    end
                    
                    %     set(h_text,'string',num2str(fix(rate_ts{j}(tt)/10)*10));
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
                    set(h_pref(kk),'xdata',projPC{which_two_dimension(1)}((kk-1)*2+1,t_ind),'ydata',projPC{which_two_dimension(2)}((kk-1)*2+1,t_ind));
                    set(h_null(kk),'xdata',projPC{which_two_dimension(1)}(kk*2,t_ind),'ydata',projPC{which_two_dimension(2)}(kk*2,t_ind));
                end
                
                %     set(h_text,'string',num2str(fix(rate_ts{j}(tt)/10)*10));
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
%                     %     set(h_text,'string',num2str(fix(rate_ts{j}(tt)/10)*10));
%                     set(h_timeline,'xdata',PCA_B_times(curr_t_ind)*[1 1]);
%                 end
        end
    end


end

















