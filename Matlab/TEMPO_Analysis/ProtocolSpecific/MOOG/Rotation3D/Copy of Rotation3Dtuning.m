function Rotation3Dtuning(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_tri = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_tri ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_tri) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');


condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition_num = munique(condition_num');

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));
% added by Katsu 111606
spon_std = std(temp_spike_rates(spon_found))
spon_ste = spon_std/sqrt(length(find(spon_found)))

% plot option, if regular plot, set to 0, if lambert plot, set to 1
% lamber_plot = 0;  % regular plot with elevation in a linear step
lamber_plot = 1; % lambert plot with elevation in a sin-transformed way


% -------------------------------------------------------------------------
% % ANOVA to test tuning significance
% resp_mat_anova = [];
% for k=1: length(unique_condition_num)
%     n = 0;
%     for i=1:length(unique_azimuth)
%         for j=1:length(unique_elevation)
%             select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
%             if (sum(select) > 0)
%                 n = n+1;
%                 resp_mat_anova{k}(:,n) = spike_rates(select)';
%             end
%         end
%     end
%     [p_anova, table, stats] = anova1(resp_mat_anova{k},[],'off');
%     P_anova(k) = p_anova;
%     anova_table{k} = table;
%     F_val(k) = anova_table{k}(2,5);
% end
% F_val = cell2mat(F_val);


% -------------------------------------------------------------------------
%ANOVA modified by Aihua, it does not require whole trials, it does not matter if trial stopped during repetition
trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-14) * length(unique_condition_num) + 1;
repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);

% first parse raw data into repetitions, including null trials
for q = 1:repetitions
   azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   elevation_rep{q} = temp_elevation(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   condition_num_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%% 111706Katsu TEMPORALY ANOVA+Mean std +ste Analysis
% resp_mat_anova90 = [];% Forward
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% %      
%      %unique_azimuth(3)=90 deg, unique_elevation(3)=0 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(3) & elevation_rep{q}==unique_elevation(3) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% %                   
%                    resp_mat_anova90{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999];
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anova90{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anova90, table90, stats] = anova1(resp_mat_anova90{k},[],'off');
%    P_anova90(k) = p_anova90
%    anova_table90{k} = table90;
%    F_val90(k) = anova_table90{k}(2,5);
% 
% clear select_cardinal;
% select_cardinal=logical(azimuth==unique_azimuth(3) & elevation==unique_elevation(3) & condition_num==unique_condition_num(k) );
% if (sum(select_cardinal) > 0)
%     resp90(k)=mean(spike_rates(select_cardinal))
%     std90(k)=std(spike_rates(select_cardinal))
%     ste90(k)=std(spike_rates(select_cardinal))/sqrt(repetitions)
% end
% end
% F_val90 = cell2mat(F_val90);
% 
% resp_mat_anova270 = [];% Backward
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% %     
%      %unique_azimuth(7)=270 deg, unique_elevation(3)=0 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(7) & elevation_rep{q}==unique_elevation(3) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% %                   
%                    resp_mat_anova270{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999];
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anova270{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anova270, table270, stats] = anova1(resp_mat_anova270{k},[],'off');
%    P_anova270(k) = p_anova270
%    anova_table270{k} = table270;
%    F_val270(k) = anova_table270{k}(2,5);
%    
% clear select_cardinal;
% select_cardinal=logical(azimuth==unique_azimuth(7) & elevation==unique_elevation(3) & condition_num==unique_condition_num(k) );
% if (sum(select_cardinal) > 0)
%     resp270(k)=mean(spike_rates(select_cardinal))
%     std270(k)=std(spike_rates(select_cardinal))
%     ste270(k)=std(spike_rates(select_cardinal))/sqrt(repetitions)
% end
% end
% F_val270 = cell2mat(F_val270);
% 
% 
% 
% resp_mat_anovaleft = []; %Left
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% 
%      %unique_azimuth(5)=180 deg, unique_elevation(3)=0 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(5) & elevation_rep{q}==unique_elevation(3) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% 
%                    resp_mat_anovaleft{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999];
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anovaleft{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anovaleft, tableleft, stats] = anova1(resp_mat_anovaleft{k},[],'off');
%    P_anovaleft(k) = p_anovaleft
%    anova_tableleft{k} = tableleft;
%    F_valleft(k) = anova_tableleft{k}(2,5);
% 
% clear select_cardinal;
% select_cardinal=logical(azimuth==unique_azimuth(5) & elevation==unique_elevation(3) & condition_num==unique_condition_num(k) );
% if (sum(select_cardinal) > 0)
%     respleft(k)=mean(spike_rates(select_cardinal))
%     stdleft(k)=std(spike_rates(select_cardinal))
%     steleft(k)=std(spike_rates(select_cardinal))/sqrt(repetitions)
% end
% end
% F_valleft = cell2mat(F_valleft);
% 
% resp_mat_anovaright = []; %Right
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% %    
%      %unique_azimuth(1)=0 deg, unique_elevation(3)=0 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(1) & elevation_rep{q}==unique_elevation(3) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% 
%                    resp_mat_anovaright{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999];
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anovaright{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anovaright, tableright, stats] = anova1(resp_mat_anovaright{k},[],'off');
%    P_anovaright(k) = p_anovaright
%    anova_tableright{k} = tableright;
%    F_valright(k) = anova_tableright{k}(2,5);
% 
% clear select_cardinal;
% select_cardinal=logical(azimuth==unique_azimuth(1) & elevation==unique_elevation(3) & condition_num==unique_condition_num(k) );
% if (sum(select_cardinal) > 0)
%     respright(k)=mean(spike_rates(select_cardinal))
%     stdright(k)=std(spike_rates(select_cardinal))
%     steright(k)=std(spike_rates(select_cardinal))/sqrt(repetitions)
% end
% end
% F_valright = cell2mat(F_valright);
% 
% 
% 
% resp_mat_anovaup = []; %Up
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% 
%      %unique_azimuth(1)=0 deg, unique_elevation(1)=-90 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(1) & elevation_rep{q}==unique_elevation(1) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% 
%                    resp_mat_anovaup{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999];
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anovaup{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anovaup, tableup, stats] = anova1(resp_mat_anovaup{k},[],'off');
%    P_anovaup(k) = p_anovaup
%    anova_tableup{k} = tableup;
%    F_valup(k) = anova_tableup{k}(2,5);
% 
% clear select_cardinal;
% select_cardinal=logical(azimuth==unique_azimuth(1) & elevation==unique_elevation(1) & condition_num==unique_condition_num(k) );
% if (sum(select_cardinal) > 0)
%     respup(k)=mean(spike_rates(select_cardinal))
%     stdup(k)=std(spike_rates(select_cardinal))
%     steup(k)=std(spike_rates(select_cardinal))/sqrt(repetitions)
% end
% end
% F_valup = cell2mat(F_valup);
% 
% resp_mat_anovadown = []; % Down
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% %       
%      %unique_azimuth(1)=0 deg, unique_elevation(5)=90 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(1) & elevation_rep{q}==unique_elevation(5) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% 
%                    resp_mat_anovadown{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999];
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anovadown{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anovadown, tabledown, stats] = anova1(resp_mat_anovadown{k},[],'off');
%    P_anovadown(k) = p_anovadown
%    anova_tabledown{k} = tabledown;
%    F_valdown(k) = anova_tabledown{k}(2,5);
% 
% clear select_cardinal;
% select_cardinal=logical(azimuth==unique_azimuth(1) & elevation==unique_elevation(5) & condition_num==unique_condition_num(k) );
% if (sum(select_cardinal) > 0)
%     respdown(k)=mean(spike_rates(select_cardinal))
%     stddown(k)=std(spike_rates(select_cardinal))
%     stedown(k)=std(spike_rates(select_cardinal))/sqrt(repetitions)
% end
% end
% F_valdown = cell2mat(F_valdown);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Also, write out some summary data to a cumulative summary file
% buff = sprintf('%s\t %4.4f\t  %4.4f\t   %4.4f\t   %4.4f\t   %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t   %4.4f\t  %4.4f\t   %4.4f\t   %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t   %4.4f\t  %4.4f\t   %4.4f\t   %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t   %4.4f\t  %4.4f\t   %4.4f\t   %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t   %4.4f\t  %4.4f\t   %4.4f\t   %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t   %4.4f\t  %4.4f\t   %4.4f\t   %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t   %4.4f\t  %4.4f\t   %4.4f\t   %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t   %4.4f\t  %4.4f\t   %4.4f\t   %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t  %4.4f\t   %4.4f\t  ', ...
%     FILE, spon_resp,  spon_std,  spon_ste, P_anovaleft, respleft, stdleft, steleft, P_anovaright, respright, stdright, steright,...
%      P_anovaup, respup, stdup,  steup, P_anovadown, respdown, stddown, stedown,...
%      P_anova90, resp90, std90,  ste90, P_anova270, resp270, std270, ste270);
%    
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Katsu_Ro6Cardinal_ste_060807.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t SPon\t Std\t  Ste\t  leftVeb_anovaP\t leftVis_anovaP\t leftVeb_mean\t leftVis_mean\t leftVeb_std\t leftVis_std\t leftVeb_ste\t leftVis_ste\t  rightVeb_anovaP\t rightVis_anovaP\t rightVeb_mean\t rightVis_mean\t rightVeb_std\t rightVis_std\t rightVeb_ste\t rightVis_ste\t upVeb_anovaP\t upVis_anovaP\t upVeb_mean\t upVis_mean\t upVeb_std\t upVis_std\t upVeb_ste\t upVis_ste\t downVeb_anovaP\t downVis_anovaP\t downVeb_mean\t downVis_mean\t downVeb_std\t downVis_std\t downVeb_ste\t downVis_ste\t 90Veb_anovaP\t 90Vis_anovaP\t 90Veb_mean\t 90Vis_mean\t 90Veb_std\t 90Vis_std\t   90Veb_ste\t 90Vis_ste\t  270Veb_anovaP\t 270Vis_anovaP\t 270Veb_mean\t 270Vis_mean\t 270Veb_std\t 270Vis_std\t  270Veb_ste\t 270Vis_ste\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);

%% 060807 Commented by Katsu

%%%%%%%%%%%%%%%%%%%%%%%%%% 090506Katsu TEMPORALY ROLL Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resp_mat_anova90 = [];
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% %        n = 0;
%      %unique_azimuth(3)=90 deg, unique_elevation(3)=0 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(3) & elevation_rep{q}==unique_elevation(3) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% %                    n = n+1;
% %                    resp_mat_anova{k}(q,n) = spike_rates_rep{q}(select_rep{q})';
%                    resp_mat_anova90{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999]
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anova90{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anova90, table90, stats] = anova1(resp_mat_anova90{k},[],'off');
%    P_anova90(k) = p_anova90
%    anova_table90{k} = table90;
%    F_val90(k) = anova_table90{k}(2,5);
% end
% F_val90 = cell2mat(F_val90);
% 
% resp_mat_anova270 = [];
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% %        n = 0;
%      %unique_azimuth(7)=270 deg, unique_elevation(3)=0 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(7) & elevation_rep{q}==unique_elevation(3) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% %                    n = n+1;
% %                    resp_mat_anova{k}(q,n) = spike_rates_rep{q}(select_rep{q})';
%                    resp_mat_anova270{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999]
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anova270{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anova270, table270, stats] = anova1(resp_mat_anova270{k},[],'off');
%    P_anova270(k) = p_anova270
%    anova_table270{k} = table270;
%    F_val270(k) = anova_table270{k}(2,5);
% end
% F_val270 = cell2mat(F_val270);
% 
% % Also, write out some summary data to a cumulative summary file
% buff = sprintf('%s\t %4.2f\t  %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t   %2.4f\t', ...
%     FILE, spon_resp, F_val90, P_anova90, F_val270, P_anova270);
%     %      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot,  p{:} , resp_std, gain, F_val, P_anova);
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Katsu_Roll_Rotation_ANOVA.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t SPon\t 90Veb_F_anova\t 90Vis_F_anova\t 90Veb_anovaP\t 90Vis_anovaP\t 270Veb_F_anova\t 270Vis_F_anova\t 270Veb_anovaP\t 270Vis_anovaP\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%              090706Katsu TEMPORALY       Yaw
%%%%%%%%%%%%%%%%%%%%%%%%%              Pitch Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resp_mat_anovaleft = [];
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% 
%      %unique_azimuth(1)=0 deg, unique_elevation(1)=-90 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(1) & elevation_rep{q}==unique_elevation(1) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% 
%                    resp_mat_anovaleft{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999];
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anovaleft{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anovaleft, tableleft, stats] = anova1(resp_mat_anovaleft{k},[],'off');
%    P_anovaleft(k) = p_anovaleft
%    anova_tableleft{k} = tableleft;
%    F_valleft(k) = anova_tableleft{k}(2,5);
% end
% F_valleft = cell2mat(F_valleft);
% 
% resp_mat_anovaright = [];
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% %        n = 0;
%      %unique_azimuth(1)=0 deg, unique_elevation(5)=90 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(1) & elevation_rep{q}==unique_elevation(5) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% 
%                    resp_mat_anovaright{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999]
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anovaright{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anovaright, tableright, stats] = anova1(resp_mat_anovaright{k},[],'off');
%    P_anovaright(k) = p_anovaright
%    anova_tableright{k} = tableright;
%    F_valright(k) = anova_tableright{k}(2,5);
% end
% F_valright = cell2mat(F_valright);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resp_mat_anovaup = [];
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% 
%      %unique_azimuth(1)=0 deg, unique_elevation(3)=0 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(1) & elevation_rep{q}==unique_elevation(3) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% 
%                    resp_mat_anovaup{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999]
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anovaup{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anovaup, tableup, stats] = anova1(resp_mat_anovaup{k},[],'off');
%    P_anovaup(k) = p_anovaup
%    anova_tableup{k} = tableup;
%    F_valup(k) = anova_tableup{k}(2,5);
% end
% F_valup = cell2mat(F_valup);
% 
% resp_mat_anovadown = [];
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
% %        n = 0;
%      %unique_azimuth(5)=180 deg, unique_elevation(3)=0 deg
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(5) & elevation_rep{q}==unique_elevation(3) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
% 
%                    resp_mat_anovadown{k}(q,2) = spike_rates_rep{q}(select_rep{q})';
%                end  
%    end
%    clear select_spo;spo=[-9999]
%    for q=1:1:repetitions
%                select_spo{q} = logical( azimuth_rep{q}==spo & elevation_rep{q}==spo);
%                if (sum(select_spo{q}) > 0)
% %                   
%                    resp_mat_anovadown{k}(q,1) = spike_rates_rep{q}(select_spo{q})';
%                end  
%    end
% 
% 
%    [p_anovadown, tabledown, stats] = anova1(resp_mat_anovadown{k},[],'off');
%    P_anovadown(k) = p_anovadown
%    anova_tabledown{k} = tabledown;
%    F_valdown(k) = anova_tabledown{k}(2,5);
% end
% F_valdown = cell2mat(F_valdown);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Also, write out some summary data to a cumulative summary file
% buff = sprintf('%s\t %4.2f\t  %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t   %2.4f\t', ...
%     FILE, spon_resp, F_valleft, P_anovaleft, F_valright, P_anovaright,  F_valup, P_anovaup, F_valdown, P_anovadown);
%    
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Katsu_YawPitch_Rotation_ANOVA.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t SPon\t leftVeb_F_anova\t leftVis_F_anova\t leftVeb_anovaP\t leftVis_anovaP\t rightVeb_F_anova\t rightVis_F_anova\t rightVeb_anovaP\t rightVis_anovaP\t upVeb_F_anova\t upVis_F_anova\t upVeb_anovaP\t upVis_anovaP\t downVeb_F_anova\t downVis_F_anova\t downVeb_anovaP\t downVis_anovaP\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Comment 090506 %%%%%%%%%%%Uncomment 090606

%%%%%%%%%Comment 090706 %%%%%%%%%%%Uncomment 090806

%%%%%%%%%Comment 111406 %%%%%%%%%%%Uncomment 111406
resp_mat_anova = [];
for k=1: length(unique_condition_num)
   clear select_rep;
   for q=1:1:repetitions
       n = 0;
       for i=1:length(unique_azimuth)
           for j=1:length(unique_elevation)
               select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==unique_elevation(j) & condition_num_rep{q}==unique_condition_num(k) );
               if (sum(select_rep{q}) > 0)
                   n = n+1;
                   resp_mat_anova{k}(q,n) = spike_rates_rep{q}(select_rep{q})';
               end
           end
       end
   end
   [p_anova, table, stats] = anova1(resp_mat_anova{k},[],'off');
   P_anova(k) = p_anova;
   anova_table{k} = table;
   F_val(k) = anova_table{k}(2,5);
end
F_val = cell2mat(F_val);


resp_mat_anova_hor = [];
for k=1: length(unique_condition_num)
   clear select_rep;
   for q=1:1:repetitions
       n = 0;
       for i=1:length(unique_azimuth) 
           select_rep_hor{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==0 & condition_num_rep{q}==unique_condition_num(k) );
           if (sum(select_rep_hor{q}) > 0)
               n = n+1;
               resp_mat_anova_hor{k}(q,n) = spike_rates_rep{q}(select_rep_hor{q})';
           end       
       end
   end
   [p_anova_hor, table_hor, stats_hor] = anova1(resp_mat_anova_hor{k},[],'off');
   P_anova_hor(k) = p_anova_hor;
   anova_table_hor{k} = table_hor;
   F_val_hor(k) = anova_table_hor{k}(2,5);
end
F_val_hor = cell2mat(F_val_hor);



%% ADD CODE HERE FOR PLOTTING
resp_mat = [];
for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=1: length(unique_condition_num)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(spike_rates(select));
                resp_mat_vector(k, j, i) = mean(spike_rates(select));
                for t = 1 : length(spike_rates(select));              % this is to calculate response matrix based on each trial
                    spike_temp = spike_rates(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(spike_rates(select));     % calculate std between trials for later DSI usage
                resp_mat_ste(k, j, i) = resp_mat_std(k, j, i)/ sqrt(length(find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
            else
%                resp_mat_trial{k}(t, j, i) = 0; % From the begining
%                (071306) this sentence was commented by Katsu
%Following was wrong, I think from 07/10/06 by Katsu
%                 resp_mat(k, j, i) = 0;
%                 resp_mat_vector(k, j, i) = resp_mat_vector(k, j, 1);
%                 resp_mat_std(k, j, i) = 0;
%                 resp_mat_ste(k, j, i) = 0;
%So corrected by Katsu 07/15/06
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,i) =0; % for vector sum only % once this value was 1, not i %%011707 by Katsu
                resp_mat_std(k, j, i) = 0;
                resp_mat_ste(k, j, i) = 0;
            end
        end        
    end
end

% % creat a real 3-D based plot where the center correspond to forward and
% % both lateral edges correspond to backward
% resp_mat_tran = [];
% for i=1:( length(unique_azimuth)+1 )        % add a azimuth 360 to make it circularlly continuously
%     for j=1:length(unique_elevation)
%         for k=1:length(unique_condition_num)
%             if ( j == 1 | j==5 )                          % fill NaN point with data
%                 resp_mat_tran(k, j, i) = resp_mat(k, j, 1);   
%             else
%                 if (i < 8 )
%                     resp_mat_tran(k, j, i) = resp_mat(k,j, 8-i);
%                 elseif(i==8)
%                     resp_mat_tran(k, j, i) = resp_mat(k,j, 8);
%                 else
%                     resp_mat_tran(k, j, i) = resp_mat(k,j, 7);
%                 end
%             end
%         end        
%     end
% end
resp_mat_tran(:,:,1) = resp_mat(:,:,7);
resp_mat_tran(:,:,2) = resp_mat(:,:,6);
resp_mat_tran(:,:,3) = resp_mat(:,:,5);
resp_mat_tran(:,:,4) = resp_mat(:,:,4);
resp_mat_tran(:,:,5) = resp_mat(:,:,3);
resp_mat_tran(:,:,6) = resp_mat(:,:,2);
resp_mat_tran(:,:,7) = resp_mat(:,:,1);
resp_mat_tran(:,:,8) = resp_mat(:,:,8);
resp_mat_tran(:,:,9) = resp_mat_tran(:,:,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %output some data  %Aihua
% switch(SpikeChan)
%     case 1
%         dum_name=FILE(1:end-4);
%     case 4
%         dum_name=[FILE(1:end-4),'_MU'];
%     case 5
%         dum_name=[FILE(1:end-4),'n1']
%     case 6
%         dum_name=[FILE(1:end-4),'n2'];
%     case 7
%         dum_name=[FILE(1:end-4),'n3'];
% end
% 
% save_name = 'resp_mat_tran'; % 
% spacer = ' ';
% eval([dum_name '= resp_mat_tran;'])
% OutPath=['C:\Aihua\z_TempOutputs\Rot\'];%OutPath=['Z:\Users\Aihua\PIVC_analysis\datatemp\Rot\'];
% cd(OutPath);
% if (exist([save_name '.mat'], 'file') == 0)    %file does not yet exist 
%     eval(['save ' OutPath save_name spacer dum_name]);   
% else 
%     eval(['save ' OutPath save_name spacer dum_name ' -APPEND']);  
% end
% 
% save_name = 'resp_mat_std'; % 
% spacer = ' ';
% eval([dum_name '= resp_mat_std;'])
% % OutPath=['Z:\Users\Aihua\PIVC_analysis\datatemp\Rot\'];
% % cd(OutPath);
% if (exist([save_name '.mat'], 'file') == 0)    %file does not yet exist 
%     eval(['save ' OutPath save_name spacer dum_name]);   
% else 
%     eval(['save ' OutPath save_name spacer dum_name ' -APPEND']);  
% end
% 
% save_name = 'resp_mat_ste'; % 
% spacer = ' ';
% eval([dum_name '= resp_mat_ste;'])
% % OutPath=['Z:\Users\Aihua\PIVC_analysis\datatemp\Rot\'];
% % cd(OutPath);
% if (exist([save_name '.mat'], 'file') == 0)    %file does not yet exist 
%     eval(['save ' OutPath save_name spacer dum_name]);   
% else 
%     eval(['save ' OutPath save_name spacer dum_name ' -APPEND']);  
% end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate maximum and minimum firing rate
max_res = max(max(max(resp_mat)));
% max_vi_ve = max(max(max(resp_mat(2:3,:,:))));
min_res = min(min(min(resp_mat_tran)));

vector_num = length(unique_azimuth) * (length(unique_elevation)-2) + 2;
repeat = floor( length(spike_rates) / vector_num );

% Define figure
xoffset=0;
yoffset=0;
figure(2);clf;
set(2,'Position', [5,15 980,650], 'Name', 'Rotation_3D');
orient landscape;
%set(0, DefaultAxesXTickMode, 'manual', DefaultAxesYTickMode, 'manual', 'DefaultAxesZTickMode', 'manual');
axis off;
% for cosine plot(Added by Aihua, 01-24-2007)
azi_cos = [1,2,3,4,5,6,7,8,9];
ele_sin = [-1,-0.707,0,0.707,1];


for k=1: length(unique_condition_num) 
    
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.4;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
    
    if lamber_plot ==1
         contourf( azi_cos, ele_sin, squeeze( resp_mat_tran(k,:,:)) );
    elseif lamber_plot ==0
         contourf( squeeze( resp_mat_tran(k,:,:)) ); 
    end
    
    % set the same scale for visual and combined conditions but here assuming vestibular response is always smaller than that in visual and
    % combined conditions
%     if ( k==2 | k==3 )
%     caxis([min_res, max_res]);
%    caxis([0, 18]);
%     end
% if ( k==1 )           % Katsu for paper settle for m3c294
% caxis([25, 80]);
% end
% if ( k==2 )
%caxis([10, 100]);
% end
% % % adjust this section to change scaling
% if ( k==1 )           % Katsu for paper settle for m3c296
% caxis([10, 100]);
% end
% if ( k==2 )
% caxis([40, 120]);
% end
    colorbar;
    % make 0 correspond to rightward and 180 correspond to leftward
    set(gca, 'ydir' , 'reverse');
    set(gca, 'xtick', [] );
    set(gca, 'ytick', [] );    
    title( h_title{k} );
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Add by Katsu 102406
%             text(5,1,'x270','FontSize',14);
%             text(7,2,'x315','FontSize',14);
%             text(7,3,'x0','FontSize',14);
%             text(7,4,'x45','FontSize',14);
%             text(5,5,'x90','FontSize',14);
%             text(2,4,'135x','FontSize',14);
%             text(2,3,'180x','FontSize',14);
%             text(2,2,'225x','FontSize',14);
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % plot 1-D for mean respond as a function of elevation
    axes('position',[0.06+xoffset 0.54+yoffset 0.04 0.24]);
    for j=1:length(unique_elevation)
        y_elevation_mean(1,j)=mean(resp_mat_tran(k,j,:));
        y_elevation_std(1,j) =std( spike_rates([find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )]) );
        y_elevation_ste(1,j) =y_elevation_std(1,j)/ sqrt(length(find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
    end
    if lamber_plot == 1
        x_elevation=[-1,-0.707,0,0.707,1]; %uncomment for cosine plot
    elseif lamber_plot ==0
        x_elevation=unique_elevation;%05/22/06 Katsu changed not cosin axis
    end
    errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'ko-');%-----------Temporaly disable
    
    xlabel('Elevation');
    view(90,90);
    set(gca, 'xtick',unique_elevation);
            xlim([-90, 90]);
    if lamber_plot == 1
       xlim([-1, 1]);
       set(gca, 'XTickMode','manual');
       set(gca, 'xtick',[-1,-0.707,0,0.707,1]);
       set(gca, 'xticklabel','-90|-45|0|45|90'); 
    elseif lamber_plot ==0
       xlim([-90, 90]);
    end
    ylim([min(y_elevation_mean(1,:))-max(y_elevation_ste(1,:)), max(y_elevation_mean(1,:))+max(y_elevation_ste(1,:))]);%----------Now disable temporaly

% plot 1-D for mean respond as a function of azimuth
    axes('position',[0.11+xoffset 0.46+yoffset 0.274 0.06]);
    for i=1:(length(unique_azimuth) )
        y_azimuth_mean(1,i)=mean(resp_mat_tran(k,:,i));
        y_azimuth_std(1,i) =std( spike_rates([find( (azimuth==unique_azimuth(i))&(condition_num==unique_condition_num(k)) )]) );
        y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (azimuth==unique_azimuth(i))&(condition_num==unique_condition_num(k)) )) );    
    end
    y_azimuth_mean(1,9) = mean(resp_mat_tran(k,:,1));
    for i=1:( length(unique_azimuth)+1 )
        if (i < 8)        
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8-i);
        elseif (i == 8)
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8);
        else
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,7);
        end
    end
    x_azimuth=1:(length(unique_azimuth)+1);
    errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'ko-');%----------Now disable temporaly
%     errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'k-');% Katsu for paper settle for m3c294
    xlim( [1, length(unique_azimuth)+1] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); 
    xlabel('Azimuth');
    ylim([min(y_azimuth_mean(1,:))-max(y_azimuth_ste(1,:)), max(y_azimuth_mean(1,:))+max(y_azimuth_ste(1,:))]);%----------Now disable temporaly

    xoffset=xoffset+0.48;
    
    % calculate min and max firing rate, standard deviation, HTI, Vectorsum
    Min_resp(k) = min( min( resp_mat_tran(k,:,:)) );
    Max_resp(k) = max( max( resp_mat_tran(k,:,:)) );
    resp_std(k) = sum( sum(resp_mat_std(k,:,:)) ) / vector_num;  % notice that do not use mean here, its (length(unique_azimuth)*length(unique_elevation)-14) vectors intead of 40
    M=squeeze(resp_mat_vector(k,:,:));     % notice that here DSI should use resp_temp without 0 value set manually
    % this part is to calculate vestibular gain
    resp_onedim{k} = [M(1,1),M(2,:),M(3,:),M(4,:),M(5,1)]';     % hard-code temperarilly    
    N=squeeze(resp_mat_vector(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually 
    [Azi, Ele, Amp] = vectorsum(N);
    Vec_sum{k}=[Azi, Ele, Amp];
    % Heading Tuning Index
    r(k) = HTI(M,spon_resp);   % call HTI function    
end
% calculate vestibular gain by a*ves=comb-vis
if (length(unique_stim_type) == 3)
   [bb,bint,rr,rint,stats] = regress( (resp_onedim{3}-resp_onedim{2}), [ones(vector_num,1),(resp_onedim{1}-spon_resp)] );    % with offset
   gain = bb(2);
else
   gain = NaN;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % linear sum model Modified by Yong (for Katsu) 06/19/07
% % R(comb) = a1 x R(vest) + a2 x R(visu) + a3  ; This is the model
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Note; This code should be run under 3 combined condition only !!!!
% %%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
% OPTIONS = optimset('fmincon');
% OPTIONS = optimset('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
% 
% yy1 = @(x)sum( (resp_onedim{1}*x(1)+resp_onedim{2}*x(2)+x(3)-resp_onedim{3}).^2 );  %w1*ve+w2*vi
% es1 = [0.5,0.5,0];
% LB1 = [-5,-5,0];
% UB1 = [5,5,100];
% 
% v1 = fmincon(yy1,es1,A,B,Aeq,Beq,LB1,UB1, NONLCON, OPTIONS); % fminsearch   
% v1
% pred = resp_onedim{1}*v1(1)+resp_onedim{2}*v1(2)+v1(3);
% % figure;% to make sure fit (green) is good or not
% % plot(resp_onedim{1},'k-');
% % hold on;
% % plot(resp_onedim{2},'r-');
% % plot(resp_onedim{3},'g-');
% % plot(pred,'b-');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %-------------------------------------------------------------------
%check significance of HTI and calculate p value, do bootstrap at the same time to test value varience
perm_num=1000;
bin = 0.005;
spike_rates_perm = [];
for n=1: perm_num
    % this is permuted based on trials
    for k=1:length(unique_condition_num)   
        spike_rates_pe{k} = spike_rates( find( condition_num==unique_condition_num(k) ) );
        spike_rates_pe{k} = spike_rates_pe{k}( randperm(length(spike_rates_pe{k})) );
    end

    % put permuted data back to spike_rates
    spike_rates_perm(length(spike_rates))=0;
    for k=1:length(unique_condition_num) 
        ii = find(stim_type == unique_stim_type(k));
        spike_rates_perm(ii) = spike_rates_pe{k};
    end
    
    % re-creat a matrix similar as resp_mat              
    resp_vector_perm = [];
    for i=1:length(unique_azimuth)
        for j=1:length(unique_elevation)
            for k=1:length(unique_condition_num)
                select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
                if (sum(select) > 0)
                    resp_mat_perm(k,j,i) = mean(spike_rates_perm(select));
                    resp_mat_perm_std(k,j,i) = std(spike_rates_perm(select));
                else
                    resp_mat_perm(k,j,i) = 0;
                    resp_mat_perm_std(k,j,i) = 0;
                end
            end        
        end
    end
    
    % re-calculate HTI now
    for k=1: length(unique_condition_num)
 %       resp_perm_std(k) = sum( sum(resp_mat_perm_std(k,:,:)) ) / vector_num; 
        M_perm=squeeze(resp_mat_perm(k,:,:));
        r_perm(k,n) = HTI(M_perm, spon_resp);                                  
    end
    % do bootstrap now
    % first to bootstap raw data among trils 
    repetition = 5;
    for k=1:length(unique_stim_type)
        for i=1:length(unique_azimuth)
            for j=1:length(unique_elevation)
                    select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (stim_type==unique_stim_type(k)) );
                    if (sum(select) > 0)
                        spike_select = spike_rates(select);
                        for b=1:repetition    % use 5 repetitions temporarilly, should doesn't matter whether have 5 repetitions or not actually
                            spike_select = spike_select( randperm(length(spike_select)) );
                            spike_bootstrap(b) = spike_select(1);   % always take the first one element
                        end 
                        resp_mat_boot(k, j, i) = mean(spike_bootstrap);
                    else
                        resp_mat_boot(k, j, i) = 0;
                    end       
             end
         end
     end
     % now recalculate values
     for k=1: length(unique_stim_type)   
         Mb=squeeze(resp_mat_boot(k,:,:));
         r_boot(k,n) = HTI(Mb,spon_resp) ; 
         % also calculate the significant different angles between preferred headings, for only vestibualr condition
         [Azi_boot, Ele_boot, Amp_boot] = vectorsum(Mb);
         Vec_sum_boot{k}=[Azi_boot, Ele_boot, Amp_boot];
         Angle_boot(n)=(180/3.14159) * acos( sin(Vec_sum{1}(2)*3.14159/180) * sin(Vec_sum_boot{1}(2)*3.14159/180)  +  cos(Vec_sum_boot{1}(2)*3.14159/180) * sin(Vec_sum_boot{1}(1)*3.14159/180) * cos(Vec_sum{1}(2)*3.14159/180) * sin(Vec_sum{1}(1)*3.14159/180) + cos(Vec_sum_boot{1}(2)*3.14159/180) * cos(Vec_sum_boot{1}(1)*3.14159/180) * cos(Vec_sum{1}(2)*3.14159/180) * cos(Vec_sum{1}(1)*3.14159/180) );
         azi_boot(1,n) = Vec_sum_boot{1}(1);
         ele_boot(1,n) = Vec_sum_boot{1}(2);
     end  
end
% now calculate p value or significant test
x_bin = 0 : bin : 1;
for k = 1 : length(unique_condition_num)
    hist_perm(k,:) = hist( r_perm(k,:), x_bin );  % for permutation
    hist_boot(k,:) = hist( r_boot(k,:), x_bin );  % for bootstrap
    [hist_boot_angle(k,:),x_angle] = hist( Angle_boot(:), 200 );  % for bootstrap, set 200 bins temporarilly
    bin_sum = 0;
    n = 0;
    while ( n < (r(k)/bin) )
          n = n+1;
          bin_sum = bin_sum + hist_perm(k, n);
          p{k} = (perm_num - bin_sum)/ perm_num;    % calculate p value for HTI
    end 
    bin_sum = 0;
    n = 0;
    while ( bin_sum < 0.025*sum( hist_boot(k,:)) )   % define confidential value to be 0.05, now consider one side only which is 0.025 of confidence
          n = n+1;
          bin_sum = bin_sum + hist_boot(k, n);      
          HTI_boot(k) = r(k) - n * bin ;    % calculate what HTI value is thought to be significant different
    end 
%     bin_sum = 0;
%     n = 0;
%     while ( bin_sum < 0.975*sum( hist_boot_angle(k,:)) )   % define confidential value to be 0.05, now consider one side only which is 0.025 of confidence
%           n = n+1;
%           bin_sum = bin_sum + hist_boot_angle(k, n);      
%           Angle_boot(k) = x_angle(n) ;    
%     end 
end
%%------------------------------------------------------------------
% DDI Direction discrimination index   by Katsu 05/18/06
%--------------------------------------------------------------------
%spike_rates_sqrt = sqrt(spike_rates);That is not original, it makes value bigger by 
each_stim_trials=repetitions*(length(unique_azimuth)*length(unique_elevation)-14) % (length(unique_azimuth)*length(unique_elevation)-14) trajectory  -90+90 -45 0 + 45 *0 45 90 .......
%
SSE_term = [];
for k=1:length(unique_stim_type)
    n=0;
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)
            clear select;
            select=logical((azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (stim_type==unique_stim_type(k)));
            if (sum(select) > 0)
                   n = n+1;
                   SSE_term(k,j,i)=sum((spike_rates(select)-mean(spike_rates(select))).^2);
               else
                   SSE_term(k,j,i)=0;
            end
         end
     end
%      SSE_azimth_sum(k,j)=sum(SSE_term(k,j,:))  %it is not correct
%      SSE_total(k)=sum(SSE_azimth_sum(k,:))
     SSE_total(k)=sum(sum(SSE_term(k,:,:)));
     max_min_term(k)=(Max_resp(k)-Min_resp(k))/2;
     var_term(k)=sqrt(SSE_total(k)/(each_stim_trials-n));
     DDI(k)=max_min_term(k)/(max_min_term(k)+var_term(k))
 end
 %---------------------------------------------------------------------------
% Now show vectorsum, DSI, p and spontaneous at the top of figure--------and DDI by Katsu 05/18/06
axes('position',[0.05,0.85, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_condition_num)] );
h_spon = num2str(spon_resp);
text(0, length(unique_condition_num), FILE);
text(10,length(unique_condition_num),'Protocol          Spon    Minimum   Maximum    Azi      Ele       Amp     Std        HTI          HTIerr         p          F-val        p-ANOVA       DDI');
for k=1:length(unique_condition_num) 
    h_text{k}=num2str( [spon_resp, Min_resp(k), Max_resp(k), Vec_sum{k}, resp_std(k), r(k), HTI_boot(k), p{k}, F_val(k), P_anova(k), DDI(k)], 4);
    text(0,length(unique_condition_num)-k,h_title{k});
    text(10,length(unique_condition_num)-k,'Rotation');
    text(20,length(unique_condition_num)-k, h_text{k} );
end

axis off;

% % % % % OUTPUT FILE
%---------------------------------------------------------------------------------------
% %Also, write out some summary data to a cumulative summary file
% buff = sprintf('%s\t %4.2f\t  %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t'); %buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...');
%      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std, gain, F_val, P_anova, DDI );
% sprint_txt = ['%s']; 
% for i = 1 : 200 % this should be large enough to cover all the data that need to be exported
%      sprint_txt = [sprint_txt, ' %4.3f'];    
% end
% % buff = sprintf(sprint_txt, FILE, resp_mat(1,3,:), resp_mat(2,3,:)  );%
% %outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Tuning_mean.dat'];original
% %outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Tuning_cah.dat'];%Aihua
% %outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Tuning_syed.dat'];%syed
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Darkness_Katsu.dat'];%00430Katsu
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_DDI_Katsu.dat'];%051806Katsu
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Darkness_DDI_Katsu.dat'];%051906Katsu
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Que_Labyrhintectomy_DDI_Katsu.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
% %     fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t');%with combined
%     fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
%--------------------------------------------------------------------------
%  Katsu Que Labyrinthectomy
%--------------------------------------------------------------------------
% sprint_txt = ['%s'];
% for i = 1 : 3000 % this should be large enough to cover all the data that need to be exported
%      sprint_txt = [sprint_txt, ' %4.3f'];    
% end
% buff = sprintf(sprint_txt, FILE, azi_boot, ele_boot  );
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_dark1.dat'];
% buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
%      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std, gain, F_val, P_anova, DDI, var_term);
%outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_DDI_cah.dat'];
% outfile = ['Z:\Users\Syed Chowdhury\syed2\Tempo_Out\3D_Rotation_Tuning_Syed.dat'];
% outfile = ['C:\Aihua\z_TempOutputs\Rotation_cah.dat'];%Aihua
% outfile = ['Z:\Users\Ayanna\tempo_output\3D_Rotation_Tuning_Ayanna.dat'];


% % % %saves data to my z:\users folder ab 20 Sept. 2007.
% % % printflag = 0;
% % % if (exist(outfile, 'file') == 0)    %file does not yet exist
% % %     printflag = 1;
% % % end
% % % fid = fopen(outfile, 'a');
% % % if (printflag)
% % % %     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t');
% % % fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t Veb_var_term\t Vis_var_term\t');
% % %     fprintf(fid, '\r\n');
% % % end
% % % fprintf(fid, '%s', buff);
% % % fprintf(fid, '\r\n');
% % % fclose(fid);

output_specifications  %%Tunde 9/28/07
%---------------------------------------------------------------------------------------
%-----------03/06/06----------------------------------------------------------------------------
%-------------Katsu Running Batch File-------Vest and Visu condition-------------------------------------------------------------------
%---------------Change Filename------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
% % Also, write out some summary data to a cumulative summary file
% buff = sprintf('%s\t %4.2f\t  %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t   %2.4f\t', ...
%      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std, gain, F_val, P_anova);
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Katsu_Rotation_Tuning_Sum.dat'];
% 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_anova\t Vis_anova\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
%---------------------------------------------------------------------------------------
% %---------------------------------------------------------------------------------------
%-----------03/07/06----------------------------------------------------------------------------
%-------------Katsu Running Batch File-----This is only Vestibular Condition---------------------------------------------------------------------
%---------------Change Filename------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
% %Also, write out some summary data to a cumulative summary file of vestibular output only!!!
% buff = sprintf('%s\t %4.2f\t  %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t   %2.4f\t', ...
%      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std, gain, F_val, P_anova, DDI);
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\syed_Rotation_Tuning_Sum.dat'];
% outfile=['F:\Syed_Analysis\syed_Rotation_Tuning_Sum.dat']
% 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t SPon\t Veb_min\t Veb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Veb_HTI\t Veb_HTIerr\t Veb_P\t Veb_std\t gain\t Veb_F_anova\t Veb_anova\t Veb_DDI');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% %---------------------------------------------------------------------------------------

% figure(2);saveas(gcf,['C:\Aihua\z_TempOutputs\figures\' FILE(1:end-4) '.png'],'png')
% saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '.fig'])
% close (2);
% SaveTrials(FILE,BegTrial,EndTrial,P_anova_hor);

% [StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
% MOOG_Rotation_PSTH_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% % MOOG_Rotation_TuningStep_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_Rotation_corr_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_Rotation_sep(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

return;