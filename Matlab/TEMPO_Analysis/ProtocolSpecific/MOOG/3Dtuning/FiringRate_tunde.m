function FiringRate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);

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


condition = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition = munique(condition');

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));

%% ADD CODE HERE FOR PLOTTING
if length(unique_azimuth)==10
    unique_azimuth0=[0:45:315]';
else
    unique_azimuth0=unique_azimuth;
end
resp_mat = [];
for i=1:length(unique_azimuth0)
    for j=1:length(unique_elevation)
        for k=1: length(unique_condition)
            %select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition==unique_condition(k)) );
            select = logical( (azimuth==unique_azimuth0(i)) & (elevation==unique_elevation(j)) & (condition==unique_condition(k)) );%AC 12-11-2007
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(spike_rates(select));
                resp_std(k, j, i) = std(spike_rates(select));
                resp_std_root(k, j, i) = std(sqrt(spike_rates(select)));
            else
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_std(k, j, i) = resp_std(k,j,1);
                resp_std_root(k, j, i) = resp_std_root(k, j, 1);
            end
        end        
    end
end

% % z-score data 
% for k = 1:length(unique_stim_type)
%     for i = 1:length(unique_azimuth)
%         for j=1:length(unique_elevation)
%             select = logical( (azimuth==unique_azimuth(i)) & (stim_type == unique_stim_type(k)) & (elevation==unique_elevation(j)) ) ;  
%             select2 = logical( (azimuth==unique_azimuth(1)) & (stim_type == unique_stim_type(k)) & (elevation==unique_elevation(j)) ) ;  
%             if (sum(select) > 0) 
%                 z_dist = spike_rates(select);
%                 z_dist = (z_dist - mean(z_dist))/std(z_dist);
%                 Z_Spikes(select) = z_dist;
%             end
%         end
%     end
% end

% calculate anova
% trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-2*(length(unique_azimuth)-1)) * length(unique_condition) + 1;
% repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);
repetitions = floor( (EndTrial-(BegTrial-1)) / 79);

resp_mat_anova = [];
for k=1: length(unique_condition) 
    n=0;
    for i=1:length(unique_azimuth)
%         for j=1:length(unique_elevation)
%             select_rep = find( azimuth==unique_azimuth(i) & elevation==unique_elevation(j) & condition==unique_condition(k) );
%             select_rep_horizontal = find( azimuth==unique_azimuth(i) & elevation==0 & condition==unique_condition(k) );
%             if (length(select_rep) > 0)    
%                 n = n+1;            
%                 for q=1:repetitions
%                    resp_mat_anova{k}(q,n) = spike_rates(select_rep(q));
%                    resp_mat_anova_horizontal{k}(q,n) = spike_rates(select_rep_horizontal(q));
%                 end
%             end
%         end
        trial_select = logical( (azimuth==unique_azimuth(i)) & elevation==0 & (condition==unique_condition(k)) );
        for jj = 1 : repetitions; 
            spike_temp = spike_rates(trial_select);   
            resp_horizontal_trial{k}(jj, i) = spike_temp( jj );  
        end        
   end
%    [p_anova, table, stats] = anova1(resp_mat_anova{k},[],'off');
%    P_anova(k) = p_anova;
%    [p_anova, table, stats] = anova1(resp_mat_anova_horizontal{k},[],'off');
%    P_anova_horizontal(k) = p_anova;
end
% deal different conditions
if length(unique_condition) ==1
   visfind = 1;
elseif length(unique_condition) ==2
   visfind = 2; 
elseif length(unique_condition) ==3
   visfind = 2; 
end
%ves(1:26) = [resp_mat(1,1,1),squeeze(resp_mat(1,2,:))',squeeze(resp_mat(1,3,:))',squeeze(resp_mat(1,4,:))',resp_mat(1,5,1)];
% vis(1:26) = [resp_mat(visfind,1,1),squeeze(resp_mat(visfind,2,:))',squeeze(resp_mat(visfind,3,:))',squeeze(resp_mat(visfind,4,:))',resp_mat(visfind,5,1)];
% com(1:26) = [resp_mat(3,1,1),squeeze(resp_mat(3,2,:))',squeeze(resp_mat(3,3,:))',squeeze(resp_mat(3,4,:))',resp_mat(3,5,1)];
%ves_std(1:26) = [resp_std(1,1,1),squeeze(resp_std(1,2,:))',squeeze(resp_std(1,3,:))',squeeze(resp_std(1,4,:))',resp_std(1,5,1)];
% vis_std(1:26) = [resp_std(visfind,1,1),squeeze(resp_std(visfind,2,:))',squeeze(resp_std(visfind,3,:))',squeeze(resp_std(visfind,4,:))',resp_std(visfind,5,1)];
%com_std(1:26) = [resp_std(3,1,1),squeeze(resp_std(3,2,:))',squeeze(resp_std(3,3,:))',squeeze(resp_std(3,4,:))',resp_std(3,5,1)];
%vesout = resp_mat(1, 3, :);
%visout = resp_mat(2, 3, :);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 5000 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
% buff = sprintf(sprint_txt, FILE, Z_Spikes );
%buff = sprintf(sprint_txt, FILE, spon_resp, vis, vis_std, P_anova(visfind)  );
%buff = sprintf(sprint_txt, FILE,  resp_mat(1,3,:),resp_mat(2,3,:),resp_std(1,3,:),resp_std(2,3,:),resp_std_root(1,3,:),resp_std_root(2,3,:),P_anova_horizontal(1),P_anova_horizontal(2)  );
buff = sprintf(sprint_txt, FILE,  resp_mat(1,3,:),resp_std(1,3,:),resp_std_root(1,3,:),resp_mat(2,3,:),resp_std(2,3,:),resp_std_root(2,3,:));
% buff = sprintf(sprint_txt, FILE,repetitions,resp_horizontal_trial{1}(:,:),resp_horizontal_trial{2}(:,:),resp_horizontal_trial{3}(:,:) );
% outfile = ['Z:\Users\Yong\Paper_Tunde_FisherInformation\GlobalTuningYongTrial.dat'];
outfile = ['Z:\Data\Tempo\Batch Files\Tunde\VIP_CELLS_AIHUA\dat_output\batch_VIP_VestiSigHor_3D.dat'];
% outfile = ['Z:\Data\Tempo\Batch Files\Tunde\VIP_CELLS_AIHUA\dat_output\batch_VIP_VisSigHor_3D.dat'];
% outfile = ['Z:\Data\Tempo\Batch Files\Tunde\VIP_CELLS_AIHUA\dat_output\test.dat'];
%outfile = ['Z:\Data\Tempo\Batch Files\Tunde\dat_files\Yong_3D_horizontalplane.dat'];
%outfile = ['C:\Aihua\z_TempOutputs\VesOnly_Aihua.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%output_specifications 
%---------------------------------------------------------------------------------------
return;