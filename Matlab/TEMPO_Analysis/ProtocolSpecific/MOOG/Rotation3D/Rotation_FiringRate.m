function Rotation_FiringRate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
eye_data = data.eye_data(:,:,~null_trials & select_trials);
spike_data = data.spike_data(:,:,~null_trials & select_trials);

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

select1 = logical( (azimuth==0)  & (elevation==90 ) & (stim_type ==1) );
select2 = logical( (azimuth==0)  & (elevation==-90 ) & (stim_type ==1) );
% for i=1:length(spike_rates)    
%     temp = find(spike_data(2,:,i)==1);
%     temp2 = temp(temp>1000 ); 
%     syncpulse_index(i) = temp2(1); % align with first syncpulse   
%     if syncpulse_index(i)<1200 %sometimes it is strange that syncpulse does not happen in a proper time
%         a6(i,:) = eye_data(5,round(syncpulse_index(i)/5) : round(syncpulse_index(i)/5)+399,i); % spike data has 5 times higher resolution than eye data
%     else
%        a6(i,:) = eye_data(5,201 : 600, i);  
%     end
% end
% acc1(1,:) = median(a6(select1,:));
% acc2(1,:) = median(a6(select2,:));
% acc= mean([-acc1,acc2]);
% acc = acc-mean(acc(1:50));
% figure;
% plot(acc,'b-');


%% ADD CODE HERE FOR PLOTTING
resp_mat = [];
for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=1: length(unique_condition)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition==unique_condition(k)) );
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(spike_rates(select));
                resp_std(k, j, i) = std(spike_rates(select));
            else
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_std(k, j, i) = resp_std(k,j,1);
            end
        end        
    end
end

ves(1:26) = [resp_mat(1,1,1),squeeze(resp_mat(1,2,:))',squeeze(resp_mat(1,3,:))',squeeze(resp_mat(1,4,:))',resp_mat(1,5,1)];
% vis(1:26) = [resp_mat(2,1,1),squeeze(resp_mat(2,2,:))',squeeze(resp_mat(2,3,:))',squeeze(resp_mat(2,4,:))',resp_mat(2,5,1)];
% com(1:26) = [resp_mat(3,1,1),squeeze(resp_mat(3,2,:))',squeeze(resp_mat(3,3,:))',squeeze(resp_mat(3,4,:))',resp_mat(3,5,1)];
ves_std(1:26) = [resp_std(1,1,1),squeeze(resp_std(1,2,:))',squeeze(resp_std(1,3,:))',squeeze(resp_std(1,4,:))',resp_std(1,5,1)];
% vis_std(1:26) = [resp_std(2,1,1),squeeze(resp_std(2,2,:))',squeeze(resp_std(2,3,:))',squeeze(resp_std(2,4,:))',resp_std(2,5,1)];
% com_std(1:26) = [resp_std(3,1,1),squeeze(resp_std(3,2,:))',squeeze(resp_std(3,3,:))',squeeze(resp_std(3,4,:))',resp_std(3,5,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 500 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
% buff = sprintf(sprint_txt, FILE, ves, ves_std );
buff = sprintf(sprint_txt, FILE, ves );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\analysis\VesOnly_syed.dat'];
% outfile = ['C:\Aihua\z_TempOutputs\VesOnly_Aihua.dat'];
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
output_specifications 
%---------------------------------------------------------------------------------------
return;