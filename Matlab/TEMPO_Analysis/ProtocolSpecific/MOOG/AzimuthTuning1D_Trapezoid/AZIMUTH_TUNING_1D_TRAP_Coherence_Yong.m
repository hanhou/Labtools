% Coherence analysis --YG, 07/30/2010
% %-----------------------------------------------------------------------------------------------------------------------
function AZIMUTH_TUNING_1D_TRAP_Coherence_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_spike_data = data.spike_data(:,:); 
temp_lfp_data10000 = data.plfp_data(:,:); 
temp_lfp_data = data.plfp_data(:,1:2:end);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 

azimuth = temp_azimuth(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';

SpikeChan(1)=1;
SpikeChan(2)=3;
% construct lfp signal
% strobe=dlmread('m21c12r1strobe.txt'); % strobe signal
% load ADtemp;
% load ADtimestamp;
% findindex=find(strobe(:,2)==12);
% lfp_data1 =[];
% lfp_data2 =[];
% for i=1:length(findindex)    
%     timestamp=strobe(findindex(i)-2,1); % time stamp for the stimulus onset (code 04)
%     lfp_index=find(ADtimestamp<=timestamp);
%     lfp_data1 = [lfp_data1 ADtemp(lfp_index(end)-1000+1:lfp_index(end)+4000,1)'];
%     lfp_data2 = lfp_data1;
%  %   lfp_data2 = [lfp_data2 ADtemp(lfp_index(end)-1000+1:lfp_index(end)+4000,2)'];
% end
% temp_lfp_data(1,:) = lfp_data1;
% temp_lfp_data(2,:) = lfp_data2;
% find spontaneous trials which azimuth,elevation,stim_type=-9999
spon_found = find(null_trials==1);     
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( :, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 99;
end

for k=1:length(unique_stim_type)
    for i=1:length(unique_azimuth)
        select = logical( azimuth==unique_azimuth(i) & stim_type==unique_stim_type(k) );
        repetition(i,k) = length(find(select==1));
    end
end  
repetition_min = min(min(repetition));

for c=1:length(SpikeChan)    
    % spikes
    spike_data(1,:) = temp_spike_data( SpikeChan(c), find(temp_spike_data(1,:)~=99) );
    spike_data(1, find(spike_data>1) ) = 1; % something is absolutely wrong  
    % lfp    
    lfpp_data(1,:) = temp_lfp_data(c,find(temp_spike_data(1,:)~=99) );
    
    for k=1:length(unique_stim_type)
        for i=1:length(unique_azimuth)
            selectfind = find( azimuth==unique_azimuth(i) & stim_type==unique_stim_type(k) );
            for j=1: repetition_min
                resp_spike{c,i,k}(:,j) = single(spike_data(1,1+(selectfind(j)-1)*5000:2500+(selectfind(j)-1)*5000)); % only use 2.5 secs
                resp_lfp{c,i,k}(:,j) = lfpp_data(1,1+(selectfind(j)-1)*5000:2500+(selectfind(j)-1)*5000); % only use 2.5 secs
            end
        end
    end    
    repetition_min = min(min(repetition));
end

% now run coherence between spike and lfp
params.tapers=[10 19];
params.Fs=1000;
params.fpass=[0 100];
params.pad=0;
params.err=[2 0.05];
params.trialave=1;
movingwin=[0.25,0.01]; % 250ms window sliding with 10ms

% define figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', '1D Direction Tuning trapezoid');
orient landscape;
title(FILE);
set(0, 'DefaultAxesXTickMode', 'auto', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'manual');
axis off;

for k=1: length(unique_stim_type) 
    for i=1:length(unique_azimuth) 
        i
        %[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpb(lfp_signal,spike_signal,params);
        [C,phi,S12,S1,S2,t,f]=cohgramcpb(resp_lfp{1,i,k},resp_spike{2,i,k},movingwin,params); % lfp-spike against time
    %    [C,t,f]=mtspecgramc(resp_lfp{1,i,k},movingwin,params); % lfp spectrum against time
        %[C{k}(:,i),f]=mtspectrumc(resp_lfp{1,i,k}(1021:2000),params); % lfp spectrum 
        
  %      axes('position',[0.05+0.12*(i-1) 0.65-0.3*(k-1) 0.1 0.15]);                          
%         contourf(t,f,C');        
%         xlim([t(1) t(end)]);
%         ylim([f(1),f(end)]); 
      %  caxis([0 0.2]);
    end     
end

return;

