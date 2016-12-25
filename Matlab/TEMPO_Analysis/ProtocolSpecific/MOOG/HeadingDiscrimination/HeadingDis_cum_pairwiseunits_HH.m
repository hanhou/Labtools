% % isolate more than 2 units from single unit recording data by offline
% spikesorting, analyze clustering structure and noise correlation among
% units --YG, 03/08
% %-----------------------------------------------------------------------------------------------------------------------
function HeadingDis_cum_pairwiseunits_HH(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_conflict_angle = data.moog_params(CONFLICT_ANGLE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_heading = data.moog_params(HEADING,:,MOOG);

%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_azimuth);
if sum(temp_conflict_angle)>=0 % exist such variable
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (temp_conflict_angle==0) ); 
else          % for other protocols, there is no such a variable
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
motion_coherence = temp_motion_coherence(~null_trials & select_trials);
conflict_angle = temp_conflict_angle(~null_trials & select_trials);
num_sigmas = temp_num_sigmas(~null_trials & select_trials);
heading = temp_heading(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_motion_coherence = munique(motion_coherence');
unique_conflict_angle = munique(conflict_angle');
unique_num_sigmas = munique(num_sigmas');
unique_heading = munique(heading');

if sum(temp_conflict_angle)>=0 % exist such variable
    Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial | temp_conflict_angle~=0 );
else
    Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial );
end
temp_spike_data = data.spike_data(:, :);
for i = 1 : length(Discard_trials)
    temp_spike_data( :, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 99;
end
StartEventBin(1)=996;

% extract channel information
channelnum_temp = size(temp_spike_data);
channelnum = channelnum_temp(1,1); % how many channels
channelcount = 0;
%SpikeChan = 1; % define the first channel you want to start here
for c = 1: channelnum
    temp(1,:) = temp_spike_data( c, find(temp_spike_data(1,:)~=99) );
    spikesum(c) = sum(temp(1,:));
    if c>=SpikeChan & spikesum(c)>20 & c~=2 % all those channels later than the first channel, but exclude the second synpulse channel
        channelcount = channelcount+1;
        channel_analyze(channelcount) = c; % the final channels that need to analyze
    end    
end

for c = 1 : channelcount  
    w=1;
    % count spikes
    spike_data(1,:) = temp_spike_data( channel_analyze(c), find(temp_spike_data(1,:)~=99) );
    spike_data(1, find(spike_data>10) ) = 1; % something is absolutely wrong  

    for ss =  1 : length(azimuth) % ss marks the index of trial
%            spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+offset+5000*(ss-1) : StartEventBin(1)+offset+duration+5000*(ss-1)) ) ; 
              %spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+500+5000*(ss-1) : StartEventBin(1)+500+115+1000+5000*(ss-1)) ) ; 
              spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+60+5000*(ss-1) : StartEventBin(1)+60+1000+5000*(ss-1)) ) ;   % HH20130906
              disp('Hardcoded: spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+60+5000*(ss-1) : StartEventBin(1)+60+1000+5000*(ss-1)) )');
              keyboard;
%              spike_rates2(ss) = sum( spike_data(1,StartEventBin(1)+1501+115+5000*(ss-1) : StartEventBin(1)+2000+115+5000*(ss-1)) ) ;
%              spike_rates(ss) = spike_rates1(ss)+spike_rates2(ss);
%           spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+(w-1)*50+5000*(ss-1) : StartEventBin(1)+115+(w-1)*50+500+5000*(ss-1)) ) ; 
%        spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+501+115+5000*(ss-1) : StartEventBin(1)+501+115+1000+5000*(ss-1)) ) ; 
    end
    spike_rates_channel(c,:) = spike_rates;  
    spike_rates_z = spike_rates; % no high-pass filter

    % creat basic matrix represents each response vector
  %  resp = [];
    for cc = 1:length(unique_motion_coherence)
        for k=1:length(unique_stim_type)
            for i=1:length(unique_heading)
                select = logical( (heading==unique_heading(i))  & (stim_type==unique_stim_type(k)) & (motion_coherence==unique_motion_coherence(cc)) );
                resp{c,k}(cc,i) = mean(spike_rates(select));                
                % z-score data for spike count correlation analysis
                z_dist = spike_rates_z(select);
                if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
                   z_dist = (z_dist - mean(z_dist))/std(z_dist);
                else
                    z_dist = 0;
                end
                Z_Spikes(select) = z_dist;  
            end 
            if cc>1 & k==1
               temp = polyfit(unique_heading, resp{c,k}(1,:)',1);
               var(c,k) = resp{c,k}(1,floor((1+length(unique_heading))/2));
            else
                temp = polyfit(unique_heading, resp{c,k}(cc,:)',1);  
                var(c,k) = resp{c,k}(cc,floor((1+length(unique_heading))/2));
            end
            slope(c,k) = temp(1);
        end
    end
    Z_Spikes_channel(c,:) = Z_Spikes;
end

% now analyze noise correlation between pairs
% first compute all possible pairs 
channelcount_temp = channelcount-1;
channelcount_loop = 0;
while channelcount_temp>=1
    channelcount_loop = channelcount_loop + channelcount_temp;
    channelcount_temp = channelcount_temp - 1;
end

outloop = 1;
insideloop = outloop+1;
for i = 1:channelcount_loop % all possible pairs
    % remove slow fluctuations at every 20 trials
    ztemp = floor(length(spike_rates)/20);
    zz1 = Z_Spikes_channel(outloop,:);
    zz2 = Z_Spikes_channel(insideloop,:);
    z1all=[];
    z2all=[];
    for zz=1:ztemp            
        if zz<ztemp
            z1=(zz1(1+(zz-1)*20:20+(zz-1)*20)-mean(zz1(1+(zz-1)*20:20+(zz-1)*20)))/std(zz1(1+(zz-1)*20:20+(zz-1)*20));
            z2=(zz2(1+(zz-1)*20:20+(zz-1)*20)-mean(zz2(1+(zz-1)*20:20+(zz-1)*20)))/std(zz2(1+(zz-1)*20:20+(zz-1)*20));
        else
            z1=(zz1(1+(zz-1)*20:end)-mean(zz1(1+(zz-1)*20:end)))/std(zz1(1+(zz-1)*20:end));
            z2=(zz2(1+(zz-1)*20:end)-mean(zz2(1+(zz-1)*20:end)))/std(zz2(1+(zz-1)*20:end));
        end
        z1all=[z1all z1];
        z2all=[z2all z2];
        z1all(z1all>3)=3;  % cutoff between -3 and 3
        z1all(z1all<-3)=-3; 
        z2all(z2all<-3)=-3;  % cutoff between -3 and 3
        z2all(z2all>3)=3;
    end
    % noise correlation with all stimuli conditions included
    [rr,pp] = corrcoef(z1all,z2all);  
    noise_r(i) = rr(1,2)
    noise_p(i) = pp(1,2)

    % separated between stimuli conditions
    for cc = 1:length(unique_motion_coherence)
        for k=1:length(unique_stim_type) % ananlyze noise correlation in different conditions, if find no difference, combine later
     %       select_stim = logical( stim_type==unique_stim_type(k) & (motion_coherence==unique_motion_coherence(cc)) );
            select_stim = logical( stim_type==unique_stim_type(k) );
            % noise correlation with stim type separated
            Z_Spikes1 = z1all(select_stim);
            Z_Spikes2 = z2all(select_stim);

            [rr,pp] = corrcoef(Z_Spikes1,Z_Spikes2);     
    %        noise_r_stim{i}(cc,k)= rr(1,2);
            noise_r_stim{i}(k)= rr(1,2);
        end
    end

    if insideloop < channelcount % more to run
       insideloop = insideloop+1;
    else
       outloop = outloop+1;    
       insideloop = outloop+1;
    end 
end
noise_r;
noise_r_stim{1};
figure(2);
for i = 1:channelcount
    subplot(channelcount,1,i);
    plot(spike_rates_channel(i,:), 'o-');
    xlim([1 length(spike_rates)]);
end

% HH20130905
if ishandle(3) close(3); end
set(figure(3),'color','w'); 
plot(Z_Spikes1,Z_Spikes2,'o'); hold on;
plot(Z_Spikes1(heading>0),Z_Spikes2(heading>0),'ro');
plot(Z_Spikes1(heading==0),Z_Spikes2(heading==0),'go');

set(gca,{'xtickmode' 'ytickmode'},{'auto' 'auto'});
title(['NC = ' num2str(noise_r(1)) ', p = ' num2str(noise_p(1))]);
axis square;
min_Z = min([Z_Spikes1 Z_Spikes2]);
max_Z = max([Z_Spikes1 Z_Spikes2]);
axis([min_Z max_Z min_Z max_Z]);
hold on; line([min_Z max_Z], [min_Z max_Z]);

%% output data
sprint_txt = ['%s\t'];
for i = 1 : 5000
     sprint_txt = [sprint_txt, ' %1.3f\t'];    
end

%for c = 1 : channelcount_loop  
% for c = 1 : channelcount  
%     if length(unique_stim_type) >=2
% %        buff = sprintf(sprint_txt,raster_max{c}(1,:),raster_max{c}(2,:));         
%  %        buff = sprintf(sprint_txt,FILE,noise_r(c),noise_p(c),noise_r_stim(c,:),noise_p_stim(c,:),corrcoef_r_unit(c,:),corrcoef_p_unit(c,:),corrcoef_r_stim,corrcoef_p_stim); 
% %         buff = sprintf(sprint_txt,FILE,corrcoef_r_stim,corrcoef_p_stim);
% %         buff = sprintf(sprint_txt,FILE,noise_rr(c,:), noise_rr_stim{c}(:,1),noise_rr_stim{c}(:,2));  
%         buff = sprintf(sprint_txt,FILE,slope(c,:), var(c,:));
%     else
%    %      buff = sprintf(sprint_txt,FILE,noise_r(c),noise_p(c),noise_r_stim(c,:),noise_p_stim(c,:),corrcoef_r_unit(c,:),corrcoef_p_unit(c,:));        
%  %       buff = sprintf(sprint_txt,FILE,noise_pp(c,1),noise_rr(c,:), noise_r_stim(c,:)); 
% %        buff = sprintf(sprint_txt,raster_max{c}(1,:));
%         buff = sprintf(sprint_txt,FILE,slope(c,:), var(c,:));
%     end
%  %   outfile = ['Z:\Users\Yong\noisecorr_1000msnew.dat'];
%  %   outfile = ['Z:\Users\Yong\noisecorr_psth.dat'];
% %    outfile = ['Z:\Users\Yong\noisecorr_dynamic200msslide50msvesvis.dat'];    
% %     outfile = ['C:\Documents and Settings\yong\My Documents\work\noisecorr_discrimination.dat'];
%  %   outfile = ['Z:\Users\Yong\noisecorr_dynamic500msnew.dat'];
%     outfile = ['Z:\Users\Yong\noisecorr_discrimstimtype.dat'];
% %   outfile = ['Z:\Users\Yong\noisecorr_congruency.dat'];
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t');
%         fprintf(fid, '\r\n');
%     end
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);  
% end

return;

