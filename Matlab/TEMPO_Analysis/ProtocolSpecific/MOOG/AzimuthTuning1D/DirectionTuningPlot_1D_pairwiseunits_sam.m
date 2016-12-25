% isolate more than 2 units from single unit recording data by offline
% spikesorting, analyze spike timing correlation among different
% units --SAM, 05/08
% %-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_1D_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE)
%--------------------------------------------------------------------------
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_interocular_dist = data.moog_params(INTEROCULAR_DIST,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);

%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_azimuth);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
motion_coherence = temp_motion_coherence(~null_trials & select_trials);
interocular_dist = temp_interocular_dist(~null_trials & select_trials);
num_sigmas = temp_num_sigmas(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_motion_coherence = munique(motion_coherence');
unique_interocular_dist = munique(interocular_dist');
unique_num_sigmas = munique(num_sigmas');

Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
temp_spike_data = data.spike_data(:, :);
for i = 1 : length(Discard_trials)
    temp_spike_data( :, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 9999;
end
StartEventBin(1)=996;

repetition = floor( length(azimuth) / (length(unique_azimuth)*length(unique_stim_type)) ); % take minimum repetition
repetition
% extract channel information
channelnum_temp = size(temp_spike_data);
channelnum = channelnum_temp(1,1); % how many channels
channelcount = 0;
%SpikeChan = 6; % define the first channel you want to start here
for c = 1: channelnum
    temp(1,:) = temp_spike_data( c, find(temp_spike_data(1,:)~=9999) );
    spikesum(c) = sum(temp(1,:));
    if c>=SpikeChan & spikesum(c)>20 & c~=2 % all those channels later than the first channel, but exclude the second synpulse channel
        channelcount = channelcount+1;
        channel_analyze(channelcount) = c; % the final channels that need to analyze
    end    
end
% channelcount = 3;
% channel_analyze = [4 5 6];

% for w=1:31    
    for c = 1 : channelcount         
        % count spikes
        spike_data(1,:) = temp_spike_data( channel_analyze(c), find(temp_spike_data(1,:)~=9999) );
        spike_data(1, find(spike_data>10) ) = 1; % something is absolutely wrong  

        for ss =  1 : length(azimuth) % ss marks the index of trial
%            spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+offset+5000*(ss-1) : StartEventBin(1)+offset+duration+5000*(ss-1)) ) ; 
%              spike_rates1(ss) = sum( spike_data(1,StartEventBin(1)+115+5000*(ss-1) : StartEventBin(1)+500+115+5000*(ss-1)) ) ; 
%              spike_rates2(ss) = sum( spike_data(1,StartEventBin(1)+1501+115+5000*(ss-1) : StartEventBin(1)+2000+115+5000*(ss-1)) ) ;
%             spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+(w-1)*50+5000*(ss-1) : StartEventBin(1)+115+(w-1)*50+500+5000*(ss-1)) ) ; 
           spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+501+115+5000*(ss-1) : StartEventBin(1)+501+115+1000+5000*(ss-1)) ) ; 
        end
        spike_rates_channel(c,:) = spike_rates;
          
        %now remove the slow drift responsiveness effet, comment out this part
        %if you don't want to, this is for z-scored only     
%         temp_spike_rates_z = FIR_Filter(spike_rates, 20, 100, 'high', 20, 0);
%         spike_rates_z = temp_spike_rates_z(1:end);
        spike_rates_z = spike_rates; % no high-pass filter

        % creat basic matrix represents each response vector
        resp = [];
        for k=1:length(unique_stim_type)   
            resp_trial_temp =[];
            resp_trial_group = [];
            for i=1:length(unique_azimuth)
                select = logical( (azimuth==unique_azimuth(i))  & (stim_type==unique_stim_type(k)) );
                select_find = find( (azimuth==unique_azimuth(i))  & (stim_type==unique_stim_type(k)) );
                
                spike_temp = spike_rates(select);    
                raster = [];
                for jj=1:length(spike_temp)
                    raster(jj,:) = spike_data(1,StartEventBin(1)+115+5000*(select_find(jj)-1) : StartEventBin(1)+115+2000+5000*(select_find(jj)-1));
                end
                resp_trial_temp = [resp_trial_temp, spike_temp];
                resp_trial_group_temp =[];
                resp_trial_group_temp(1:length(spike_temp)) = i;
                resp_trial_group = [resp_trial_group,resp_trial_group_temp];  
                
                resp(i, k) = mean(spike_rates(select));        
                resp_std(i,k) = std(spike_rates(select));                
                resp_err(i,k) = std(spike_rates(select)) / sqrt(repetition); 
                raster_mean(i,:) = mean(raster(:,:)); % mean raster
                raster_sum(i) = sum(raster_mean(i,500:1500));                

                % z-score data for spike count correlation analysis
                z_dist = spike_rates_z(select);
                if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
                   z_dist = (z_dist - mean(z_dist))/std(z_dist);
                else
                    z_dist = 0;
                end
                Z_Spikes(select) = z_dist;  
            end  
            resp_stdmean{c}(k) = mean(resp_std(:,k));
            resp_f = resp(:,k);
            resp_f(find(resp_f<=1)) = 1;
%             ff{c}(w,k) = mean( resp_std(:,k).^2 ./ resp_f );
            resp_trial{k}(:, 1) = resp_trial_temp;
            resp_trial{k}(:, 2) = resp_trial_group;
            raster_maxx = find(raster_sum==max(raster_sum));
            raster_max_temp = raster_mean(raster_maxx(1),:);
            for j = 1:40                 
                raster_max{c}(k,j) = sum( raster_max_temp(1,1+(j-1)*50:50+(j-1)*50) ); 
            end
        end
        Z_Spikes_channel(c,:) = Z_Spikes;

        % vectorsum and calculate preferred direction
        % vectors must be symetric, otherwise there will be a bias both on
        % preferred direction and the length of resultant vector
        % the following is to get rid off non-symetric data, hard code temporally
        if length(unique_azimuth) >8
            resp_s(1,:) = resp(1,:);
            resp_s(2,:) = resp(2,:);
            resp_s(3,:) = resp(4,:);
            resp_s(4,:) = resp(6,:);
            resp_s(5,:) = resp(7,:);
            resp_s(6,:) = resp(8,:);
            resp_s(7,:) = resp(9,:);
            resp_s(8,:) = resp(10,:);
            unique_azimuth_s(1:8) = [0,45,90,135,180,225,270,315];
            unique_elevation_s(1:8) = 0;
        else
            resp_s(:,:) = resp(:,:);
            unique_azimuth_s = unique_azimuth;
            unique_elevation_s(1:length(unique_azimuth)) = 0;
        end
        resp_pair{c}(:,:) = resp(:,:);
        resp_err_pair{c}(:,:) = resp_err(:,:);
        resp_pair_horizontalplane{c}(:,:) = resp_s(:,:);

        % preferred direction and p value
        for k = 1: length(unique_stim_type)
            [az(c,k), el(c,k), amp(c,k)] = vectorsumAngle(resp_s(:,k), unique_azimuth_s, unique_elevation_s);
            p_1D(c,k) = anovan(resp_trial{k}(:,1),{resp_trial{k}(:,2)},'display','off');        
        end 

        % congruency, data need to be at least >= 2 stim type
        if length(unique_stim_type)>=2 
            [rr,pp] = corrcoef(resp_pair{c}(:,1),resp_pair{c}(:,2)); % temporarily hard coded
            corrcoef_r_stim(c) = rr(1,2);
            corrcoef_p_stim(c) = pp(1,2);
        else
            corrcoef_r_stim(c) = 99;
            corrcoef_p_stim(c) = 99;
        end
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
            z2all(z2all<-3)=-3;  % cutoff between -3 and 3
        end
        % noise correlation with all stimuli conditions included
        [rr,pp] = corrcoef(z1all,z2all);  
        noise_r(i) = rr(1,2);
        noise_p(i) = pp(1,2);
%         noise_rr(i,w) = rr(1,2); % for output
%         noise_pp(i,w) = pp(1,2);
        
        % separated between stimuli conditions
        for k=1:length(unique_stim_type) % ananlyze noise correlation in different conditions, if find no difference, combine later
            Z_Spikes1 = [];
            Z_Spikes2 = [];
            select_stim = logical( stim_type==unique_stim_type(k) );
            % noise correlation with stim type separated
            Z_Spikes1(k,:) = z1all(select_stim);
            Z_Spikes2(k,:) = z2all(select_stim);

            [rr,pp] = corrcoef(Z_Spikes1(k,:),Z_Spikes2(k,:)); 
            noise_r_stim(i,k) = rr(1,2);
%             noise_rr_stim{i}(w,k)= rr(1,2);
            noise_p_stim(i,k) = pp(1,2);             

            % this is only the regular correlation between two tuning curves
            [rr,pp] = corrcoef(resp_pair{outloop}(:,k),resp_pair{insideloop}(:,k));
            corrcoef_r_unit(i,k) = rr(1,2);
            corrcoef_p_unit(i,k) = pp(1,2);
            if length(unique_stim_type)>=2
                [rr,pp] = corrcoef([resp_pair{outloop}(:,1);resp_pair{outloop}(:,2)],[resp_pair{insideloop}(:,1);resp_pair{insideloop}(:,2)]);
                corrcoef_r_twounit(i) = rr(1,2);
                corrcoef_p_twounit(i) = pp(1,2);
            else
                corrcoef_r_twounit(i) = corrcoef_r_unit(i,k);
                corrcoef_p_twounit(i) = corrcoef_p_unit(i,k);
            end
        end
        
        if insideloop < channelcount % more to run
           insideloop = insideloop+1;
        else
           outloop = outloop+1;    
           insideloop = outloop+1;
        end 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %show file name and some values in text
%     axes('position',[0.05,0.85, 0.9,0.1] );
%     xlim( [0,100] );
%     ylim( [0,length(unique_stim_type)*(channelcount-1)] );
%     text(0, length(unique_stim_type)*(channelcount-1), FILE);
%     text(15,length(unique_stim_type)*(channelcount-1),'preferred           p                  noisecorrelation');
%     count = 0;
%     for k=1:length(unique_stim_type)
%         for c=1:2
%             count = count+1;
%             text(0,length(unique_stim_type)*(channelcount-1)-count, num2str(unique_stim_type(k)));
%             text(15,length(unique_stim_type)*(channelcount-1)-count, num2str(az(c,k)) );
%             text(25,length(unique_stim_type)*(channelcount-1)-count, num2str(p_1D(c,k)) );
%             text(35,length(unique_stim_type)*(channelcount-1)-count, num2str(noise_r_stim(k)) );
%         end
%     end
%     axis off;
    %% ---------------------------------------------------------------------------------------
    % Also, write out some summary data to a cumulative summary file
% end
sprint_txt = ['%s\t'];
for i = 1 : 5000
     sprint_txt = [sprint_txt, ' %1.3f\t'];    
end

%for c = 1 : channelcount_loop  
for c = 1 : channelcount  
    if length(unique_stim_type) >=2
        buff = sprintf(sprint_txt,FILE, channel_analyze(c), az(c,1), az(c,2) , p_1D(c,1), p_1D(c,2), corrcoef_r_stim(c),corrcoef_p_stim(c) ); 
    else
        buff = sprintf(sprint_txt,FILE, channel_analyze(c), az(c,1), 99 , p_1D(c,1), 99, corrcoef_r_stim(c),corrcoef_p_stim(c) );
    end
    outfile = ['Z:\Users\Sam\Results\tuning.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t channel\t VesPref\t VisPref\t VesP\t VisP\t congruency_r\t congruency_p');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);  
end

return;