% % isolate more than 2 units from single unit recording data by offline
% spikesorting, analyze clustering structure and noise correlation among
% units --YG, 03/08
% %-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_1D_pairwiseunits(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
    temp_spike_data( :, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 99;
end
StartEventBin(1)=996;

repetition = floor( length(azimuth) / (length(unique_azimuth)*length(unique_stim_type)) ); % take minimum repetition
repetition
% extract channel information
channelnum_temp = size(temp_spike_data);
channelnum = channelnum_temp(1,1); % how many channels
channelcount = 0;

%SpikeChan = 5; % define the first channel you want to start here
for c = 1: channelnum
    temp(1,:) = temp_spike_data( c, find(temp_spike_data(1,:)~=99) );
    spikesum(c) = sum(temp(1,:));
    if c>=SpikeChan & spikesum(c)>20 & c~=2 % all those channels later than the first channel, but exclude the second synpulse channel
        channelcount = channelcount+1;
        channel_analyze(channelcount) = c; % the final channels that need to analyze
    end    
end

for w=1:1    
    for c = 1 : channelcount         
        % count spikes
        spike_data(1,:) = temp_spike_data( channel_analyze(c), find(temp_spike_data(1,:)~=99) );
        spike_data(1, find(spike_data>10) ) = 1; % something is absolutely wrong  

        for ss =  1 : length(azimuth) % ss marks the index of trial
%            spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+offset+5000*(ss-1) : StartEventBin(1)+offset+duration+5000*(ss-1)) ) ; 

%              spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+5000*(ss-1) : StartEventBin(1)+500+115+5000*(ss-1)) ) ; 
%              spike_rates2(ss) = sum( spike_data(1,StartEventBin(1)+1501+115+5000*(ss-1) : StartEventBin(1)+2000+115+5000*(ss-1)) ) ;
 %           spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+(w-1)*50+5000*(ss-1) : StartEventBin(1)+115+(w-1)*50+500+5000*(ss-1)) ) ; 
           %spike_rates(ss) = sum(spike_data(1,StartEventBin(1)+501+115+5000*(ss-1) : StartEventBin(1)+501+115+1000+5000*(ss-1)) ) ; 
           spike_rates(ss) = sum(spike_data(1,StartEventBin(1)+60+5000*(ss-1) : StartEventBin(1)+60+1000+5000*(ss-1)) ) ;   % HH20130906
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
                repeats_temp(k,i) = length(select_find);
                
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
                resp_mean{c,k}(i,w) = mean(spike_rates(select));        
                resp_variance{c,k}(i,w) = (std(spike_rates(select)))^2;
                resp_sse(i,k) = sum( (spike_rates(select)-mean(spike_rates(select))).^2 );                
                
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
            ff{c}(w,k) = mean( resp_std(:,k).^2 ./ resp_f );
            ffvariance{c,k}(:,w) = resp_std(:,k).^2;
            ffmean{c,k}(:,w) = resp_f;
            resp_trial{k}(:, 1) = resp_trial_temp;
            resp_trial{k}(:, 2) = resp_trial_group;
            raster_maxx = find(raster_sum==max(raster_sum));
            raster_max_temp = raster_mean(raster_maxx(1),:);  
            DDI(c,k) = (max(resp(:,k))-min(resp(:,k)))/( max(resp(:,k))-min(resp(:,k))+2*sqrt( sum(resp_sse(:,k))/(sum(repeats_temp(k,:))-length(unique_azimuth)) ) );
            for j = 1:40                 
                raster_max{c}(k,j) = sum( raster_max_temp(1,1+(j-1)*50:50+(j-1)*50) ); 
   %             raster_max{c}(k,j) = sum( raster_max_temp(1,1+(j-1)*50:500+(j-1)*50) );
            end              
        end
        repeats = min(min(repeats_temp(:,:)));
        Z_Spikes_channel(c,:) = Z_Spikes;
        
        % vectorsum and calculate preferred direction
        % vectors must be symetric, otherwise there will be a bias both on
        % preferred direction and the length of resultant vector
        % the following is to get rid off non-symetric data, hard code temporally
        if length(unique_azimuth) >8
            findtemp=find(unique_azimuth==0); resp_s(1,:) = resp(findtemp,:);
            findtemp=find(unique_azimuth==45); resp_s(2,:) = resp(findtemp,:);
            findtemp=find(unique_azimuth==90); resp_s(3,:) = resp(findtemp,:);
            findtemp=find(unique_azimuth==135); resp_s(4,:) = resp(findtemp,:);
            findtemp=find(unique_azimuth==180); resp_s(5,:) = resp(findtemp,:);
            findtemp=find(unique_azimuth==225); resp_s(6,:) = resp(findtemp,:);
            findtemp=find(unique_azimuth==270); resp_s(7,:) = resp(findtemp,:);
            findtemp=find(unique_azimuth==315); resp_s(8,:) = resp(findtemp,:);
        else
            resp_s(:,:) = resp(:,:);
        end
        unique_azimuth_s(1:8) = [0,45,90,135,180,225,270,315];
        unique_elevation_s(1:8) = 0;  
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
        end
        % method 2 to remove slow fluctuation
        for zz=1:length(spike_rates)
            startbin = zz-10;
            endbin = zz+10;
            if startbin<1
                startbin=1;
            end
            if endbin>length(spike_rates)
                endbin=length(spike_rates);
            end
            z1all(zz) = Z_Spikes_channel(outloop,zz)-mean(Z_Spikes_channel(outloop,startbin:endbin));
            z2all(zz) = Z_Spikes_channel(insideloop,zz)-mean(Z_Spikes_channel(insideloop,startbin:endbin));
        end        
        
        % don't remove slow fluctuation
%         z1all = zz1;
%         z2all = zz2;

        z1all(z1all>3)=3;  % cutoff between -3 and 3
        z1all(z1all<-3)=-3;
        z2all(z2all<-3)=-3;  % cutoff between -3 and 3
        z2all(z2all>3)=3;
        % noise correlation with all stimuli conditions included
        [rr,pp] = corrcoef(z1all,z2all);  
        noise_r(i) = rr(1,2);
        noise_p(i) = pp(1,2);
        noise_rr(i,w) = rr(1,2); % for output
        noise_pp(i,w) = pp(1,2);
        
        % separated between stimuli conditions
        for k=1:length(unique_stim_type) % ananlyze noise correlation in different conditions, if find no difference, combine later
            Z_Spikes1 = [];
            Z_Spikes2 = [];            
            select_stim = logical( stim_type==unique_stim_type(k) );
            % noise correlation with stim type separated
            Z_Spikes1{i}(k,:) = z1all(select_stim);
            Z_Spikes2{i}(k,:) = z2all(select_stim);

            [rr,pp] = corrcoef(Z_Spikes1{i}(k,:),Z_Spikes2{i}(k,:)); 
            noise_r_stim(i,k) = rr(1,2);
            noise_rr_stim{i}(w,k)= rr(1,2);
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
    
    disp('NC       NC_p       SC        SC_p');
    for i = 1:length(noise_r)
        fprintf('  %g\t %g\t  %g\t %g\n',noise_r(i), noise_p(i), corrcoef_r_unit(i), corrcoef_p_unit(i));
    end
    
%     % examine stability of firing rate over time
    figure(20);
    for i = 1:channelcount
        subplot(channelcount,1,i);
        plot(spike_rates_channel(i,:), 'o-');
        xlim([1 length(spike_rates)]);
    end

%     Define figure
    outloop = 1;
    insideloop = outloop+1;
    for i = 1:channelcount_loop % all possible pairs
        set(figure(i+30),'color','w');
        set(i+30,'Position', [5,15 700,600], 'Name', '1D Direction Tuning');
        orient landscape;
        set(0, 'DefaultAxesXTickMode', 'auto', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'auto');
        
        f{1,1}='bo-'; f{2,1}='g+--'; 
%         f{1,2}='ro-'; f{2,2}='r+--'; 
%         f{1,3}='go-'; f{2,3}='g+--'; 

        for k=1: length(unique_stim_type)
            subplot(2,2,k)        
            set(errorbar(unique_azimuth, resp_pair{outloop}(:,k), resp_err_pair{outloop}(:,k), f{1,k} ),'linewidth',2); % the 1st unit
            hold on;
            set(errorbar(unique_azimuth, resp_pair{insideloop}(:,k), resp_err_pair{insideloop}(:,k), f{2,k} ),'linewidth',2); % the 2nd unit
            ylim_ = ylim;
            line([az(outloop,k) az(outloop,k)],ylim_,'color',f{1,k}(1),'linewidth',2);
            line([az(insideloop,k) az(insideloop,k)],ylim_,'color',f{2,k}(1),'linewidth',2);
            
            if k==1
               ylabel('spikes/s');
            end
            xlabel('azimuth');
            xlim( [min(unique_azimuth) max(unique_azimuth)] );
            set(gca, 'xtick',[unique_azimuth]);
            title( [FILE '   ' num2str(unique_stim_type(k)) ', ' num2str(p_1D(outloop,k)) ', ' num2str(p_1D(insideloop,k)) ', SC = ' num2str(corrcoef_r_unit(i,k)) ', p = ' num2str(corrcoef_p_unit(i))] );                 
            set(gca,'xtick',0:90:360);
            
            % noise correlation
            subplot(2,2,2)
            
            colors = colormap(lines);
            % HH20130905 For different stimulus angle
            for azis = 1 : length(unique_azimuth)
                select_stim = logical( stim_type==unique_stim_type(k) & azimuth == unique_azimuth(azis));
                % plot(Z_Spikes_channel(outloop,select_stim),Z_Spikes_channel(insideloop,select_stim), 'o','color',colors(azis,:));
                plot(Z_Spikes_channel(outloop,select_stim),Z_Spikes_channel(insideloop,select_stim), 'ko','markerfacecolor','k','markersize',5);
                hold on;
                [r_azis,p_azis]=corrcoef(Z_Spikes_channel(outloop,select_stim),Z_Spikes_channel(insideloop,select_stim));
                r_diff_azi(azis) = r_azis(1,2);
                p_diff_azi(azis) = p_azis(1,2);
            end
            
            fprintf('%s\t %s\t %s\n',num2str(unique_azimuth'),num2str(r_diff_azi),num2str(p_diff_azi));

            para = polyfit(Z_Spikes_channel(outloop,:),Z_Spikes_channel(insideloop,:),1);
            xmin = min([Z_Spikes_channel(outloop,:),Z_Spikes_channel(insideloop,:)]);
            xmax = max([Z_Spikes_channel(outloop,:),Z_Spikes_channel(insideloop,:)]);
            xlim([xmin, xmax]);
            ylim([xmin, xmax]);
            axis square;
            xlabel('Z-scored response (Cell 1)');
            ylabel('Z-scored response (Cell 2)');
            hold on;
            
            plot([xmin,xmax],polyval(para,[xmin,xmax]),'k--');
            %plot([xmin,xmax], [xmin, xmax],'k--');
             
            title(['NC = ' num2str(noise_r_stim(i,k)) ', p = ' num2str(noise_p_stim(i,k))]);
            
            % HH20130905
            subplot(2,2,3); box on;
            plot(unique_azimuth,r_diff_azi','o-'); hold on;
            plot(unique_azimuth(p_diff_azi<0.05),r_diff_azi(p_diff_azi<0.05),'o','markerfacecolor','b');
            xlabel('azimuth');
            ylabel('NC');
            xlim( [min(unique_azimuth) max(unique_azimuth)] );
            ylim_ = ylim;
            line([az(outloop,k) az(outloop,k)],ylim_,'color',f{1,k}(1));
            line([az(insideloop,k) az(insideloop,k)],ylim_,'color',f{2,k}(1));
            title(['NC = ' num2str(noise_r_stim(i,k)) ]);

            
        end
        
        if insideloop < channelcount % more to run
            insideloop = insideloop+1;
        else
           outloop = outloop+1;    
        end 
    end

    % %SU vs. MU------------------------------------------------------------------------
    % %cross correlogram to show whether the SU has been removed cleanly
    % CorrSUsubtractBefore  = data.corr{1};
    % CorrSUsubtractAfter  = data.corr{2};
    % %------------------------------------------------------------------------

    %% ---------------------------------------------------------------------------------------
    % Also, write out some summary data to a cumulative summary file
end
sprint_txt = ['%s\t'];
for i = 1 : 5000
     sprint_txt = [sprint_txt, ' %1.3f\t'];    
end

% for c = 1 : channelcount_loop  
%  % for c = 1 : channelcount
% %        buff = sprintf(sprint_txt,FILE, p_1D(c,1), p_1D(c,2), raster_max{c}(1,:),raster_max{c}(2,:)); 
% %        buff = sprintf(sprint_txt,FILE, ff{c}(:,1), ff{c}(:,2) );
%   %       buff = sprintf(sprint_txt,FILE,noise_r(:),noise_p(c),noise_r_stim(c,:),noise_p_stim(c,:),corrcoef_r_unit(c,:),corrcoef_p_unit(c,:),corrcoef_r_stim,corrcoef_p_stim);        
% %         outfile = ['Z:\Users\Yong\noisecorr_1000msnew.dat'];
% %        buff = sprintf(sprint_txt,FILE,noise_rr(c,:), noise_rr_stim{c}(:,1),noise_rr_stim{c}(:,2));  
%   %       buff = sprintf(sprint_txt,FILE,p_1D(c,1),p_1D(c,2),resp_pair_horizontalplane{c}(:,1),resp_pair_horizontalplane{c}(:,2)); 
%  %        buff = sprintf(sprint_txt,FILE,repeats); 
%  % buff = sprintf(sprint_txt,FILE, raster_max{c}(1,:), raster_max{c}(2,:));  
%   %  outfile = ['Z:\Users\Sheng\noiseves.dat'];
%     buff = sprintf(sprint_txt,FILE,noise_r); 
%     outfile = ['Z:\Users\Yong\noise.dat'];
%     
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

