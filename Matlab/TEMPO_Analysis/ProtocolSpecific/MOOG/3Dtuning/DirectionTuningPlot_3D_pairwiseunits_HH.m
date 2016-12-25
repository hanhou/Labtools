function DirectionTuningPlot_3D_pairwiseunits_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_duration = data.moog_params(DURATION,:,MOOG);

%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_azimuth);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
duration = temp_duration(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_duration = munique(duration');

Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
temp_spike_data = data.spike_data(:, :);
for i = 1 : length(Discard_trials)
    temp_spike_data( :, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 99;
end
StartEventBin(1)=996;

repetition = floor( length(azimuth) / (26*length(unique_stim_type)) );

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
            %             spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+5000*(ss-1) : StartEventBin(1)+500+115+5000*(ss-1)) ) ;
            %              spike_rates2(ss) = sum( spike_data(1,StartEventBin(1)+1501+115+5000*(ss-1) : StartEventBin(1)+2000+115+5000*(ss-1)) ) ;
            %           spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+(w-1)*50+5000*(ss-1) : StartEventBin(1)+115+(w-1)*50+500+5000*(ss-1)) ) ;
            %            spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+501+115+5000*(ss-1) : StartEventBin(1)+501+115+1000+5000*(ss-1)) ) ;
            if unique_duration == 2000
                %                 % use the middle 1 second
                %       %          spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+5000*(ss-1) : StartEventBin(1)+615+5000*(ss-1)) ) ;
                spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+615+5000*(ss-1) : StartEventBin(1)+1615+5000*(ss-1)) ) ;
            elseif unique_duration == 1000
                %                 % use the whole 1 second
                spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+115+5000*(ss-1) : StartEventBin(1)+1115+5000*(ss-1)) ) ;
            end
        end
        spike_rates_channel(c,:) = spike_rates;
        
        % remove slow drift of responsiveness over time by a high pass filter
        %         temp_spike_rates_z = FIR_Filter(spike_rates, 20, 100, 'high', 20, 0);
        %         spike_rates_z = temp_spike_rates_z(1:end);
        spike_rates_z = spike_rates; % no high-pass filter
        
        % create basic matrix represents each response vector
        resp = [];
        for k=1:length(unique_stim_type)
            pc = 0;
            resp_trial_temp =[];
            resp_trial_group = [];
            resp_trial_temp_horizontalplane =[];
            resp_trial_group_horizontalplane = [];
            for j=1:length(unique_elevation)
                for i=1:length(unique_azimuth)
                    select = logical( (azimuth==unique_azimuth(i))  & (elevation==unique_elevation(j))  & (stim_type==unique_stim_type(k)) );
                    select_find = find( (azimuth==unique_azimuth(i))  & (elevation==unique_elevation(j))  & (stim_type==unique_stim_type(k)) );
                    select_horizontalplane = logical( (azimuth==unique_azimuth(i))  & (elevation==0)  & (stim_type==unique_stim_type(k)) );
                    if (sum(select) > 0)
                        pc = pc+1;
                        repeats_temp(k,pc) = length(select);
                        spike_temp = spike_rates(select);
                        raster = [];
                        for jj=1:length(spike_temp)
                            raster(jj,:) = spike_data(1,StartEventBin(1)+115+5000*(select_find(jj)-1) : StartEventBin(1)+115+2000+5000*(select_find(jj)-1));
                        end
                        resp_trial_temp = [resp_trial_temp, spike_temp];
                        resp_trial_group_temp =[];
                        resp_trial_group_temp(1:length(spike_temp)) = pc;
                        resp_trial_group = [resp_trial_group,resp_trial_group_temp];
                        % for pvalue in the horizontal plane
                        if unique_elevation(j)==0
                            spike_temp_horizontalplane = spike_rates(select_horizontalplane);
                            resp_trial_temp_horizontalplane = [resp_trial_temp_horizontalplane, spike_temp_horizontalplane];
                            resp_trial_group_temp_horizontalplane =[];
                            resp_trial_group_temp_horizontalplane(1:length(spike_temp_horizontalplane)) = i;
                            resp_trial_group_horizontalplane = [resp_trial_group_horizontalplane,resp_trial_group_temp_horizontalplane];
                        end
                        
                        resp(k,j,i) = mean(spike_rates(select));
                        resp_std(k,j,i) = std(spike_rates(select));
                        resp_std_horizontal(i,k) = std(spike_rates(select_horizontalplane));
                        resp_s(i,k) = mean(spike_rates(select_horizontalplane));
                        resp_26(pc,k) = mean(spike_rates(select));
                        resp_std_26(pc,k) = std(spike_rates(select));
                        resp_mean{c,k}(pc,w) = mean(spike_rates(select));
                        resp_variance{c,k}(pc,w) = (std(spike_rates(select)))^2;
                        resp_sse(pc,k) = sum( (spike_rates(select)-mean(spike_rates(select))).^2 );
                        
                        raster_mean(pc,:) = mean(raster(:,:)); % mean raster
                        raster_sum(pc) = sum(raster_mean(pc,500:1500));
                        % z-score data for spike count correlation analysis
                        z_dist = spike_rates_z(select);
                        if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
                            z_dist = (z_dist - mean(z_dist))/std(z_dist);
                        else
                            z_dist = 0;
                        end
                        Z_Spikes(select) = z_dist;
                    else
                        resp(k,j,i) = resp(k,j,1);         % for later vectorsum use
                        resp_std(k,j,i) = 0;
                    end
                end
            end
            resp_stdmean{c}(k) = mean(resp_std_horizontal(:,k));
            %    resp_f = resp_s(:,k);
            resp_f = resp_26(:,k);
            resp_f(find(resp_f<=1)) = 1;
            ff{c}(w,k) = mean( resp_std_26(:,k).^2 ./ resp_f );
            ffvariance{c,k}(:,w) = resp_std_26(:,k).^2;
            ffmean{c,k}(:,w) = resp_f;
            resp_pair{c}(1,k) = resp(k,1,1);
            resp_pair{c}(2:9,k) = resp(k,2,:);
            resp_pair{c}(10:17,k) = resp(k,3,:);
            resp_pair{c}(18:25,k) = resp(k,4,:);
            resp_pair{c}(26,k) = resp(k,5,1);
            resp_pair_horizontalplane{c}(:,k) = resp_s(:,k);
            resp_pair_ste_horizontalplane{c}(:,k) = resp_std_horizontal(:,k)/sqrt(repetition);
            raster_maxx = find(raster_sum==max(raster_sum));
            raster_max_temp = raster_mean(raster_maxx(1),:);
            for j = 1:40
                raster_max{c}(k,j) = sum( raster_max_temp(1,1+(j-1)*50:50+(j-1)*50) );
                %             raster_max{c}(k,j) = sum( raster_max_temp(1,1+(j-1)*50:500+(j-1)*50) );
            end
            resp_trial{k}(:, 1) = resp_trial_temp;
            resp_trial{k}(:, 2) = resp_trial_group;
            resp_trial_horizontalplane{k}(:, 1) = resp_trial_temp_horizontalplane;
            resp_trial_horizontalplane{k}(:, 2) = resp_trial_group_horizontalplane;
            
            N{c,k} = squeeze(resp(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually
            [azi(c,k), ele(c,k), amp] = vectorsum(N{c,k});
            % anova1
            p_1D(c,k) = anovan(resp_trial{k}(:,1),{resp_trial{k}(:,2)},'display','off');
            p_1D_horizontalplane(c,k) = anovan(resp_trial_horizontalplane{k}(:,1),{resp_trial_horizontalplane{k}(:,2)},'display','off');
        end
        Z_Spikes_channel(c,:) = Z_Spikes;
        repeats = min(min(repeats_temp(:,:)));
        
        % congruency between stim type
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
            select_stim = logical( stim_type==unique_stim_type(k) );
            Z_Spikes1 = [];
            Z_Spikes2 = [];
            % noise correlation with stim type separated
            Z_Spikes1 = z1all(select_stim);
            Z_Spikes2 = z2all(select_stim);
            zzz1{i,k}=Z_Spikes1;
            zzz2{i,k}=Z_Spikes2;
            
            [rr,pp] = corrcoef(Z_Spikes1,Z_Spikes2);
            noise_r_stim(i,k) = rr(1,2)
            noise_rr_stim{i}(w,k)= rr(1,2);
            noise_p_stim(i,k) = pp(1,2)
            
            % this is only the regular correlation between two tuning curves
            [rr,pp] = corrcoef(resp_pair{outloop}(:,k),resp_pair{insideloop}(:,k));
            corrcoef_r_unit(i,k) = rr(1,2)
            corrcoef_p_unit(i,k) = pp(1,2)
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
    
    %     % Define figure
    azi_cos = [1,2,3,4,5,6,7,8];
    ele_sin = [-1,-0.707,0,0.707,1];
    outloop = 1;
    insideloop = outloop+1;
    for i = 1:channelcount_loop % all possible pairs
        figure(i+1);
        set(i+1,'Position', [5,15 750,650], 'Name', '3D Direction Tuning');
        orient landscape;
        set(0, 'DefaultAxesXTickMode', 'auto', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'auto');
        
        for k=1: length(unique_stim_type)
            subplot(3,3,k)
            contourf( azi_cos, ele_sin, squeeze(N{outloop,k}) );
            color bar;
            if k==1
                ylabel('elevation');
            end
            xlabel('azimuth');
            set(gca, 'ydir' , 'reverse');
            title( [FILE num2str(unique_stim_type(k)) '     ' num2str(p_1D(outloop,k)) ] );
            
            subplot(3,3,k+3)
            contourf( azi_cos, ele_sin, squeeze(N{insideloop,k}) );
            color bar;
            if k==1
                ylabel('elevation');
            end
            xlabel('azimuth');
            set(gca, 'xtick',[unique_azimuth]);
            set(gca, 'ydir' , 'reverse');
            title( [FILE num2str(unique_stim_type(k)) '     ' num2str(p_1D(insideloop,k)) ] );
            
            % noise correlation
            subplot(3,3,k+6)
            select_stim = logical( stim_type==unique_stim_type(k) );
            plot(Z_Spikes_channel(outloop,select_stim),Z_Spikes_channel(insideloop,select_stim), 'o');
            xmin = min([Z_Spikes_channel(outloop,select_stim),Z_Spikes_channel(insideloop,select_stim)]);
            xmax = max([Z_Spikes_channel(outloop,select_stim),Z_Spikes_channel(insideloop,select_stim)]);
            xlim([xmin, xmax]);
            ylim([xmin, xmax]);
            hold on;
            plot([xmin,xmax], [xmin, xmax],'-');
            title( num2str([noise_r_stim(i,k) noise_p_stim(i,k)]) );
        end
        
        if insideloop < channelcount % more to run
            insideloop = insideloop+1;
        else
            outloop = outloop+1;
        end
    end
end
sprint_txt = ['%s\t'];
for i = 1 : 5000
    sprint_txt = [sprint_txt, ' %1.3f\t'];
end




% for c = 1 : channelcount_loop
% %  for c = 1 : channelcount
% %        buff = sprintf(sprint_txt,FILE, p_1D(c,1), p_1D(c,2), raster_max{c}(1,:),raster_max{c}(2,:));
% %        buff = sprintf(sprint_txt,FILE, ff{c}(:,1), ff{c}(:,2) );
%   %       buff = sprintf(sprint_txt,FILE,noise_r(:),noise_p(c),noise_r_stim(c,:),noise_p_stim(c,:),corrcoef_r_unit(c,:),corrcoef_p_unit(c,:),corrcoef_r_stim,corrcoef_p_stim);
% %         outfile = ['Z:\Users\Yong\noisecorr_1000msnew.dat'];
% %        buff = sprintf(sprint_txt,FILE,noise_rr(c,:), noise_rr_stim{c}(:,1),noise_rr_stim{c}(:,2));
%   %       buff = sprintf(sprint_txt,FILE,p_1D(c,1),p_1D(c,2),resp_pair_horizontalplane{c}(:,1),resp_pair_horizontalplane{c}(:,2));
%  %        buff = sprintf(sprint_txt,FILE,repeats);
% % buff = sprintf(sprint_txt,FILE, raster_max{c}(1,:), raster_max{c}(2,:));
%   %  outfile = ['Z:\Users\Sheng\noiseves.dat'];
%     buff = sprintf(sprint_txt,FILE,noise_r);
%     outfile = ['Z:\Users\Yong\noise.dat'];
% %    outfile = ['Z:\Users\Yong\tuningsimilaritytwounit.dat'];
% %    outfile = ['Z:\Users\Yong\noisecorr_dynamic200msslide50msvesvis.dat'];
% %     buff = sprintf(sprint_txt,FILE,noise_rr(c,:));
%  %    outfile = ['Z:\Users\Yong\noisecorr_baselinecombined.dat'];
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

%     outloop = 1;
%     insideloop = outloop+1;
%     for i = 1:channelcount_loop % all possible pairs
%         buff = sprintf(sprint_txt,FILE,length(spike_rates),unique_stim_type,noise_r(i),noise_p(i),noise_r_stim(i,:),noise_p_stim(i,:),corrcoef_r_unit(i,:),corrcoef_p_unit(i,:),corrcoef_r_stim(outloop), corrcoef_r_stim(insideloop));
%     %    buff = sprintf(sprint_txt,FILE,repetition, length(spike_rates));
%     %buff = sprintf(sprint_txt, FILE, spon_resp(1),DDI_all(1,:),DDI_forward(1,:),HTI_modified(1,:),az(1,:),p_1D(1,:) );
%        buff = sprintf(sprint_txt,FILE,noise_rr(i,:),noise_pp(i,:));
% %         if w==1
% %            outfile = ['Z:\Users\Yong\noisecorr2000ms.dat'];
% %         elseif w==2
% %             outfile = ['Z:\Users\Yong\noisecorr1500ms.dat'];
% %         elseif w==3
% %             outfile = ['Z:\Users\Yong\noisecorr1000ms.dat'];
% %         elseif w==4
% %             outfile = ['Z:\Users\Yong\noisecorr500ms.dat'];
% %         elseif w==5
% %             outfile = ['Z:\Users\Yong\noisecorr200ms.dat'];
% %         else
% %             outfile = ['Z:\Users\Yong\noisecorr100ms.dat'];
% %         end
%         outfile = ['Z:\Users\Yong\noisecorr_shift200ms_extra.dat'];
%   %      outfile = ['Z:\Users\Yong\noisecorr_windowlength100ms.dat'];
%         printflag = 0;
%         if (exist(outfile, 'file') == 0)    %file does not yet exist
%             printflag = 1;
%         end
%         fid = fopen(outfile, 'a');
%         if (printflag)
%             fprintf(fid, 'FILE\t');
%             fprintf(fid, '\r\n');
%         end
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%         fclose(fid);
%
%         if insideloop < channelcount % more to run
%            insideloop = insideloop+1;
%         else
%            outloop = outloop+1;
%            insideloop = outloop+1;
%         end
%     end
%end

return;