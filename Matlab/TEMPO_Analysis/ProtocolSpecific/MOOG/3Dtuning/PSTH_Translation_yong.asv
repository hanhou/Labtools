%----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	Yong, 04/07/09
%-----------------------------------------------------------------------------------------------------------------------

function PSTH_Translation_yong(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
% SpikeChan = 3;
% plfp_chan = 2;

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_spike_data = data.spike_data(SpikeChan,:);
temp_spike_rates = data.spike_rates(SpikeChan, :);  

% if abs(sum(sum(sum(data.plfp_data)))) >0
%    temp_plfp_data(1,1:5000,1:length(temp_azimuth)) = data.plfp_data(plfp_chan,1:2:10000,:); % sample 1k Hz
% else
%    temp_plfp_data(1,1:5000,1:length(temp_azimuth)) = 0;
% end

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_trials ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
% stim_duration = length(temp_spike_data)/length(temp_azimuth);
% spike_data = data.spike_data(1, ((BegTrial-1)*stim_duration+1):EndTrial*stim_duration);
spike_rates= temp_spike_rates(~null_trials & select_trials);
% plfp_data(:,:) = squeeze(temp_plfp_data(1,:,~null_trials & select_trials));
% plfp_data_null(:,:) = squeeze(temp_plfp_data(1,:,null_trials));
% plfp_null(1,:) = median(plfp_data_null(:,:),2);
spike_data_null(:,:) = squeeze(data.spike_data(SpikeChan,:,null_trials));
spike_null = mean(spike_data_null,2);
% notice that this bad_trials is the number without spon trials 

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition_num = munique(condition_num');

% add parameters here
% timebin for plot PSTH
timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);
for i=1:100    
    spike_null_psth(1,i) = sum(spike_null(1+(i-1)*50:50+(i-1)*50));
end

% find spontaneous trials which azimuth,elevation,stim_type=-9999
spon_found = find(null_trials==1);     

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 99;
end
spike_data(1,:) = temp_spike_data( 1, find(temp_spike_data(1,:)~=99) );
spike_data(1, find(spike_data>10) ) = 1; % something is absolutely wrong 

% count spikes from raster data (spike_data)
time_step=1;
for k=1: length(unique_condition_num)
    count = 0;
    for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)
            temp_count = []; % initialize 
            temp_raster = [];             
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );            
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                count=count+1;
                resp{k}(j,i) = mean(spike_rates(select));
                act_found = find( select==1 );
                % count spikes per timebin on every same condition trials
                for repeat=1:length(act_found) 
                    for n=1:(x_length)
                        temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                        time_step=time_step+timebin;
                    end
                    time_step=1; 
                    temp_raster(repeat,:) = spike_data(1, frequency*(act_found(repeat)-1)+1 : frequency*(act_found(repeat)) ); % raw raster, i.e. bin=1ms
                end
                count_y_trial{i,j,k}(:,:) = temp_count;  % each trial's PSTH 
                count_y_trial_raster{i,j,k}(:,:) = temp_raster; % each trial's raster
                % get the average of the total same conditions if repetion is > 1
           %     if (length(act_found) > 1);
                dim=size(temp_count);
                if dim(1) > 1;
                   count_y{i,j,k} = mean(temp_count);
                else
                   count_y{i,j,k}= temp_count;     % for only one repetition cases
                end   
                max_count(k,count) = max(count_y{i,j,k});
%                 plfp{i,j,k}=median(plfp_data(:,select),2);
%                 max_count_plfp(k,count) = max(plfp{i,j,k});
%                 min_count_plfp(k,count) = min(plfp{i,j,k});
             else
                resp{k}(j,i) = 0; 
                count_y{i,j,k}=count_y{1,j,k};
                count_y_trial{i,j,k}(:,:) = count_y_trial{1,j,k}(:,:);
                count_y_trial_raster{i,j,k}(:,:) = count_y_trial_raster{1,j,k}(:,:);
     %           plfp{i,j,k} = plfp{1,j,k};
             end   
             % normalize count_y
             if max(count_y{i,j,k})~=0;
                count_y_norm{i,j,k}=count_y{i,j,k} / max(count_y{i,j,k});
             else
                count_y_norm{i,j,k}=0;
             end
        end
    end  
    % now find the peak
    [row_max, col_max] = find( resp{k}(:,:)==max(max(resp{k}(:,:))) );
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_m{k}=row_max(1);
    col_m{k}=col_max(1);
    if max(count_y{col_max(1), row_max(1), k})~=0;
     %  count_y_max{k} = count_y{col_max(1), row_max(1), k} / max(count_y{col_max(1), row_max(1), k}); % normalized
       count_y_max{k} = count_y{col_max(1), row_max(1), k};
      %count_y_max{k} = count_y{4, 3, k};
    else
       count_y_max{k} =0;
    end

    unique_azimuth(col_m{k})
    unique_elevation(row_m{k})
end

% transform to plot format
for k=1: length(unique_condition_num)
    pc = 1;
    for j=1:length(unique_elevation)
        for i = 1:length(unique_azimuth)+1  
            if (i < length(unique_azimuth))                                 
                count_y_trial_raster_transform{i,j,k}(:,:)=count_y_trial_raster{length(unique_azimuth)-i,j,k}(:,:);
                count_y_transform{i,j,k}=count_y{length(unique_azimuth)-i,j,k};   
%                plfp_transform{i,j,k}=plfp{length(unique_azimuth)-i,j,k};  
            elseif(i==length(unique_azimuth))
                count_y_trial_raster_transform{i,j,k}(:,:)=count_y_trial_raster{i,j,k}(:,:); 
                count_y_transform{i,j,k}=count_y{i,j,k}(:,:);  
 %               plfp_transform{i,j,k}=plfp{i,j,k}(:,:);
            else
                count_y_trial_raster_transform{i,j,k}(:,:)=count_y_trial_raster{length(unique_azimuth)-1,j,k}(:,:); 
                count_y_transform{i,j,k}=count_y{length(unique_azimuth)-1,j,k}(:,:); 
%                plfp_transform{i,j,k}=plfp{length(unique_azimuth)-1,j,k}(:,:); 
            end 
            count_y_transform26{k}(:,pc) = count_y_transform{i,j,k}(1,22:62);
            pc = pc+1;
        end
    end
    if col_m{k}<=7
        col_m_transform{k}=8 - col_m{k};
    else 
        col_m_transform{k}=8;
    end
end

% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [(StartEventBin(1,1)+115)/timebin, (StartEventBin(1,1)+115)/timebin];
x_stop =  [(StopEventBin(1,1)+115)/timebin,  (StopEventBin(1,1)+115)/timebin];
% define figure
% now plot
xtime_plfp = 1:0.02:100.99; % 10000 data points for LFP
stimulus_duration = StopEventBin(1) - StartEventBin(1); % maybe 2 seconds, 1 second or others...
for k=1: length(unique_condition_num)
    figure(k+1);
    set(k+1,'Position', [5,5 1000,700], 'Name', 'PSTH');
    orient landscape;
    title(FILE);
    axis off;
    xoffset=0;
    yoffset=0;
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.42;
        xoffset = 0;
    end
    % output some text 
    axes('position',[0 0 1 0.9]); 
    xlim([-50,50]);
    ylim([-50,50]);
    text(-30+xoffset*100,52+yoffset*110, h_title{k} );
    text(-30+xoffset*100+10,52+yoffset*110, num2str(SpikeChan) );
    %text(-47,-40, 'Azim: 270       225       180        135        90        45        0        315        270');  
    temp=fliplr(unique_azimuth');unique_azimuth_plot=[temp(2:end) temp(1:2)];clear temp
    text(-47,-40, ['Azim:' num2str(unique_azimuth_plot)]);  
    text(25,-40, 'Translation');
    axis off;
    hold on;
    for i=1:length(unique_azimuth)+1                    % aizmuth 270 are plotted two times in order to make circular data
        for j=1:length(unique_elevation)
            if i==1 | (i>1 & j>1 &j<5)
                axes('position',[0.11*i-0.09+xoffset (0.92-0.15*j)+yoffset 0.1 0.12]); 
                bar( x_time,count_y_transform{i,j,k}(1,:));  
%                [AX,H1,H2]=plotyy( x_time,count_y_transform{i,j,k}(1,:), xtime_plfp, plfp_transform{i,j,k},'plot');  
                hold on;
            %    plot(xtime_plfp, plfp_transform{i,j,k}, 'r-');
                plot( x_start, [0,max(max_count(k,:))], 'r-');
                plot( x_stop,  [0,max(max_count(k,:))], 'r-');
                set( gca, 'xticklabel', ' ' );
                % set the same scale for all plot
%                 xlim(AX(1),[10,20+20*(StopEventBin(1,1)-StartEventBin(1,1))/1000+10]); % from 0 to 500ms after stimulus offset
%                 xlim(AX(2),[10,20+20*(StopEventBin(1,1)-StartEventBin(1,1))/1000+10]);
%                 ylim(AX(1),[0,max(max_count(k,:))]);
%                 ylim(AX(2),[min(min_count_plfp(k,:)),max(max_count_plfp(k,:))]);
            elseif i==5 & j==5
                axes('position',[0.1*i-0.05+xoffset (0.92-0.15*j)+yoffset 0.08 0.12]); 
%                [AX,H1,H2]=plotyy( x_time,spike_null_psth, xtime_plfp, plfp_null,'plot');  
                hold on;
            %    plot(xtime_plfp, plfp_transform{i,j,k}, 'r-');
                plot( x_start, [0,max(max_count(k,:))], 'r-');
                plot( x_stop,  [0,max(max_count(k,:))], 'r-');
                set( gca, 'xticklabel', ' ' );
                % set the same scale for all plot
%                 xlim(AX(1),[10,20+20*(StopEventBin(1,1)-StartEventBin(1,1))/1000+10]); % from 0 to 500ms after stimulus offset
%                 xlim(AX(2),[10,20+20*(StopEventBin(1,1)-StartEventBin(1,1))/1000+10]);
%                 ylim(AX(1),[0,max(max_count(k,:))]);
%                 ylim(AX(2),[min(min_count_plfp(k,:)),max(max_count_plfp(k,:))]);
            end
        end    
    end 
%    xoffset=xoffset+0.46;    
end

aa = 0;
figure(5);
set(5,'Position', [5,5 1000,700], 'Name', 'PSTH');
orient landscape;
title(FILE);
axis off;
% now plot
% col_m_transform{3}=3
% row_m{3}=3
for k=1: length(unique_condition_num) 
    for r = 1:length( count_y_trial_raster_transform{col_m_transform{k},row_m{k},k}(:,1) ) 
        % raster for each trial
        axes('position',[0.3*(k-1)+0.05 0.9-r*0.03 0.25 0.015]);                                 
        bar( count_y_trial_raster_transform{col_m_transform{k},row_m{k},k}(r,:) );  
        hold on;
        plot( [StartEventBin(1) StartEventBin(1)], [0 1], 'r-');
        plot( [StopEventBin(1) StopEventBin(1)], [0 1], 'r-');
        axis off;
        set( gca, 'xticklabel', ' ' );
        set( gca, 'yticklabel', ' ' );
        xlim([1,StopEventBin(1)+1000]);
        ylim([0,1]);
    end              
    axes('position',[0.3*(k-1)+0.05 0.05 0.25 0.25]);  
    bar( x_time,count_y_transform{col_m_transform{k},row_m{k},k}(1,:)); 
    hold on;
%     plot( x_start, y_marker, 'r-');
%     plot( x_stop,  y_marker, 'r-');
    set( gca, 'xticklabel', ' ' );
    % set the same scale for all plot
    xlim([0,x_length*(stimulus_duration+1000+1000)/frequency]);
    ylim([0,max(max_count(k,:))]);
end
% % Brief frequency analysis for Tanya; 9-18-09
% for i=1:length(unique_azimuth)
%     for j=1:length(unique_elevation)
%         [f, amp, resp_phase] = FT(x_time(1:62)*0.05, count_y{i,j,k}(1,1:62), length(x_time(1:62)), 1, 0);
%         max_freq(i,j) = f(find(amp==max(amp)));
%         max_amp(i,j) = max(amp);
%     end
% end
% save(['C:\moog\' FILE '_maxfreq.mat']);
% data_x = x_time*0.05; data_y=count_y{i,j,k}(1,:); FFT_PTS=length(x_time); DC_remove = 1; plot_flag = 1;

%---------------------------------------------------------------------------------------
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s'];
for i = 1 : x_length * 3
     sprint_txt = [sprint_txt, ' %1.2f'];    
end
%buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_45(1,:)/8, count_y_45(2,:)/8,count_y_90(1,:)/8, count_y_90(2,:)/8,count_y_135(1,:)/8, count_y_135(2,:)/8, count_y_180(1,:), count_y_180(2,:));  
%buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2} );  % for 2 conditions
if length(unique_stim_type)==1
    buff = sprintf(sprint_txt, FILE, count_y_max{1} );
else
    buff = sprintf(sprint_txt, FILE, count_y_max{1}, count_y_max{2} );
end
% for 1 conditions
%buff = sprintf(sprint_txt, FILE, count_trial_beg,count_trial_end ); 

%outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\DirectionTuning3D_PSTH_Tanya.dat'];
outfile = ['Z:\Users\Yong\FEF\psth.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
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

return;

