function DirectionTuningPlot_3D_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
%---------------------------------------------------------------------------
% data get and chop
% Path_Defs;
% ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
% 
%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_spike_data_5 = data.spike_data(5,:);
temp_spike_data_6 = data.spike_data(6,:);

%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_azimuth);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');

% find spontaneous trials which azimuth,elevation,stim_type=-9999
spon_found = find(null_trials==1);

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data_5)/length(temp_azimuth);%find length of stimulus duration by dividing the total spike data vector by the number of trials
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)%mark trials to discard
    temp_spike_data_5( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
    temp_spike_data_6( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data(1,:) = temp_spike_data_5( temp_spike_data_5~=9999 );%transfer non marked to new vector
spike_data(1, find(spike_data(1,:)>100) ) = 1; % something is absolutely wrong
spike_data(2,:) = temp_spike_data_6( temp_spike_data_6~=9999 );
spike_data(2, find(spike_data(2,:)>100) ) = 1; % something is absolutely wrong

% extract channel information
% channelnum_temp = size(data.spike_data);
% channelnum = channelnum_temp(1,1); % how many channels
% trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-2*(length(unique_azimuth)-1)) * length(unique_stim_type) + 1;
% repetitions = floor( length(find(select_trials==1)) / trials_per_rep );

%for non 90 elevations
for a=1:length(unique_stim_type)
    for b=2:length(unique_elevation)-1
        for c=1:length(unique_azimuth)
            select = find( (stim_type==unique_stim_type(a))  & (elevation==unique_elevation(b))  & (azimuth==unique_azimuth(c)) );
            for d=1:length(select)
                spikes_chan5_less90{a,b-1,c}(d,:) = spike_data(1, (select(d)-1)*5000+1 : select(d)*5000);
                spikes_chan6_less90{a,b-1,c}(d,:) = spike_data(2, (select(d)-1)*5000+1 : select(d)*5000);
                spikes_chan5_less90_2s{a,b-1,c}(d,:) = spikes_chan5_less90{a,b-1,c}(d,996:3006);
                spikes_chan6_less90_2s{a,b-1,c}(d,:) = spikes_chan6_less90{a,b-1,c}(d,996:3006);
            end
        end
    end
end
% for a=1:length(unique_stim_type)
%     for b=1:length(unique_elevation)-2
%         for c=1:length(unique_azimuth)
%             if isempty(spikes_chan5_less90_2s{a,b,c})==1
%                 spikes_chan5_less90_2s{a,b,c}=zeros(2,2011);
%             end
%             if isempty(spikes_chan6_less90_2s{a,b,c})==1
%                 spikes_chan6_less90_2s{a,b,c}=zeros(2,2011);
%             end
%         end
%     end
% end
%for 90 elevations
for a=1:length(unique_stim_type)
    for b=1:length(unique_elevation)
        if unique_elevation(b)==-90 || unique_elevation(b)==90
            for c=1:length(unique_azimuth)
                select = find( (stim_type==unique_stim_type(a))  & (elevation==unique_elevation(b))  & (azimuth==unique_azimuth(c)) );
                for d=1:length(select)
                    spikes_chan5_90{a,b}(d,:) = spike_data(1, (select(d)-1)*5000+1 : select(d)*5000);
                    spikes_chan6_90{a,b}(d,:) = spike_data(2, (select(d)-1)*5000+1 : select(d)*5000);
                    spikes_chan5_90_2s_temp{a,b}(d,:) = spikes_chan5_90{a,b,c}(d,996:3006);
                    spikes_chan6_90_2s_temp{a,b}(d,:) = spikes_chan6_90{a,b,c}(d,996:3006);
                end
            end
        end
    end
end
for a=1:length(unique_stim_type)
    spikes_chan5_90_2s_temp{a,2}= spikes_chan5_90_2s_temp{a,end};
    spikes_chan6_90_2s_temp{a,2}= spikes_chan6_90_2s_temp{a,end};
end
for a=1:length(unique_stim_type)
    for b=1:2
        spikes_chan5_90_2s{a,b}=spikes_chan5_90_2s_temp{a,b};
        spikes_chan6_90_2s{a,b}=spikes_chan6_90_2s_temp{a,b};
    end
end
% for a=1:length(unique_stim_type)
%     for b=1:2
%         if isempty(spikes_chan5_90_2s{a,b})==1
%             spikes_chan5_90_2s{a,b}=zeros(2,2011);
%         end
%         if isempty(spikes_chan6_90_2s{a,b})==1
%             spikes_chan6_90_2s{a,b}=zeros(2,2011);
%         end
%     end
% end
%--------------------------------------------------------------------------
% % find raw xcorr and shift across every cell
for a=1:length(unique_stim_type)%iterate across stim type
    for b=1:length(unique_elevation)-2%iterate across elevations less the 90s
        for c=1:length(unique_azimuth)%iterate across azimuth directions
            if isempty(spikes_chan5_less90_2s{a,b,c})~=1
                for d = 1 : length(spikes_chan5_less90_2s{a,b,c}(:,1))%iterate across trials
                    if d<length(spikes_chan5_less90_2s{a,b,c}(:,1))%calculate xcorr across channel within cells(stimulus conditions) and correct with next trial correlation
                        if sum(spikes_chan5_less90_2s{a,b,c}(d,:))>0 && sum(spikes_chan6_less90_2s{a,b,c}(d,:))>0 && sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))>0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))*sum(spikes_chan6_less90_2s{a,b,c}(d,:))/4)^(1/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d+1,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))*sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))/4)^(1/2);
                        elseif sum(spikes_chan5_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') ;
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d+1,:),50,'unbiased') ;
                        elseif sum(spikes_chan6_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))/2)^(1/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d+1,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))/2)^(1/2);
                        elseif sum(spikes_chan5_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(d,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased');
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d+1,:),50,'unbiased') / (sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))/2)^(1/2);
                        elseif sum(spikes_chan5_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan6_less90_2s{a,b,c}(d,:))/2)^(1/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d+1,:),50,'unbiased') ;
                        elseif sum(spikes_chan5_less90_2s{a,b,c}(d,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan6_less90_2s{a,b,c}(d,:))/2)^(1/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d+1,:),50,'unbiased') / (sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))/2)^(1/2);
                        elseif sum(spikes_chan6_less90_2s{a,b,c}(d,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))/2)^(1/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d+1,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))*sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))/4)^(1/2);
                        elseif sum(spikes_chan6_less90_2s{a,b,c}(d+1,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))*sum(spikes_chan6_less90_2s{a,b,c}(d,:))/4)^(1/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d+1,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))/2)^(1/2);
                        end
                    elseif d==length(spikes_chan5_less90_2s{a,b,c}(:,1))%boundary condition for last trial in cell, correlation correctetd against first trial in opposite channel
                        if sum(spikes_chan5_less90_2s{a,b,c}(d,:))>0 && sum(spikes_chan6_less90_2s{a,b,c}(d,:))>0 && sum(spikes_chan6_less90_2s{a,b,c}(1,:))>0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))*sum(spikes_chan6_less90_2s{a,b,c}(d,:))/4)^(1/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(1,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))*sum(spikes_chan6_less90_2s{a,b,c}(1,:))/4)^(1/2);
                        elseif sum(spikes_chan5_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(1,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') ;
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(1,:),50,'unbiased') ;
                        elseif sum(spikes_chan6_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(1,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(1,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))/2);
                        elseif sum(spikes_chan5_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(1,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan6_less90_2s{a,b,c}(d,:))/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(1,:),50,'unbiased');
                        elseif sum(spikes_chan5_less90_2s{a,b,c}(d,:))==0 && sum(spikes_chan6_less90_2s{a,b,c}(d,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased');
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(1,:),50,'unbiased') / (sum(spikes_chan6_less90_2s{a,b,c}(1,:))/2);
                        elseif sum(spikes_chan5_less90_2s{a,b,c}(d,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan6_less90_2s{a,b,c}(d,:))/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(1,:),50,'unbiased') / (sum(spikes_chan6_less90_2s{a,b,c}(1,:))/2);
                        elseif sum(spikes_chan6_less90_2s{a,b,c}(d,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(1,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))*sum(spikes_chan6_less90_2s{a,b,c}(1,:))/4)^(1/2);
                        elseif sum(spikes_chan6_less90_2s{a,b,c}(1,:))==0
                            xcorr_cells_less90_init{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(d,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))*sum(spikes_chan6_less90_2s{a,b,c}(d,:))/4)^(1/2);
                            xcorr_cells_less90_shift{a,b,c}(d,:) = xcorr(spikes_chan5_less90_2s{a,b,c}(d,:),spikes_chan6_less90_2s{a,b,c}(1,:),50,'unbiased') / (sum(spikes_chan5_less90_2s{a,b,c}(d,:))/2);
                        end
                    end
                end
                size_less90_shift=size(xcorr_cells_less90_shift{a,b,c});
                if size_less90_shift(1)>1
                    xcorr_cells_less90_shift_mean{a,b,c}=mean(xcorr_cells_less90_shift{a,b,c});%takes mean xcorr shift within stim conditions
                elseif size_less90_shift(1)==1
                    xcorr_cells_less90_shift_mean{a,b,c}=xcorr_cells_less90_shift{a,b,c};
                end
            end
        end
    end
    for b=1:2%iterate across elevation 90s
        for d = 1 : length(spikes_chan5_90_2s{a,b}(:,1))%iterate across trials
            if d<length(spikes_chan5_90_2s{a,b}(:,1))%calculate xcorr across channel within cells(stimulus conditions) and correct with next trial correlation
                if sum(spikes_chan5_90_2s{a,b}(d,:))>0 && sum(spikes_chan6_90_2s{a,b}(d,:))>0 && sum(spikes_chan6_90_2s{a,b}(d+1,:))>0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))*sum(spikes_chan6_90_2s{a,b}(d,:))/4)^(1/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d+1,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))*sum(spikes_chan6_90_2s{a,b}(d+1,:))/4)^(1/2);
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(d+1,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') ;
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d+1,:),50,'unbiased') ;
                elseif sum(spikes_chan6_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(d+1,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d+1,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))/2);
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(d,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') ;
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d+1,:),50,'unbiased') / (sum(spikes_chan6_90_2s{a,b}(d+1,:))/2);
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(d+1,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan6_90_2s{a,b}(d,:))/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d+1,:),50,'unbiased') ;
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan6_90_2s{a,b}(d,:))/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d+1,:),50,'unbiased') / (sum(spikes_chan6_90_2s{a,b}(d+1,:))/2);
                elseif sum(spikes_chan6_90_2s{a,b}(d,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d+1,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))*sum(spikes_chan6_90_2s{a,b}(d+1,:))/4)^(1/2);
                elseif sum(spikes_chan6_90_2s{a,b}(d+1,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))*sum(spikes_chan6_90_2s{a,b}(d,:))/4)^(1/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d+1,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))/2);
                end
            elseif d==length(spikes_chan5_90_2s{a,b}(:,1))%boundary condition for last trial in cell, correlation correctetd against first trial in opposite channel
                if sum(spikes_chan5_90_2s{a,b}(d,:))>0 && sum(spikes_chan6_90_2s{a,b}(d,:))>0 && sum(spikes_chan6_90_2s{a,b}(1,:))>0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))*sum(spikes_chan6_90_2s{a,b}(d,:))/4)^(1/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))*sum(spikes_chan6_90_2s{a,b}(1,:))/4)^(1/2);
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))>0 && sum(spikes_chan6_90_2s{a,b}(d,:))>0 && sum(spikes_chan6_90_2s{a,b}(1,:))>0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') ;
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') ;
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(1,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') ;
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') ;
                elseif sum(spikes_chan6_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(1,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))/2);
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(d,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') ;
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') / (sum(spikes_chan6_90_2s{a,b}(1,:))/2);
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))==0 && sum(spikes_chan6_90_2s{a,b}(1,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan6_90_2s{a,b}(d,:))/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') ;
                elseif sum(spikes_chan5_90_2s{a,b}(d,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan6_90_2s{a,b}(d,:))/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') / (sum(spikes_chan6_90_2s{a,b}(1,:))/2);
                elseif sum(spikes_chan6_90_2s{a,b}(d,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))*sum(spikes_chan6_90_2s{a,b}(1,:))/4)^(1/2);
                elseif sum(spikes_chan6_90_2s{a,b}(1,:))==0
                    xcorr_cells_90_init{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(d,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))*sum(spikes_chan6_90_2s{a,b}(d,:))/4)^(1/2);
                    xcorr_cells_90_shift{a,b}(d,:) = xcorr(spikes_chan5_90_2s{a,b}(d,:),spikes_chan6_90_2s{a,b}(1,:),50,'unbiased') / (sum(spikes_chan5_90_2s{a,b}(d,:))/2);
                end
            end
        end
        size_90_shift = size(xcorr_cells_90_shift{a,b});
        if size_90_shift>1
            xcorr_cells_90_shift_mean{a,b}=mean(xcorr_cells_90_shift{a,b});%takes mean xcorr shift within stim conditions
        elseif size_90_shift==1
            xcorr_cells_90_shift_mean{a,b}=xcorr_cells_90_shift{a,b}
        end
    end
end

%divide up cells by stim type
%initials
for b=1:length(unique_elevation)-2%iterate across elevations less the 90s
    for c=1:length(unique_azimuth)%iterate across azimuth directions
            xcorr_cells_less90_vestib_init{b,c}=xcorr_cells_less90_init{1,b,c};
            xcorr_cells_less90_visual_init{b,c}=xcorr_cells_less90_init{2,b,c};
            xcorr_cells_less90_combin_init{b,c}=xcorr_cells_less90_init{3,b,c};
    end
end
for b=1:2%iterate across 90s elevation
        xcorr_cells_90_vestib_init{b}=xcorr_cells_90_init{1,b};
        xcorr_cells_90_visual_init{b}=xcorr_cells_90_init{2,b};
        xcorr_cells_90_combin_init{b}=xcorr_cells_90_init{3,b};
end
%remove empty cells
veindex=1;
viindex=1;
coindex=1;
for b=1:length(unique_elevation)-2
    for c=1:length(unique_azimuth)
        if isempty(xcorr_cells_less90_vestib_init{b,c})==0
            size_less90_vestib=size(xcorr_cells_less90_vestib_init{b,c});
            if size_less90_vestib(1)>1
                for z=1:size_less90_vestib(1)
                    xcorr_less90_vestib_init(b,veindex,:)=xcorr_cells_less90_vestib_init{b,c}(z,:);
                    veindex=veindex+1;
                end
            else
                xcorr_less90_vestib_init(b,veindex,:)=xcorr_cells_less90_vestib_init{b,c};
                veindex=veindex+1;
            end
        end
        if isempty(xcorr_cells_less90_visual_init{b,c})==0
            size_less90_visual=size(xcorr_cells_less90_visual_init{b,c});
            if size_less90_visual(1)>1
                for z=1:size_less90_visual(1)
                    xcorr_less90_visual_init(b,veindex,:)=xcorr_cells_less90_visual_init{b,c}(z,:);
                    viindex=viindex+1;
                end
            else
                xcorr_less90_visual_init(b,viindex,:)=xcorr_cells_less90_visual_init{b,c};
                viindex=viindex+1;
            end
        end
        if isempty(xcorr_cells_less90_combin_init{b,c})==0
            size_less90_combin=size(xcorr_cells_less90_combin_init{b,c});
            if size_less90_combin(1)>1
                for z=1:size_less90_combin(1)
                    xcorr_less90_combin_init(b,veindex,:)=xcorr_cells_less90_combin_init{b,c}(z,:);
                    coindex=coindex+1;
                end
            else
                xcorr_less90_combin_init(b,veindex,:)=xcorr_cells_less90_combin_init{b,c};
                coindex=coindex+1;
            end
        end 
    end
end
veindex=1;
viindex=1;
coindex=1;
for b=1:2
    if isempty(xcorr_cells_90_vestib_init{b})==0
        size_90_vestib=size(xcorr_cells_90_vestib_init{b});
        if size_90_vestib(1)>1
            for z=1:size_90_vestib(1)
                xcorr_90_vestib_init(veindex,:)=xcorr_cells_90_vestib_init{b}(z,:);
                veindex=veindex+1;
            end
        else
            xcorr_90_vestib_init(veindex,:)=xcorr_cells_90_vestib_init{b};
            veindex=veindex+1;
        end
    end
    if isempty(xcorr_cells_90_visual_init{b})==0
        size_90_visual=size(xcorr_cells_90_visual_init{b});
        if size_90_visual(1)>1
            for z=1:size_90_visual(1)
                xcorr_90_visual_init(viindex,:)=xcorr_cells_90_visual_init{b}(z,:);
                viindex=viindex+1;
            end
        else
            xcorr_90_visual_init(viindex,:)=xcorr_cells_90_visual_init{b};
            viindex=viindex+1;
        end
    end
    if isempty(xcorr_cells_90_combin_init{b})==0
        size_90_combin=size(xcorr_cells_90_combin_init{b});
        if size_90_combin(1)>1
            for z=1:size_90_combin(1)
                xcorr_90_combin_init(coindex,:)=xcorr_cells_90_combin_init{b}(z,:);
                coindex=coindex+1;
            end
        else
            xcorr_90_combin_init(coindex,:)=xcorr_cells_90_combin_init{b};
            coindex=coindex+1;
        end
    end
end
%shifts
for b=1:length(unique_elevation)-2%iterate across elevations less the 90s
    for c=1:length(unique_azimuth)%iterate across azimuth directions
            xcorr_cells_less90_shift_vestib_mean{b,c}=xcorr_cells_less90_shift_mean{1,b,c};
            xcorr_cells_less90_shift_visual_mean{b,c}=xcorr_cells_less90_shift_mean{2,b,c};
            xcorr_cells_less90_shift_combin_mean{b,c}=xcorr_cells_less90_shift_mean{3,b,c};
    end
end
for b=1:2%iterate across 90s elevation
        xcorr_cells_90_shift_vestib_mean{b}=xcorr_cells_90_shift_mean{1,b};
        xcorr_cells_90_shift_visual_mean{b}=xcorr_cells_90_shift_mean{2,b};
        xcorr_cells_90_shift_combin_mean{b}=xcorr_cells_90_shift_mean{3,b};
end
%remove empty cells
veindex=1;
viindex=1;
coindex=1;
for b=1:length(unique_elevation)-2
    for c=1:length(unique_azimuth)
        if isempty(xcorr_cells_less90_shift_vestib_mean{b,c})==0
            xcorr_less90_vestib_shift(b,veindex,:)=xcorr_cells_less90_shift_vestib_mean{b,c};
            veindex=veindex+1;
        end
        if isempty(xcorr_cells_less90_shift_visual_mean{b,c})==0
            xcorr_less90_visual_shift(b,viindex,:)=xcorr_cells_less90_shift_visual_mean{b,c};
            viindex=viindex+1;
        end
        if isempty(xcorr_cells_less90_shift_combin_mean{b,c})==0
            xcorr_less90_combin_shift(b,coindex,:)=xcorr_cells_less90_shift_combin_mean{b,c};
            coindex=coindex+1;
        end         
    end
end
veindex=1;
viindex=1;
coindex=1;
for b=1:2
    if isempty(xcorr_cells_90_shift_vestib_mean{b})==0
        xcorr_90_vestib_shift(veindex,:)=xcorr_cells_90_shift_vestib_mean{b};
        veindex=veindex+1;
    end
    if isempty(xcorr_cells_90_shift_visual_mean{b})==0
        xcorr_90_visual_shift(viindex,:)=xcorr_cells_90_shift_visual_mean{b};
        viindex=viindex+1;
    end
    if isempty(xcorr_cells_90_shift_combin_mean{b})==0
        xcorr_90_combin_shift(coindex,:)=xcorr_cells_90_shift_combin_mean{b};
        coindex=coindex+1;
    end
end
%--------------------------------------------------------------------------
% make shift predictor based on stim type

for b=1:length(unique_elevation)-2
    xcorr_shift_vestib_mean_temp(b,:)=mean(xcorr_less90_vestib_shift(b,:,:));
    xcorr_shift_visual_mean_temp(b,:)=mean(xcorr_less90_visual_shift(b,:,:));
    xcorr_shift_combin_mean_temp(b,:)=mean(xcorr_less90_combin_shift(b,:,:));
end
size_xcorr_90_ve_shift=size(xcorr_90_vestib_shift);
size_xcorr_90_vi_shift=size(xcorr_90_visual_shift);
size_xcorr_90_co_shift=size(xcorr_90_combin_shift);
if size_xcorr_90_ve_shift(1)>1
    xcorr_shift_vestib_mean_temp(end,:) = mean(xcorr_90_vestib_shift);
else
    xcorr_shift_vestib_mean_temp(end,:) = xcorr_90_vestib_shift;
end
if size_xcorr_90_vi_shift(1)>1
    xcorr_shift_visual_mean_temp(end,:) = mean(xcorr_90_visual_shift);
else
    xcorr_shift_visual_mean_temp(end,:) = xcorr_90_visual_shift;
end
if size_xcorr_90_co_shift(1)>1
    xcorr_shift_combin_mean_temp(end,:) = mean(xcorr_90_combin_shift);
else
    xcorr_shift_combin_mean_temp(end,:) = xcorr_90_combin_shift;
end
xcorr_shift_vestib_mean = mean(xcorr_shift_vestib_mean_temp);
xcorr_shift_visual_mean = mean(xcorr_shift_visual_mean_temp);
xcorr_shift_combin_mean = mean(xcorr_shift_combin_mean_temp);
%--------------------------------------------------------------------------
% compute xcorr true by vestibular/visual/combined conditions
size_init_xcorr_less90_vestib = size(xcorr_less90_vestib_init);
size_init_xcorr_less90_visual = size(xcorr_less90_visual_init);
size_init_xcorr_less90_combin = size(xcorr_less90_combin_init);
size_init_xcorr_90_vestib = size(xcorr_90_vestib_init);
size_init_xcorr_90_visual = size(xcorr_90_visual_init);
size_init_xcorr_90_combin = size(xcorr_90_combin_init);

for b=1:length(unique_elevation)-2
    for z=1:size_init_xcorr_less90_vestib(2)
        xcorr_less90_vestib_corrected(b,z,:)=squeeze(xcorr_less90_vestib_init(b,z,:))-squeeze(xcorr_shift_vestib_mean)';
    end
    for z=1:size_init_xcorr_less90_visual(2)
        xcorr_less90_visual_corrected(b,z,:)=squeeze(xcorr_less90_visual_init(b,z,:))-squeeze(xcorr_shift_visual_mean)';
    end
    for z=1:size_init_xcorr_less90_combin(2)
        xcorr_less90_combin_corrected(b,z,:)=squeeze(xcorr_less90_combin_init(b,z,:))-squeeze(xcorr_shift_combin_mean)';
    end    
end
for z=1:size_init_xcorr_90_vestib(1)
    xcorr_90_vestib_corrected(z,:)=xcorr_90_vestib_init(z,:)-xcorr_shift_vestib_mean;
end
for z=1:size_init_xcorr_90_visual(1)
    xcorr_90_visual_corrected(z,:)=xcorr_90_visual_init(z,:)-xcorr_shift_visual_mean;
end
for z=1:size_init_xcorr_90_combin(1)
    xcorr_90_combin_corrected(z,:)=xcorr_90_combin_init(z,:)-xcorr_shift_combin_mean;
end
%xcorr mean
for b=1:length(unique_elevation)-2
    xcorr_less90_vestib_corrected_mean1(b,:)=mean(xcorr_less90_vestib_corrected(b,:,:));
    xcorr_less90_visual_corrected_mean1(b,:)=mean(xcorr_less90_visual_corrected(b,:,:));
    xcorr_less90_combin_corrected_mean1(b,:)=mean(xcorr_less90_combin_corrected(b,:,:));
end
xcorr_less90_vestib_corrected_mean2=mean(xcorr_less90_vestib_corrected_mean1);
xcorr_less90_visual_corrected_mean2=mean(xcorr_less90_visual_corrected_mean1);
xcorr_less90_combin_corrected_mean2=mean(xcorr_less90_combin_corrected_mean1);
xcorr_90_vestib_corrected_mean=mean(xcorr_90_vestib_corrected);
xcorr_90_visual_corrected_mean=mean(xcorr_90_visual_corrected);
xcorr_90_combin_corrected_mean=mean(xcorr_90_combin_corrected);
xcorr_vestib_corrected_temp(1,:)=xcorr_less90_vestib_corrected_mean2;
xcorr_vestib_corrected_temp(2,:)=xcorr_90_vestib_corrected_mean;
xcorr_visual_corrected_temp(1,:)=xcorr_less90_visual_corrected_mean2;
xcorr_visual_corrected_temp(2,:)=xcorr_90_visual_corrected_mean;
xcorr_combin_corrected_temp(1,:)=xcorr_less90_combin_corrected_mean2;
xcorr_combin_corrected_temp(2,:)=xcorr_90_combin_corrected_mean;
xcorr_vestib_corrected_mean=1000*mean(xcorr_vestib_corrected_temp);
xcorr_visual_corrected_mean=1000*mean(xcorr_visual_corrected_temp);
xcorr_combin_corrected_mean=1000*mean(xcorr_combin_corrected_temp);

figure(3)
clf
hold on
plot(xcorr_vestib_corrected_mean,'k')
plot(xcorr_visual_corrected_mean,'r')
plot(xcorr_combin_corrected_mean,'g')
axis normal
legend('xcorr vestibular', 'xcorr visual', 'xcorr combined')
hold off
% ---------------------------------------------------------------------------------------
% Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s\t'];
for i = 1 : 1000
    sprint_txt = [sprint_txt, ' %1.3f\t'];
end

buff = sprintf(sprint_txt, FILE, xcorr_vestib_corrected_mean, xcorr_visual_corrected_mean, xcorr_combin_corrected_mean);
%buff = sprintf(sprint_txt, FILE, noise_r_stim );

outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\pairwiseunits3DSUSU_timing.dat'];

printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t ');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);

%---------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
return;
% function DirectionTuningPlot_3D_pairwiseunits_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% Path_Defs;
% ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
% 
% %get the column of values for azimuth and elevation and stim_type
% temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
% temp_elevation = data.moog_params(ELEVATION,:,MOOG);
% temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
% 
% %get indices of any NULL conditions (for measuring spontaneous activity
% trials = 1:length(temp_azimuth);
% select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
% null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
% azimuth = temp_azimuth(~null_trials & select_trials);
% elevation = temp_elevation(~null_trials & select_trials);
% stim_type = temp_stim_type(~null_trials & select_trials);
% 
% unique_azimuth = munique(azimuth');
% unique_elevation = munique(elevation');
% unique_stim_type = munique(stim_type');
% 
% % extract channel information
% channelnum_temp = size(data.spike_rates);
% channelnum = channelnum_temp(1,1); % how many channels
% 
% trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-2*(length(unique_azimuth)-1)) * length(unique_stim_type) + 1;
% repetitions = floor( length(find(select_trials==1)) / trials_per_rep );
% %repetitions = 2;
%  mean(data.spike_rates(:, :), 2)
% % for c = 1 : channelnum - 4 %only use from 5th channel
% for c = 1 : 2 %Use only the last 2 channels -- edited AS
%    % temp_spike_rates = data.spike_rates(c+channelnum-2, :); 
%     if c==1
%       temp_spike_rates = data.spike_rates(1,BegTrial : EndTrial);
%    else
%       temp_spike_rates = data.spike_rates(3, BegTrial : EndTrial);
%    end
%     % remove slow drift of responsiveness over time by a high pass filter
%     temp_spike_rates_z = FIR_Filter(temp_spike_rates, 20, 100, 'high', 20, 0);
%     spike_rates_z = temp_spike_rates_z(select_trials);
%     
%     spike_rates = temp_spike_rates(~null_trials & select_trials);
%     spike_rates_channel(c,:) = spike_rates;
%     
%     % create basic matrix represents each response vector
%     resp = [];
%     for k=1:length(unique_stim_type)        
%         for j=1:length(unique_elevation)
%             for i=1:length(unique_azimuth)
%                 select = logical( (azimuth==unique_azimuth(i))  & (elevation==unique_elevation(j))  & (stim_type==unique_stim_type(k)) );
%                 if (sum(select) > 0)  
%                     for t = 1 : repetitions; 
%                         spike_temp = spike_rates(select);   
%                         resp_trial{k}(t,j, i) = spike_temp(t);    
%                         resp_trial_anova1{k}(t,i+(j-1)*8) = resp_trial{k}(t,j,i);    % for later anova use
%                     end
%                     resp(k,j,i) = mean(spike_rates(select));        
%                     resp_std(k,j,i) = std(spike_rates(select));        
% 
%                     % z-score data for spike count correlation analysis
%                     z_dist = spike_rates_z(select);
%                     if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
%                        z_dist = (z_dist - mean(z_dist))/std(z_dist);
%                     else
%                        z_dist = 0;
%                     end
%                     Z_Spikes(select) = z_dist;   
%                 else
%                     for t = 1 : repetitions;                         
%                         resp_trial{k}(t,j, i) = resp_trial{k}(t,j, 1);    
%                         resp_trial_anova1{k}(t,i+(j-1)*8) = resp_trial{k}(t,j,1);    % for later anova use
%                     end
%                     resp(k,j,i) = 0;         % for later vectorsum use
%                     resp_std(k,j,i) = 0;   
%                 end
%             end
%         end  
%         resp_pair{c}(k,1) = resp(k,1,1);
%         resp_pair{c}(k,2:9) = resp(k,2,:);
%         resp_pair{c}(k,10:17) = resp(k,3,:);
%         resp_pair{c}(k,18:25) = resp(k,4,:);
%         resp_pair{c}(k,26) = resp(k,5,1);
%         
%         N=squeeze(resp(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually 
%         [azi(c,k), ele(c,k), amp] = vectorsum(N);
%     end    
%     Z_Spikes_channel(c,:) = Z_Spikes;   
%     
%     % anova1
%     p_1D(c,k) =  anova1( resp_trial_anova1{k}(:,:),'','off' );          
%     % congruency between stim type 
%     if length(unique_stim_type)>=2
%         [rr,pp] = corrcoef(resp_pair{c}(1,:),resp_pair{c}(2,:));
%         corrcoef_r_congruency(c) = rr(1,2);
%         corrcoef_p_congruency(c) = pp(1,2);
%     end
% end
% 
% % now analyze noise correlation between units
% for k=1:length(unique_stim_type) % ananlyze noise correlation in different conditions
%     select_stim = logical( stim_type==unique_stim_type(k)  );
%     select_stim_1D = logical( stim_type==unique_stim_type(k) & elevation==0 );
%     % noise correlation with stim type seperated
%     [rr,pp] = corrcoef(Z_Spikes_channel(1,select_stim),Z_Spikes_channel(2,select_stim));
%     noise_r_stim(k) = rr(1,2);
%     noise_p_stim(k) = pp(1,2);  
%     
%     [rr,pp] = corrcoef(Z_Spikes_channel(1,select_stim_1D),Z_Spikes_channel(2,select_stim_1D));
%     noise_r_stim_1D(k) = rr(1,2);
%     noise_p_stim_1D(k) = pp(1,2);
%     % this is only the regular correlation between two tuning curves
%     [rr,pp] = corrcoef(resp_pair{1}(k,:),resp_pair{2}(k,:));
%     corrcoef_r_unit(k) = rr(1,2);
%     corrcoef_p_unit(k) = pp(1,2);
% end
% 
% figure(2);
% set(2,'Position', [5,15 980,650], 'Name', '1D Direction Tuning');
% orient landscape;
% set(0, 'DefaultAxesXTickMode', 'auto', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'auto');
% % subplot(2,2,1);
% plot(Z_Spikes_channel(1,stim_type==1), Z_Spikes_channel(2,stim_type==1), 'bo');
% title([FILE(1:end-4)  num2str(noise_r_stim(1))]);
% saveas(gcf,  ['C:\work\noise\'  FILE(1:end-4)],  'png');
% close(gcf);
% % subplot(2,2,2);
% % plot(Z_Spikes_channel(1,stim_type==2), Z_Spikes_channel(2,stim_type==2), 'ro');
% % title(num2str(noise_r_stim(2)));
% % subplot(2,2,3);
% % plot(Z_Spikes_channel(1,stim_type==3), Z_Spikes_channel(2,stim_type==3), 'go');
% % title(num2str(noise_r_stim(3)));
% 
% %% ---------------------------------------------------------------------------------------
% % Also, write out some summary data to a cumulative summary file
% % sprint_txt = ['%s\t'];
% % for i = 1 : 1000
% %      sprint_txt = [sprint_txt, ' %1.4f\t\t'];    
% % end
% % 
% % buff = sprintf(sprint_txt, FILE, corrcoef_r_unit, corrcoeff_p_unit, corrcoef_r_congruency, corrcoef_p_congruency, noise_r_stim, noise_p_stim, noise_r_stim_1D, noise_p_stim_1D );
% % %buff = sprintf(sprint_txt, FILE, noise_r_stim );
% % 
% % outfile = ['Z:\Users\Adhira\Spike2_rotation3D_noisecorr.dat'];
% %     
% % printflag = 0;
% % if (exist(outfile, 'file') == 0)    %file does not yet exist
% %     printflag = 1;
% % end
% % fid = fopen(outfile, 'a');
% % if (printflag)
% %     fprintf(fid, 'FILE\t ');
% %     fprintf(fid, '\r\n');
% % end
% % fprintf(fid, '%s', buff);
% % fprintf(fid, '\r\n');
% % fclose(fid);
% 
% % outfile = ['Z:\Users\Adhira\DirectionTuning3D_noisecorr_MU.dat'];
% outfile = ['C:\work\noise\noisecorr_tran1.txt']; %temporarily
% 
% fid = fopen(outfile, 'a');
% % if (printflag)
% %     fprintf(fid, 'FILE\t stim \t signal_corr\t signal_corr_p\t noise_corr\t noise_corr_p');
% %     fprintf(fid, '\r\n');
% % end
% 
% sprint_txt = [];
% for k=1:length(unique_stim_type)
%  %       buff = sprintf(sprint_txt, FILE, unique_stim_type(k), corrcoef_r_unit(k), corrcoef_p_unit(k), corrcoef_r_congruency(1), corrcoef_r_congruency(2), noise_r_stim(k), noise_p_stim(k), noise_r_stim_1D(k), noise_p_stim_1D(k));
%  %       buff = sprintf(sprint_txt, FILE, unique_stim_type(k), corrcoef_r_congruency(:),corrcoef_p_congruency(:) );
% %         buff = sprintf(sprint_txt, FILE, corrcoef_r_unit(k), corrcoef_p_unit(k),noise_r_stim(k), noise_p_stim(k), p_1D(:,k) );
% %         fprintf(fid, buff);
%          fprintf(fid, '\r\n');
%          fprintf(fid,  '%10s\t %10d\t  %5f\t %5f\t %5f\t %5f\t %5f\t %\5f\t', FILE(1:end-4), Protocol, corrcoef_r_unit(k), corrcoef_p_unit(k),noise_r_stim(k), noise_p_stim(k));
% %         fprintf(fid, '\r\n');
% end
% 
% fclose(fid);
% 
% svmat = ['C:\work\noise\' FILE(1:end-4)];
% eval(['save ' svmat ' resp_pair Z_Spikes_channel  select_stim  select_stim_1D noise_r_stim noise_p_stim noise_r_stim_1D noise_p_stim_1D corrcoef_r_unit corrcoef_p_unit p_1D']); 
% 
% svmat = ['C:\work\noise\1\' FILE(1:7) '_tran'];
% eval(['save ' svmat ' resp_pair Z_Spikes_channel  select_stim  select_stim_1D noise_r_stim noise_p_stim noise_r_stim_1D noise_p_stim_1D corrcoef_r_unit corrcoef_p_unit p_1D']); 
% 
% %---------------------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% return;