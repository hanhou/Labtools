function [CorrSepOri CorrSepRemoved]=SpikeCorrSep(select_data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,FigureIndex);
%SpikeCorrSep.m - this function plot the cross correlation between the SU and MU (seperate for different stimulus condition: vestibular, visual or combined)

%*************************************************************************%
%seperate the trials for different stimulus: visual,vestibular and combined
Path_Defs;
ProtocolDefs;
SpikeChan=1
temp_stim_type = select_data.moog_params(STIM_TYPE,:,MOOG);
switch (Protocol)
    case(DIRECTION_TUNING_3D)
        temp_azimuth = select_data.moog_params(AZIMUTH,:,MOOG);  
    case(ROTATION_TUNING_3D)
        temp_azimuth = select_data.moog_params(ROT_AZIMUTH,:,MOOG);
end

null_trials = logical( (temp_azimuth == select_data.one_time_params(NULL_VALUE)) );
temp_spike_rates = select_data.spike_rates(SpikeChan, :);  
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_tri = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_tri ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_tri) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end
stim_type = temp_stim_type(~null_trials & select_trials);
unique_stim_type = munique(stim_type');
condition_num = stim_type;
unique_condition_num = munique(condition_num');%munique(temp_stim_type')

h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
%*************************************************************************%
%SpikeChan=1;
SpikeChan2=4;
correl_total = 0;
for i=1:length(unique_condition_num)
    select_rep=find(condition_num==unique_condition_num(i));
    for j=1:length(select_rep)
        spike_data1=select_data.spike_data(SpikeChan,:,select_rep(j));
        spike_data2=select_data.spike_data(SpikeChan2,:,select_rep(j));
        [correl,lags]=xcorr(spike_data1,spike_data2,50);
        correl_total=correl+correl_total;
    end
    [val, ind] = max(correl_total);
    %choptime = lags(ind);
    figure(FigureIndex);subplot(2,2,i);
    plot(lags, correl_total,'r');title(h_title{i});    
    CorrSepRemoved(i,:)=correl_total;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot the cross-correlation between SU and MU (after the SU was removed)
% %SpikeChan=1;
% SpikeChan2=4;
% num_trials = size(select_data.spike_data);
% num_trials = num_trials(3);
% correl_total = 0;
% for i=2:num_trials
%    spike_data1 = select_data.spike_data(SpikeChan, :, i);
%    spike_data2 = select_data.spike_data(SpikeChan2, :, i);
%    [correl, lags] = xcorr(spike_data1, spike_data2, 50);%    [correl, lags] = xcorr(spike_data1, spike_data2, 50,'coeff');   
%    correl_total = correl+correl_total;
% end
% % correl_total=(correl_total)/num_trials;%Add by AHC 04-12-06
% [val, ind] = max(correl_total);
% %choptime = lags(ind);
% %hold on; 
% figure(4);clf;
% plot(lags, correl_total,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx=select_data.spike_data(1:2,:,:);
select_data.spike_data=xx;clear xx;
[select_data winDiscrim] = Packaroni(select_data, select_data.fName, 1, 0);  %get window discriminator data
[select_data shiftValue]= SpikeCorrC(select_data, 1, 1, winDiscrim);
if  ( (shiftValue == -50) & (sum( sum( select_data.spike_data(winDiscrim,:,:) ) == 0) ) )
    % remove window discrimination data if not being used
    select_data.spike_data = selct_data.spike_data(1:size(select_data.spike_data,1) - 1,:,:);
    shiftValue = 0;
end
[select_data newChans] = Packaroni(select_data, select_data.fName, 0, shiftValue);  %get rest of data

%*************************************************************************%
%plot the cross-correlation between SU and MU (before the SU was removed)
% SpikeChan=1;SpikeChan2=4;
% num_trials = size(select_data.spike_data);
% num_trials = num_trials(3);
% correl_total = 0;
% for i=2:num_trials
%    spike_data1 = select_data.spike_data(SpikeChan, :, i);
%    spike_data2 = select_data.spike_data(SpikeChan2, :, i);
%    [correl, lags] = xcorr(spike_data1, spike_data2, 50);%    [correl, lags] = xcorr(spike_data1, spike_data2, 50,'coeff');   
%    correl_total = correl+correl_total;
% end
% [val, ind] = max(correl_total);
% choptime = lags(ind);
% figure(4);hold on;
% plot(lags, correl_total);

for i=1:length(unique_condition_num)
    select_rep=find(condition_num==unique_condition_num(i));
    for j=1:length(select_rep)
        spike_data1=select_data.spike_data(SpikeChan,:,select_rep(j));
        spike_data2=select_data.spike_data(SpikeChan2,:,select_rep(j));
        [correl,lags]=xcorr(spike_data1,spike_data2,50);
        correl_total=correl+correl_total;
    end
    [val, ind] = max(correl_total);
    %choptime = lags(ind);
    figure(FigureIndex);subplot(2,2,i);hold on;
    plot(lags, correl_total);title(h_title{i});    
    CorrSepOri(i,:)=correl_total;
    legend('Before','SU Removed');
end
