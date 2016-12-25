% Trapezoidal velocity profile -GY 07/30/2010
% %-----------------------------------------------------------------------------------------------------------------------
function AZIMUTH_TUNING_1D_TRAP_Tuning_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);

%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_azimuth);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';

spon_found = find(null_trials==1); 
%spon_resp = mean(temp_spike_rates(spon_found))
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
temp_spike_data = data.spike_data(:, :);
for i = 1 : length(Discard_trials)
    temp_spike_data( :, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 99;
end

SpikeChan(1) = 1;
%SpikeChan(2) = 3;
resp = [];
for c=1:length(SpikeChan)
    spike_data(1,:) = temp_spike_data( SpikeChan(c), find(temp_spike_data(1,:)~=99) );
    spike_data(1, find(spike_data>1) ) = 1; % something is absolutely wrong  
    for ss =  1 : length(azimuth) % ss marks the index of trial
       spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+100+115+5000*(ss-1) : StartEventBin(1)+200+115+800+5000*(ss-1)) ) ; 
    end
    spike_rates = spike_rates*2; % mean firing rate
    % creat basic matrix represents each response vector
    
    for k=1:length(unique_stim_type)
        for i=1:length(unique_azimuth)
            select = logical( azimuth==unique_azimuth(i) & stim_type==unique_stim_type(k) );
            repetition(i,k) = length(find(select==1));
            resp{c}(i, k) = mean(spike_rates(select));        
            resp_std(i,k) = std(spike_rates(select));        
            resp_err{c}(i,k) = std(spike_rates(select)) / sqrt(repetition(i,k));
        end
    end    
    repetition_min = min(min(repetition));
    % compute anova1 p value
    for k=1:length(unique_stim_type)
        for i=1:length(unique_azimuth)
            select = logical( azimuth==unique_azimuth(i) & stim_type==unique_stim_type(k) );
            spike_temp = spike_rates(select);
            resp_trial{k}(1:repetition_min,i) = spike_temp(1:repetition_min);
        end
        p_1D(c,k) = anova1(resp_trial{k},'','off');
    end
end
repetition_min
% Define figure
xoffset=0;
yoffset=0;
figure(2);
set(2,'Position', [5,15 980,650], 'Name', '1D Direction Tuning trapezoid');
orient landscape;
set(0, 'DefaultAxesXTickMode', 'manual', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'manual');

% temporarily hard coded, will be probematic if there are more than 3*3 conditions
% repetitions -GY
f{1,1}='bo-'; f{1,2}='bo-'; f{1,3}='bo-'; 
f{2,1}='ro-'; f{2,2}='ro-'; f{2,3}='ro-'; 
f{3,1}='go-'; f{3,2}='go-'; f{3,3}='go-'; 
for c=1:length(SpikeChan)
    for k=1: length(unique_stim_type) 
        axes('position',[0.1+0.3*(k-1) 0.6-0.42*(c-1) 0.25 0.3]);
        errorbar(unique_azimuth, resp{c}(:,k), resp_err{c}(:,k), 'o-' );

        ylabel('spikes/s');
        xlabel('azimuth');
        xlim( [0,360] );
        title( ['chanel' num2str(SpikeChan(c)) ' ' h_title{unique_stim_type(k)} ' ' 'p' num2str(p_1D(c,k))] );
        set(gca, 'xtick', unique_azimuth);  
    end
end
%show file name and some values in text
axes('position',[0.1,0.85, 0.1,0.1] );
xlim( [0,100] );
ylim( [0,1] );
text(0, 1, FILE);

axis off;

% %% ---------------------------------------------------------------------------------------
% % Also, write out some summary data to a cumulative summary file
% sprint_txt = ['%s\t'];
% for i = 1 : 200
%      sprint_txt = [sprint_txt, ' %1.3f\t'];    
% end
% % if length(unique_stim_type)~=1
%     buff = sprintf(sprint_txt, FILE, spon_resp, az, p_1D, Dprime, DDI );
% %    buff = sprintf(sprint_txt, FILE, unique_stim_type, unique_motion_coherence, DDI );
%     outfile = [BASE_PATH 'ProtocolSpecific\MOOG\AzimuthTuning1D\DirectionTuning1D.dat'];
% % else    
% %     buff = sprintf(sprint_txt, FILE, spon_resp, az(:), amp(:), p_1D{:}, DI);
% %     outfile = [BASE_PATH 'ProtocolSpecific\MOOG\AzimuthTuning1D\DirectionTuning1D_Hui.dat'];
% % end
% 
% % buff = sprintf(sprint_txt, FILE, az(:), p_1D{:},congruency, Z_Spikes );
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\AzimuthTuning1D\DirectionTuning1D_Zscore.dat'];
%     
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
% 	if length(unique_stim_type)~=1
%         %     fprintf(fid, 'FILE\t Spon\t HTIve\t HTIvi\t VesP\t VisP\t VesPre\t VisPre\t VesSlo\t Vis\Slo\t VesMin\t VisMin\t VesMax\t VisMax\t');
%         fprintf(fid, 'FILE\t Spon\t VesPref\t VisPref\t VesP\t VisP\t VesDDI\t VisDDI\t Congruency');
% 	else
%         fprintf(fid, 'FILE\t Spon\t Az_6\t Az_7_5\t Az_9\t Amp_6\t Amp_7_5\t Amp_9\t P_6\t P_7.5\t P_9\t DI_6\t DI_7.5\t DI_9\t');
% 	end
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% %---------------------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% % SaveTrials(FILE,BegTrial,EndTrial,p_1D)
return;

