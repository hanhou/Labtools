%-----------------------------------------------------------------------------------------------------------------------
%-- CueDirecEffects.m -- RMS LFP activity and spikes during delay period and during stimulus for various cue directions.
%--	VR, 9/21/05 
%-----------------------------------------------------------------------------------------------------------------------

function CueDirecEffects(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

global cum_delay_LFP cum_delay_spikes cum_delay_LFP_bp; %to allow activity to be cumulated across trials
SAVE_GLOBAL_DATA = 1;

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,BegTrial:EndTrial,PATCH1);
unique_direction = munique(direction');
Pref_direction = data.one_time_params(PREFERRED_DIRECTION);

%get the motion coherences
coherence = data.dots_params(DOTS_COHER, BegTrial:EndTrial, PATCH1);
unique_coherence = munique(coherence');

%get the cue validity: -1=Invalid; 0=Neutral; 1=Valid; 2=CueOnly
cue_val = data.cue_params(CUE_VALIDITY,BegTrial:EndTrial,PATCH2);
unique_cue_val = munique(cue_val');
cue_val_names = {'NoCue','Invalid','Neutral','Valid','CueOnly'};

%get the cue directions
cue_direc = data.cue_params(CUE_DIREC, BegTrial:EndTrial, PATCH1);
unique_cue_direc = munique(cue_direc');

%compute cue types - 0=neutral, 1=directional, 2=cue_only
cue_type = abs(cue_val); %note that invalid(-1) and valid(+1) are directional
unique_cue_type = munique(cue_type');

%now, select trials that fall between BegTrial and EndTrial
trials = 1:size(data.dots_params,2);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

nreps = floor(length(cue_type)/60); %assumes the total number of trials is a multiple of 60. 

%get the firing rates and lfp during delay and stim for all the trials
stim_spikes = data.spike_rates(SpikeChan, BegTrial:EndTrial);
if (isempty(data.lfp_data)) %in case the lfp data wasn't saved, fill a matrix with zeros so that the other analyses can occur
    data.lfp_data = zeros(size(data.spike_data(1,:,BegTrial:EndTrial)));
    SAVE_GLOBAL_DATA = 0;
end
for i = 1:sum(select_trials)
    start_delay(i) = find(data.event_data(1,:,i+BegTrial-1) == CUE_ON_CD);
    end_delay(i) = find(data.event_data(1,:,i+BegTrial-1) == VSTIM_ON_CD);
    delay_spikes(i) = sum(data.spike_data(SpikeChan,start_delay(i):end_delay(i),i+BegTrial-1)) / length(start_delay(i):end_delay(i)) * 1000;
    %note that lfp is sampled at half the frequency as spikes, so divide bins by 2
    delay_lfp(i) = sqrt(mean( data.lfp_data(1,ceil(start_delay(i)/2):floor(end_delay(i)/2),i+BegTrial-1).^2 )); 
    start_stim(i) = ceil(end_delay(i)/2);
    end_stim(i) = floor(find(data.event_data(1,:,i+BegTrial-1) == VSTIM_OFF_CD)/2);
    stim_lfp(i) = sqrt(mean( data.lfp_data(1,start_stim(i):end_stim(i),i+BegTrial-1) .^2 ));
    
    %do the following to get the power of lfp between 50 and 150Hz 
    %(remove 120 Hz contribution as noise), 400 samples sampled at 500Hz
    band = find( (500*(0:200)./400 >= 50) & (500*(0:200)./400 <= 150) & (500*(0:200)./400 ~= 120) ); 
    lfp_stim_powerspect{i} = abs(fft(data.lfp_data(1,start_stim(i):end_stim(i),i+BegTrial-1),400)).^2 ./ 400;
    stim_lfp_bp(i) = sum(lfp_stim_powerspect{i}(band));
    lfp_delay_powerspect{i} = abs(fft(data.lfp_data(1,start_delay(i):end_delay(i),i+BegTrial-1),400)).^2 ./ 400;
    delay_lfp_bp(i) = sum(lfp_delay_powerspect{i}(band));
end

%get outcome for each trial: 0=incorrect, 1=correct
trials_outcomes = logical (data.misc_params(OUTCOME,BegTrial:EndTrial) == CORRECT);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (coherence == data.one_time_params(NULL_VALUE)) );


%keyboard

hlist = []; %list of figure handles for saving out graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first perform ANOVA of LFP during delay period using cue direction 

%get trial #s when cue is in preferred and null directions and neutral cue trials
PD_cue_trials = find( (squeeze_angle(cue_direc) == squeeze_angle(Pref_direction)) & (cue_val ~= 0));
ND_cue_trials = find( (squeeze_angle(cue_direc) ~= squeeze_angle(Pref_direction)) & (cue_val ~= 0));
neutral_cue_trials = find(cue_val == 0);

%list of all selected trials marked with -1=NullDirCues, 1=PrefDirCues, %0=remainder(neutral cues)
signed_cue_dirs = zeros(1,sum(select_trials));
signed_cue_dirs(PD_cue_trials') = 1;
signed_cue_dirs(ND_cue_trials') = -1;

cum_delay_LFP = {[] [] []}; cum_delay_spikes = {[] [] []}; cum_delay_LFP_bp = {[] [] []};
if (SAVE_GLOBAL_DATA)
    cum_delay_LFP{1} = [cum_delay_LFP{1} (delay_lfp(PD_cue_trials) - mean(delay_lfp))./std(delay_lfp)];
    cum_delay_LFP_bp{1} = [cum_delay_LFP_bp{1} (delay_lfp_bp(PD_cue_trials) - mean(delay_lfp_bp))./std(delay_lfp_bp)];
    cum_delay_spikes{1} = [cum_delay_spikes{1} (delay_spikes(PD_cue_trials) - mean(delay_spikes))./std(delay_spikes)];
    cum_delay_LFP{2} = [cum_delay_LFP{2} (delay_lfp(neutral_cue_trials) - mean(delay_lfp))./std(delay_lfp)];
    cum_delay_LFP_bp{2} = [cum_delay_LFP_bp{2} (delay_lfp_bp(neutral_cue_trials) - mean(delay_lfp_bp))./std(delay_lfp_bp)];
    cum_delay_spikes{2} = [cum_delay_spikes{2} (delay_spikes(neutral_cue_trials) - mean(delay_spikes))./std(delay_spikes)];
    cum_delay_LFP{3} = [cum_delay_LFP{3} (delay_lfp(ND_cue_trials) - mean(delay_lfp))./std(delay_lfp)];
    cum_delay_LFP_bp{3} = [cum_delay_LFP_bp{3} (delay_lfp_bp(ND_cue_trials) - mean(delay_lfp_bp))./std(delay_lfp_bp)];
    cum_delay_spikes{3} = [cum_delay_spikes{3} (delay_spikes(ND_cue_trials) - mean(delay_spikes))./std(delay_spikes)];
end

%1-way ttest to test whether PDcue > NDcue
[h_delayspikes,p_delayspikes]=ttest2(delay_spikes(PD_cue_trials),delay_spikes(ND_cue_trials),.05,'right');
[h_delaylfp,p_delaylfp]=ttest2(delay_lfp(PD_cue_trials),delay_lfp(ND_cue_trials),.05,'right');
[h_delaylfp_bp,p_delaylfp_bp]=ttest2(delay_lfp_bp(PD_cue_trials),delay_lfp_bp(ND_cue_trials),.05,'right');
[h_stimspikes,p_stimspikes]=ttest2(stim_spikes(PD_cue_trials(find(cue_val(PD_cue_trials)~=CUEONLY))),...
    stim_spikes(ND_cue_trials(find(cue_val(ND_cue_trials)~=CUEONLY))),.05,'right');  %exclude CueOnly trials from stim
[h_stimlfp,p_stimlfp]=ttest2(stim_lfp(PD_cue_trials(find(cue_val(PD_cue_trials)~=CUEONLY))),...
    stim_lfp(ND_cue_trials(find(cue_val(ND_cue_trials)~=CUEONLY))),.05,'right'); 
[h_stimlfp_bp,p_stimlfp_bp]=ttest2(stim_lfp_bp(PD_cue_trials(find(cue_val(PD_cue_trials)~=CUEONLY))),...
    stim_lfp_bp(ND_cue_trials(find(cue_val(ND_cue_trials)~=CUEONLY))),.05,'right'); 

%1-way anova to test whether there is a significant difference among PDcue-neutral-NDcue
p_1way_anova_delaylfp = anovan(delay_lfp,{signed_cue_dirs},'varnames','CueDir','display','off');
p_1way_anova_delaylfp_bp = anovan(delay_lfp_bp,{signed_cue_dirs},'varnames','CueDir','display','off');
p_1way_anova_delayspikes = anovan(delay_spikes,{signed_cue_dirs},'varnames','CueDir','display','off');
p_1way_anova_stimlfp = anovan(stim_lfp(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY))},'varnames','CueDir','display','off');
p_1way_anova_stimlfp_bp = anovan(stim_lfp_bp(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY))},'varnames','CueDir','display','off');
p_1way_anova_stimspikes = anovan(stim_spikes(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY))},'varnames','CueDir','display','off');

%2-way anova with CueDir and Signed Coherence as factors (ignore CueOnly trials)
%if there's an interaction, that would argue that the cue direction modifies the response to the visual stimulus 
%in a way not accounted for by the fact that the stimulus itself is changing. 
%(what about something like a sequential f-test to see whether knowledge of the cue_dir improves fit?)
signed_coherence = coherence.*(-1+2.*(squeeze_angle(direction)==squeeze_angle(Pref_direction)));
p_2way_anova_stimlfp = anovan(stim_lfp(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY)) signed_coherence(find(cue_val~=CUEONLY))},'model',[1 2 3],'varnames',{'CueDir';'SignedCoherence'},'display','off')
p_2way_anova_stimlfp_bp = anovan(stim_lfp_bp(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY)) signed_coherence(find(cue_val~=CUEONLY))},'model',[1 2 3],'varnames',{'CueDir';'SignedCoherence'},'display','off')
p_2way_anova_stimspikes = anovan(stim_spikes(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY)) signed_coherence(find(cue_val~=CUEONLY))},'model',[1 2 3],'varnames',{'CueDir';'SignedCoherence'},'display','off')


%% use zscores to combine across signed_coherences (do this only for stim period, of course)
stimlfp_zscores = stim_lfp; stimspikes_zscores = stim_spikes; stimlfp_bp_zscores = stim_lfp_bp;
%include cueonly trials initially to allow proper dir/coher indexing
for i = 1:length(unique_direction)
    for j = 1:length(unique_coherence)
        indices = find( (direction == unique_direction(i)) & (coherence == unique_coherence(j)) & (cue_val ~= 2) );
        stimlfp_zscores(indices) = zscore(stim_lfp(indices));
        stimlfp_bp_zscores(indices) = zscore(stim_lfp_bp(indices));
        stimspikes_zscores(indices) = zscore(stim_spikes(indices));
    end
end
p_zscore_anova_stimlfp = anovan(stimlfp_zscores(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY))},'varnames','CueDir','display','off');
p_zscore_anova_stimlfp_bp = anovan(stimlfp_bp_zscores(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY))},'varnames','CueDir','display','off');
p_zscore_anova_stimspikes = anovan(stimspikes_zscores(find(cue_val~=CUEONLY)),{signed_cue_dirs(find(cue_val~=CUEONLY))},'varnames','CueDir','display','off');
%now remove cueonly trials
% stimlfp_zscores = stimlfp_zscores(find(cue_val~=CUEONLY)); 
% stimlfp_bp_zscores = stimlfp_bp_zscores(find(cue_val~=CUEONLY)); 
% stimspikes_zscores = stimspikes_zscores(find(cue_val~=CUEONLY));
PD_cue_noCO_trials = PD_cue_trials(cue_val(PD_cue_trials)~=CUEONLY); %remove cueonly trials from list of PD/ND trials
ND_cue_noCO_trials = ND_cue_trials(cue_val(ND_cue_trials)~=CUEONLY);

%% compute mean and std dev values
mean_delay_lfp = [mean(delay_lfp(PD_cue_trials)); mean(delay_lfp(neutral_cue_trials)); mean(delay_lfp(ND_cue_trials))];
std_delay_lfp  = [std(delay_lfp(PD_cue_trials)); std(delay_lfp(neutral_cue_trials)); std(delay_lfp(ND_cue_trials))];
mean_delay_lfp_bp = [mean(delay_lfp_bp(PD_cue_trials)); mean(delay_lfp_bp(neutral_cue_trials)); mean(delay_lfp_bp(ND_cue_trials))];
std_delay_lfp_bp  = [std(delay_lfp_bp(PD_cue_trials)); std(delay_lfp_bp(neutral_cue_trials)); std(delay_lfp_bp(ND_cue_trials))];
mean_delay_spikes = [mean(delay_spikes(PD_cue_trials)); mean(delay_spikes(neutral_cue_trials)); mean(delay_spikes(ND_cue_trials))];
std_delay_spikes =  [std(delay_spikes(PD_cue_trials)); std(delay_spikes(neutral_cue_trials)); std(delay_spikes(ND_cue_trials))];
mean_stim_zlfp = [mean(stimlfp_zscores(PD_cue_noCO_trials)); mean(stimlfp_zscores(neutral_cue_trials)); mean(stimlfp_zscores(ND_cue_noCO_trials))];
std_stim_zlfp =  [std(stimlfp_zscores(PD_cue_noCO_trials)); std(stimlfp_zscores(neutral_cue_trials)); std(stimlfp_zscores(ND_cue_noCO_trials))];
mean_stim_zlfp_bp = [mean(stimlfp_bp_zscores(PD_cue_noCO_trials)); mean(stimlfp_bp_zscores(neutral_cue_trials)); mean(stimlfp_bp_zscores(ND_cue_noCO_trials))];
std_stim_zlfp_bp =  [std(stimlfp_bp_zscores(PD_cue_noCO_trials)); std(stimlfp_bp_zscores(neutral_cue_trials)); std(stimlfp_bp_zscores(ND_cue_noCO_trials))];
mean_stim_zspikes = [mean(stimspikes_zscores(PD_cue_noCO_trials)); mean(stimspikes_zscores(neutral_cue_trials)); mean(stimspikes_zscores(ND_cue_noCO_trials))];
std_stim_zspikes =  [std(stimspikes_zscores(PD_cue_noCO_trials)); std(stimspikes_zscores(neutral_cue_trials)); std(stimspikes_zscores(ND_cue_noCO_trials))];

%% make subplots showing mean values +/- std devs
hlist(1+length(hlist))=figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Cue Direction Effects');
subplot(2,6,7);
errorbar([1:3], mean_delay_lfp, std_delay_lfp./sqrt([length(PD_cue_trials) length(neutral_cue_trials) length(ND_cue_trials)]'));
xlabel('Pr Neu Nu'); set(gca,'XTickLabel',[]); title('Delay: LFP');
subplot(2,6,8);
errorbar([1:3], mean_delay_lfp_bp, std_delay_lfp_bp./sqrt([length(PD_cue_trials) length(neutral_cue_trials) length(ND_cue_trials)]'));
xlabel('Pr Neu Nu'); set(gca,'XTickLabel',[]); title('Delay: BP');
subplot(2,6,9);
errorbar([1:3], mean_delay_spikes, std_delay_spikes./sqrt([length(PD_cue_trials) length(neutral_cue_trials) length(ND_cue_trials)]'));
xlabel('Pr Neu Nu'); set(gca,'XTickLabel',[]); title('Delay: Spikes');
subplot(2,6,10);
errorbar([1:3], mean_stim_zlfp, std_stim_zlfp./sqrt([length(PD_cue_trials) length(neutral_cue_trials) length(ND_cue_trials)]'));
xlabel('Pr Neu Nu'); set(gca,'XTickLabel',[]); title('Stim: LFP');
subplot(2,6,11);
errorbar([1:3], mean_stim_zlfp_bp, std_stim_zlfp_bp./sqrt([length(PD_cue_trials) length(neutral_cue_trials) length(ND_cue_trials)]'));
xlabel('Pr Neu Nu'); set(gca,'XTickLabel',[]); title('Stim: BP');
subplot(2,6,12);
errorbar([1:3], mean_stim_zspikes, std_stim_zspikes./sqrt([length(PD_cue_trials) length(neutral_cue_trials) length(ND_cue_trials)]'));
xlabel('Pr Neu Nu'); set(gca,'XTickLabel',[]); title('Stim: Spikes');

%% ********************** PRINT INFO *****************************
%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%now, print out some specific useful info.
xpos = -10; ypos = 25;
font_size = 8;
bump_size = 5;
line = sprintf('1-tailed T-tests:');
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('  p(delay/lfp)=%5.3f, p(delay/spikes)=%5.3f, p(delay/bp)=%5.3f',p_delaylfp, p_delayspikes, p_delaylfp_bp);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('  p(stim/lfp)=%5.3f, p(stim/spikes)=%5.3f, p(stim/bp)=%5.3f', p_stimlfp, p_stimspikes, p_stimlfp_bp);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Delay period 1-way ANOVAs:');
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('  p(delay/lfp)=%5.3f, p(delay/spikes)=%5.3f, p(delay/bp)=%5.3f', p_1way_anova_delaylfp, p_1way_anova_delayspikes, p_1way_anova_delaylfp_bp);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('VStim period 2-way ANOVAs:');
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('  p(stim/LFP:me)=%5.3f, p(stim/LFP:intrxn)=%5.3f, p(stim/spike:me)=%5.3f, p(stim/spikes:intrxn)=%5.3f',...
    p_2way_anova_stimlfp(1), p_2way_anova_stimlfp(3), p_2way_anova_stimspikes(1), p_2way_anova_stimspikes(3));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('  p(stim/BP:me)=%5.3f, p(stim/BP:intrxn)=%5.3f', p_2way_anova_stimlfp_bp(1), p_2way_anova_stimlfp_bp(3));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('VStim period Z-scored 1-way ANOVAs:');
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('  p(stim/lfp)=%5.3f, p(stim/spikes)=%5.3f, p(stim/bp)=%5.3f', p_zscore_anova_stimlfp, p_zscore_anova_stimspikes, p_zscore_anova_stimlfp_bp);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;


% line = sprintf('Directions tested: %6.3f, %6.3f deg', unique_direction(1), unique_direction(2) );
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;


% %now try to plot timecourses

% PD_lfpdata = squeeze(data.lfp_data(:,:,PD_cue_trials));
% ND_lfpdata = squeeze(data.lfp_data(:,:,ND_cue_trials));
% neutral_lfpdata = squeeze(data.lfp_data(:,:,neutral_cue_trials));
% 
% PD_spikedata = squeeze(data.spike_data(SpikeChan,:,PD_cue_trials));
% ND_spikedata = squeeze(data.spike_data(SpikeChan,:,ND_cue_trials));
% neutral_spikedata = squeeze(data.spike_data(SpikeChan,:,neutral_cue_trials));
% 
% lfp_delay_t = [-400:2:0]; spike_delay_t = [-400:1:0];
% figure;
% subplot(2,2,1); hold on;
% plot(lfp_delay_t,getval(conv( sqrt(mean(PD_lfpdata(ceil(end_delay/2)-200-2:floor(end_delay/2)+1,:).^2,2)), 0.25*ones(1,4)),4:204),'b');
% plot(lfp_delay_t,getval(conv( sqrt(mean(ND_lfpdata(ceil(end_delay/2)-200-2:floor(end_delay/2)+1,:).^2,2)), 0.25*ones(1,4)),4:204),'r');
% plot(lfp_delay_t,getval(conv( sqrt(mean(neutral_lfpdata(ceil(end_delay/2)-200-2:floor(end_delay/2)+1,:).^2,2)), 0.25*ones(1,4)),4:204),'g');
% %xlabel('Time before stimulus (ms)'); 
% ylabel('rms LFP'); 
% %legend('PrefDirec Cue','NullDirec Cue','Neutral Cue');
% 
% subplot(2,2,3); hold on;
% plot(spike_delay_t,getval(conv( sum(PD_spikedata(end_delay-400-2:end_delay+1,:),2), 0.25*ones(1,4)),4:404),'b')
% plot(spike_delay_t,getval(conv( sum(ND_spikedata(end_delay-400-2:end_delay+1,:),2), 0.25*ones(1,4)),4:404),'r')
% plot(spike_delay_t,getval(conv( sum(neutral_spikedata(end_delay-400-2:end_delay+1,:),2), 0.25*ones(1,4)),4:404),'g')
% xlabel('Time before stimulus (ms)'); ylabel('firing rate (imp/s)');
% %legend('PrefDirec Cue','NullDirec Cue','Neutral Cue');
% 
% lfp_stim_t = [0:2:1000]; spike_delay_t = [0:1:1000];
% subplot(2,2,2); hold on;
% plot(lfp_stim_t,getval(conv( sqrt(mean(PD_lfpdata(start_stim-2:start_stim+500+1,:).^2,2)), 0.25*ones(1,4)),4:504),'b');
% plot(lfp_stim_t,getval(conv( sqrt(mean(ND_lfpdata(start_stim-2:start_stim+500+1,:).^2,2)), 0.25*ones(1,4)),4:504),'r');
% plot(lfp_stim_t,getval(conv( sqrt(mean(neutral_lfpdata(start_stim-2:start_stim+500+1,:).^2,2)), 0.25*ones(1,4)),4:504),'g');
% %xlabel('Time before stimulus (ms)'); ylabel('rms LFP'); 
% legend('PrefDirec Cue','NullDirec Cue','Neutral Cue'); legend(gca,'boxoff');
% 
% subplot(2,2,4); hold on;
% plot(spike_delay_t,getval(conv( sum(PD_spikedata(end_delay-6:end_delay+1000+3,:),2), 0.125*ones(1,8)),6:1006),'b')
% plot(spike_delay_t,getval(conv( sum(ND_spikedata(end_delay-6:end_delay+1000+3,:),2), 0.125*ones(1,8)),6:1006),'r')
% plot(spike_delay_t,getval(conv( sum(neutral_spikedata(end_delay-6:end_delay+1000+3,:),2), 0.125*ones(1,8)),6:1006),'g')
% xlabel('Time after stimulus (ms)'); %ylabel('firing rate (imp/s)');
% %legend('PrefDirec Cue','NullDirec Cue','Neutral Cue');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




PD_lfpdata = squeeze(data.lfp_data(:,:,PD_cue_trials));
ND_lfpdata = squeeze(data.lfp_data(:,:,ND_cue_trials));
neutral_lfpdata = squeeze(data.lfp_data(:,:,neutral_cue_trials));

PD_spikedata = squeeze(data.spike_data(SpikeChan,:,PD_cue_trials));
ND_spikedata = squeeze(data.spike_data(SpikeChan,:,ND_cue_trials));
neutral_spikedata = squeeze(data.spike_data(SpikeChan,:,neutral_cue_trials));

lfp_delay_t = [-400:2:0]; spike_delay_t = [-400:1:0];
hlist(1+length(hlist))=figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Cue Direction Effects: TimeCourses');
subplot(4,1,1); hold on;
plot(lfp_delay_t,getval(conv( sqrt(mean(PD_lfpdata(ceil(end_delay/2)-200-2:floor(end_delay/2)+1,:).^2,2)), 0.25*ones(1,4)),4:204),'b');
plot(lfp_delay_t,getval(conv( sqrt(mean(ND_lfpdata(ceil(end_delay/2)-200-2:floor(end_delay/2)+1,:).^2,2)), 0.25*ones(1,4)),4:204),'r');
plot(lfp_delay_t,getval(conv( sqrt(mean(neutral_lfpdata(ceil(end_delay/2)-200-2:floor(end_delay/2)+1,:).^2,2)), 0.25*ones(1,4)),4:204),'g');
title(FILE);
%xlabel('Time before stimulus (ms)'); 
ylabel('rms LFP'); 
legend('PrefDirec Cue','NullDirec Cue','Neutral Cue','Location','NorthOutside','Orientation','Horizontal');legend(gca,'BoxOff')

subplot(4,1,2); hold on;
plot(spike_delay_t,getval(conv( sum(PD_spikedata(end_delay-400-6:end_delay+3,:),2), 0.125*ones(1,8)),6:406),'b')
plot(spike_delay_t,getval(conv( sum(ND_spikedata(end_delay-400-6:end_delay+3,:),2), 0.125*ones(1,8)),6:406),'r')
plot(spike_delay_t,getval(conv( sum(neutral_spikedata(end_delay-400-6:end_delay+3,:),2), 0.125*ones(1,8)),6:406),'g')
xlabel('Time before stimulus (ms)'); ylabel('mean spikes/bin(1ms)');
%legend('PrefDirec Cue','NullDirec Cue','Neutral Cue');

lfp_stim_t = [0:2:1000]; spike_stim_t = [0:1:1000];
subplot(4,1,3); hold on;
plot(lfp_stim_t,getval(conv( sqrt(mean(PD_lfpdata(start_stim-2:start_stim+500+1,:).^2,2)), 0.25*ones(1,4)),4:504),'b');
plot(lfp_stim_t,getval(conv( sqrt(mean(ND_lfpdata(start_stim-2:start_stim+500+1,:).^2,2)), 0.25*ones(1,4)),4:504),'r');
plot(lfp_stim_t,getval(conv( sqrt(mean(neutral_lfpdata(start_stim-2:start_stim+500+1,:).^2,2)), 0.25*ones(1,4)),4:504),'g');
%xlabel(FILE);
%xlabel('Time before stimulus (ms)'); 
ylabel('rms LFP'); 
%legend('PrefDirec Cue','NullDirec Cue','Neutral Cue'); legend(gca,'boxoff');

subplot(4,1,4); hold on;
plot(spike_stim_t,getval(conv( mean(PD_spikedata(end_delay-6:end_delay+1000+3,find(cue_val(PD_cue_trials)~=CUEONLY)),2), 0.125*ones(1,8)),6:1006),'b')
plot(spike_stim_t,getval(conv( mean(ND_spikedata(end_delay-6:end_delay+1000+3,find(cue_val(ND_cue_trials)~=CUEONLY)),2), 0.125*ones(1,8)),6:1006),'r')
plot(spike_stim_t,getval(conv( mean(neutral_spikedata(end_delay-6:end_delay+1000+3,:),2), 0.125*ones(1,8)),6:1006),'g') %not necessary to exclude CueOnly from Neutral trials
xlabel('Time after stimulus (ms)'); ylabel('mean spikes/bin(1ms)');
%legend('PrefDirec Cue','NullDirec Cue','Neutral Cue');


%keyboard
output = 1;
if (output)
%     %------------------------------------------------------------------------
%     %write out all relevant parameters to a cumulative text file, VR 11/21/05
%     outfile = [BASE_PATH 'ProtocolSpecific\CuedDirectionDiscrim\CueDir_Effects_summary.dat'];
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Tt:dl_lfp\t Tt:dl_bp\t Tt:dl_spk\t Tt:st_lfp\t Tt:st_bp\t Tt:st_spk\t 1A:dl_lfp\t 1A:dl_bp\t 1A:dl_spk\t 1A:st_lfp\t 1A:st_bp\t 1A:st_spk\t 2A-ME:st_lfp\t 2A-In:st_lfp\t 2A-ME:st_bp\t 2A-In:st_bp\t 2A-ME:st_spk\t 2A-In:st_spk\t 1A-Zs:st_lfp\t 1A-Zs:st_bp\t 1A-Zs:st_spk\t Mean_Dl_Lfp_Pd\t Std_Dl_Lfp_Pd\t Mean_Dl_Lfp_Neu\t Std_Dl_Lfp_Neu\t Mean_Dl_Lfp_Nd\t Std_Dl_Lfp_Nd\t Mean_Dl_Bp_Pd\t Std_Dl_Bp_Pd\t Mean_Dl_Lfp_Bu\t Std_Dl_Bp_Neu\t Mean_Dl_Bp_Nd\t Std_Dl_Bp_Nd\t Mean_Dl_Spk_Pd\t Std_Dl_Spk_Pd\t Mean_Dl_Spk_Neu\t Std_Dl_Spk_Neu\t Mean_Dl_Spk_Nd\t Std_Dl_Spk_Nd\t Mean_St_Lfp_Pd\t Std_St_Lfp_Pd\t Mean_St_Lfp_Neu\t Std_St_Lfp_Neu\t Mean_St_Lfp_Nd\t Std_St_Lfp_Nd\t Mean_St_Bp_Pd\t Std_St_Bp_Pd\t Mean_St_Bp_Neu\t Std_St_Bp_Neu\t Mean_St_Bp_Nd\t Std_St_Bp_Nd\t Mean_St_Spk_Pd\t Std_St_Spk_Pd\t Mean_St_Spk_Neu\t Std_St_Spk_Neu\t Mean_St_Spk_Nd\t Std_St_Spk_Nd\t');
%         
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
%     buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t',...
%         FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%         p_delaylfp, p_delaylfp_bp, p_delayspikes, p_stimlfp, p_stimlfp_bp, p_stimspikes, ...
%         p_1way_anova_delaylfp, p_1way_anova_delaylfp_bp, p_1way_anova_delayspikes, ...
%         p_1way_anova_stimlfp, p_1way_anova_stimlfp_bp, p_1way_anova_stimspikes, ...
%         p_2way_anova_stimlfp(1), p_2way_anova_stimlfp(3), p_2way_anova_stimlfp_bp(1), p_2way_anova_stimlfp_bp(3), ...
%         p_2way_anova_stimspikes(1), p_2way_anova_stimspikes(3), ...
%         p_zscore_anova_stimlfp, p_zscore_anova_stimlfp_bp, p_zscore_anova_stimspikes, ...
%         mean_delay_lfp(1), std_delay_lfp(1), mean_delay_lfp(2), std_delay_lfp(2), mean_delay_lfp(3), std_delay_lfp(3), ...
%         mean_delay_lfp_bp(1), std_delay_lfp_bp(1), mean_delay_lfp_bp(2), std_delay_lfp_bp(2), mean_delay_lfp_bp(3), std_delay_lfp_bp(3), ...
%         mean_delay_spikes(1), std_delay_spikes(1), mean_delay_spikes(2), std_delay_spikes(2), mean_delay_spikes(3), std_delay_spikes(3), ...
%         mean_stim_zlfp(1), std_stim_zlfp(1), mean_stim_zlfp(2), std_stim_zlfp(2), mean_stim_zlfp(3), std_stim_zlfp(3), ...
%         mean_stim_zlfp(1), std_stim_zlfp_bp(1), mean_stim_zlfp_bp(2), std_stim_zlfp_bp(2), mean_stim_zlfp_bp(3), std_stim_zlfp_bp(3), ...
%         mean_stim_zspikes(1), std_stim_zspikes(1), mean_stim_zspikes(2), std_stim_zspikes(2), mean_stim_zspikes(3), std_stim_zspikes(3) );
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
    %------------------------------------------------------------------------
    %write out LFP parameters to a cumulative text file, VR 11/21/05
    outfile = [BASE_PATH 'ProtocolSpecific\CuedDirectionDiscrim\CueDir_LFP_Effects_summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Tt:dl_lfp\t Tt:st_lfp\t 1A:dl_lfp\t 1A:st_lfp\t 2A-ME:st_lfp\t 2A-In:st_lfp\t 1A-Zs:st_lfp\t Mean_Dl_Lfp_Pd\t Std_Dl_Lfp_Pd\t Mean_Dl_Lfp_Neu\t Std_Dl_Lfp_Neu\t Mean_Dl_Lfp_Nd\t Std_Dl_Lfp_Nd\t Mean_St_Lfp_Pd\t Std_St_Lfp_Pd\t Mean_St_Lfp_Neu\t Std_St_Lfp_Neu\t Mean_St_Lfp_Nd\t Std_St_Lfp_Nd\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t',...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        p_delaylfp, p_stimlfp, p_1way_anova_delaylfp, p_1way_anova_stimlfp, ...
        p_2way_anova_stimlfp(1), p_2way_anova_stimlfp(3), p_zscore_anova_stimlfp, ...
        mean_delay_lfp(1), std_delay_lfp(1), mean_delay_lfp(2), std_delay_lfp(2), mean_delay_lfp(3), std_delay_lfp(3), ...
        mean_stim_zlfp(1), std_stim_zlfp(1), mean_stim_zlfp(2), std_stim_zlfp(2), mean_stim_zlfp(3), std_stim_zlfp(3));
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %----------------------------------------------------------------------
    %write out LFP Power parameters to a cumulative text file, VR 11/21/05
    outfile = [BASE_PATH 'ProtocolSpecific\CuedDirectionDiscrim\CueDir_LFP-BP_Effects_summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Tt:dl_bp\t Tt:st_bp\t 1A:dl_bp\t 1A:st_bp\t 2A-ME:st_bp\t 2A-In:st_bp\t 1A-Zs:st_bp\t Mean_Dl_Bp_Pd\t Std_Dl_Bp_Pd\t Mean_Dl_Bp_Neu\t Std_Dl_Bp_Neu\t Mean_Dl_Bp_Nd\t Std_Dl_Bp_Nd\t Mean_St_Bp_Pd\t Std_St_Bp_Pd\t Mean_St_Bp_Neu\t Std_St_Bp_Neu\t Mean_St_Bp_Nd\t Std_St_Bp_Nd\t');
        
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t',...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        p_delaylfp_bp, p_stimlfp_bp, p_1way_anova_delaylfp_bp, p_1way_anova_stimlfp_bp, ...
        p_2way_anova_stimlfp_bp(1), p_2way_anova_stimlfp_bp(3), p_zscore_anova_stimlfp_bp, ...
        mean_delay_lfp_bp(1), std_delay_lfp_bp(1), mean_delay_lfp_bp(2), std_delay_lfp_bp(2), mean_delay_lfp_bp(3), std_delay_lfp_bp(3), ...
        mean_stim_zlfp_bp(1), std_stim_zlfp_bp(1), mean_stim_zlfp_bp(2), std_stim_zlfp_bp(2), mean_stim_zlfp_bp(3), std_stim_zlfp_bp(3));
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %------------------------------------------------------------------------
    %write out spike-count parameters to a cumulative text file, VR 11/21/05
    outfile = [BASE_PATH 'ProtocolSpecific\CuedDirectionDiscrim\CueDir_Spikes_Effects_summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Tt:dl_spk\t Tt:st_spk\t 1A:dl_spk\t 1A:st_spk\t 2A-ME:st_spk\t 2A-In:st_spk\t 1A-Zs:st_spk\t Mean_Dl_Spk_Pd\t Std_Dl_Spk_Pd\t Mean_Dl_Spk_Neu\t Std_Dl_Spk_Neu\t Mean_Dl_Spk_Nd\t Std_Dl_Spk_Nd\t Mean_St_Spk_Pd\t Std_St_Spk_Pd\t Mean_St_Spk_Neu\t Std_St_Spk_Neu\t Mean_St_Spk_Nd\t Std_St_Spk_Nd\t');
        
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t',...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        p_delayspikes, p_stimspikes, p_1way_anova_delayspikes, p_1way_anova_stimspikes, ...
        p_2way_anova_stimspikes(1), p_2way_anova_stimspikes(3), p_zscore_anova_stimspikes, ...
        mean_delay_spikes(1), std_delay_spikes(1), mean_delay_spikes(2), std_delay_spikes(2), mean_delay_spikes(3), std_delay_spikes(3), ...
        mean_stim_zspikes(1), std_stim_zspikes(1), mean_stim_zspikes(2), std_stim_zspikes(2), mean_stim_zspikes(3), std_stim_zspikes(3) );
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
end
%% ********************** SAVE FIGS *****************************
SAVE_FIGS = 0;
if SAVE_FIGS
    saveas(hlist(1), sprintf('%s_CueDirEffects.fig',FILE),'fig');
    saveas(hlist(2), sprintf('%s_CueDirEffects_TimeC.fig',FILE),'fig');
end

return