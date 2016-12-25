%-----------------------------------------------------------------------------------------------------------------------
%-- MPDepthTuningCurve.m -- Plots tuning curves for depth from motion parallax. 
%-- Started by JWN, 9/17/04
%-- Last by JWN, 4/13/05
%-----------------------------------------------------------------------------------------------------------------------
function MPDepthTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

line_types = {'bo-' 'ro-' 'go-' 'ko-' 'mo-' 'co-' 'yo-' 'bs-'};
symbols = {'bo' 'ro' 'bs' 'rs' 'bd' 'rd' 'bv' 'rv'};
colors = {'b' 'r' 'g' 'k' 'm' 'm' 'y' 'b'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types3 = {'k:' 'r:' 'g:' 'b:' 'g^-' 'b^-' 'r-^' 'k-^'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'mo-' 'co-' 'yo-' 'bs-'};
NULL_VALUE = -9999;

disp(sprintf('(MPDepthTuningCurve) Started at %s.',datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
% Get the mean firing rates for all the trials 
% spike_rates = data.spike_rates(SpikeChan, :);  % The easy way
latency = 100;
begin_time = find(data.event_data(1,:,1)==StartCode) + latency;
end_time = find(data.event_data(1,:,1)==StopCode) + latency;
raw_spikes = data.spike_data(1,begin_time:end_time,:);
spike_rates = 1000*squeeze(mean(raw_spikes))';  % The hard way

% Calculate number of trials
trials = size(MPphase,2);

set(1, 'HandleVisibility', 'off');
figure(2);  set(2, 'HandleVisibility', 'on');  clf(2);
set(2,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573],'Name', 'Depth Tuning Curves');  % Better for printing?

% Print a header
note = MPGetNote(FILE,'z:\Users\Jacob\BarracudaNotes.txt');
note = sprintf('(%d,%d) %d  %s', data.neuron_params(RF_XCTR), data.neuron_params(RF_YCTR), data.neuron_params(RF_DIAMETER), note);

clf(2);
% Plot tuning curves with error bars
titles = { sprintf('%s\n%s',note,'Motion Parallax'), 'Binocular Disparity', 'Retinal Motion', 'Congruent' };
for i=1:4  % Four blocks, only one spike channel
    subplot(5, 1, i);	% Use a different subplot for each different block.  Leave fifth subplot for MIs.
    hold on;
    indices = logical((MPtrial_types == (i-1)) & (MPdepths ~= NULL_VALUE) & MPphase == 0);  % Zero phase in blue
    if(sum(indices)>0)
        PlotTuningCurve(MPdepths(indices)', spike_rates(indices)', symbols{i*2-1}, line_types4{1},1,1);
    end
    hold on;
    indices = logical((MPtrial_types == (i-1)) & (MPdepths ~= NULL_VALUE) & MPphase == 180);  % 180 phase in red
    if(sum(indices)>0)
        PlotTuningCurve(MPdepths(indices)', spike_rates(indices)', symbols{i*2}, line_types4{2},1,1);
        % Reorganize spike_rates into depth columns for ANOVAing
        anova_data0 = zeros(sum(indices)/(num_depths-1),num_depths-1);
        anova_data180 = zeros(sum(indices)/(num_depths-1),num_depths-1);
        for j = 1:num_depths-1
            anova_data0(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == (i-1) & MPphase == 0)';    
            anova_data180(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == (i-1) & MPphase == 180)';
        end
        % Reorganize spike_rates into table for 2-way ANOVAing
        if(num_depths ~= 10)
            disp('(MPDepthTuningCurve) WARNING: 2-way ANOVA expects 9(+1) depths');
        else
            anova_dataNeg = [ anova_data0(:,4); anova_data180(:,4); anova_data0(:,3); anova_data180(:,3); anova_data0(:,2); anova_data180(:,2); anova_data0(:,1); anova_data180(:,1) ];   
            anova_dataPos = [ anova_data0(:,6); anova_data180(:,6); anova_data0(:,7); anova_data180(:,7); anova_data0(:,8); anova_data180(:,8); anova_data0(:,9); anova_data180(:,9) ];   
            anova_data2way = [ anova_dataNeg anova_dataPos ];
            psmi = anova2(anova_data2way,sum(indices)*2/(num_depths-1),'off');
            psign = psmi(1);
            pmag = psmi(2);
            pint = psmi(3);
        end
        p0 = anova1(anova_data0,([]),'off');
        if p0<0.05
            if p0<0.01
                if p0<0.001
                    stars = '***';
                else
                    stars = '**';
                end
            else
                stars = '*';
            end
        else
            stars = '';
        end
        p0 = sprintf('P = %0.4f %s',p0,stars);
        p180 = anova1(anova_data180,([]),'off');
        if p180<0.05
            if p180<0.01
                if p180<0.001
                    stars = '***';
                else
                    stars = '**';
                end
            else
                stars = '*';
            end
        else
            stars = '';
        end
        p180 = sprintf('P = %0.4f %s',p180,stars);
        children = get(gca,'children');
        Legend([children(5);children(2)],p0,p180);
        Legend(gca,'boxoff')
    end
end
YLabel('Response (spikes/sec)');
% Now add a dashed line for the spontaneous level
for i=1:4
    subplot(5, 1, i);	% Use a different subplot for each different block.  Leave fifth subplot for MIs.
    hold on;
    indices = logical((MPtrial_types == (i-1)) & (MPdepths == NULL_VALUE) & MPphase == 0);  % Zero phase in blue
    spont_x = [min(MPdepths(find(MPdepths~=-9999))) max(MPdepths(find(MPdepths~=-9999)))];
    spont_y = [mean(spike_rates(indices)) mean(spike_rates(indices))];
    plot(spont_x, spont_y,  line_types2{1});
    indices = logical((MPtrial_types == (i-1)) & (MPdepths == NULL_VALUE) & MPphase == 180);  % 180 phase in red
    spont_y = [mean(spike_rates(indices)) mean(spike_rates(indices))];
    plot(spont_x, spont_y,  line_types2{2});
    hold off;
    yl = YLim;
    YLim([0 yl(2)]);	% Set the lower limit of the Y axis to zero
    Title(titles{i});
end

% Get MIs and plot them
mis = MPGetMI(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin);
meanmis = mean(mis,3);
errmis = std(mis,0,3)/sqrt(size(mis,3));
subplot(5, 1, 5);  % Get subplot in main figure window
hold on;
for i=1:size(mis,1)
    errorbar(uMPdepths(find(uMPdepths~=NULL_VALUE)),meanmis(i,2:num_depths),errmis(i,2:num_depths),line_types5{i});
    anova_data = squeeze(mis(i,:,:))';
    p = anova1(anova_data,([]),'off');
    if p<0.05
        if p<0.01
            if p<0.001
                stars = '***';
            else
                stars = '**';
            end
        else
            stars = '*';
        end
    else
        stars = '';
    end
    ps{i} = sprintf('P = %0.4f %s',p,stars);
end
children = get(gca,'children');
switch(size(mis,1))
    case 4, Legend([children(7);children(5);children(3);children(1)],ps);
    case 3, Legend([children(5);children(3);children(1)],ps);
    case 2, Legend([children(3);children(1)],ps);
    case 1, Legend([children(1)],ps);
end
Legend(gca,'boxoff')
hold off;
XLabel('Depth (deg)');
YLabel('Amplitude (spikes/s)');
Title('Modulation Indices');
disp('(MPDepthTuningCurve) Done.');
return;