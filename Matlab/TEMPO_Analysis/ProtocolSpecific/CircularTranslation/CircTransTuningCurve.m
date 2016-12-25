%-----------------------------------------------------------------------------------------------------------------------
%-- CircTransTuningCurve.m -- Plots out responses to CW and CCW motin and runs t-test on distribution of responses.
%--	BJAP, 7/20/00
%-----------------------------------------------------------------------------------------------------------------------
function CircTransTuningCurve(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE)

Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

%analyze the Fourier transforms of one cycle of the response to each
%stimulus.
N_trials = length(data.spike_rates(SpikeChan, :));
harmonics = ComputeFourierHarmonics(data, N_trials, StartCode, StopCode, StartOffset, StopOffset);

for cond = 1:size(unique_conds,1)
    select_trials = conditions(4,:) == cond - 1;
    spike_rates{cond} = data.spike_rates(SpikeChan, select_trials);
    first_harmonics{cond} = harmonics(SpikeChan, select_trials, 1);
end

%the next few lines of code for sorting out the spike rate distributions depend on consistent ordering of conditions
cw_conds = unique_conds(:,1) == 1;
ccw_conds = unique_conds(:, 1) == -1;
j = 1;
for i = 1:length(cw_conds)
    if cw_conds(i) ~= 0
        sorted_rates{j, 1} = spike_rates{i};
        sorted_first_harmonics{j, 1} = first_harmonics{i};
        sorted_cond(j,:) = unique_conds(i,2:3);
        j = j + 1;
    end
end
j = 1;
for i = 1:length(ccw_conds)
    if ccw_conds(i) ~= 0
        sorted_rates{j, 2} = spike_rates{i};
        sorted_first_harmonics{j, 2} = first_harmonics{i};
        j = j + 1;
    end
end

null_rate = mean(spike_rates{1});
null_first_harmonic = mean(first_harmonics{1});

for i=1:size(sorted_cond, 1);
    p_value(i) = signtest(sorted_rates{i, 1}, sorted_rates{i, 2}, 0.05 );
    first_harm_p_value(i) = signtest(sorted_first_harmonics{i, 1}, sorted_first_harmonics{i, 2}, 0.05 );
end

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Circular Translation');
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

axis([0 100 0 100]);
axis('off');
xpos = -10;
ypos = 20;
font_size = 9;
bump_size = 7;

for i = 1:size(sorted_cond,1)
    line = sprintf('%s %d %s %d: Signtest p = %4.3f (rates), = %4.3f (first harmonics)', param{2}, sorted_cond(i,1), param{3}, sorted_cond(i,2), p_value(i) , first_harm_p_value(i)  );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end    

for i = 1:size(sorted_cond,1)
    cw_spike_rates(i) = mean(sorted_rates{i,1} );
    cw_std_err(i) = std(sorted_rates{i,1} ) / sqrt( length(sorted_rates{1}) );
    ccw_spike_rates(i) = mean(sorted_rates{i,2} );
    ccw_std_err(i) = std(sorted_rates{i,2} ) / sqrt( length(sorted_rates{1}) );

    cw_first_harmonic(i) = mean(sorted_first_harmonics{i,1} );
    cw_harmonic_std_err(i) = std(sorted_first_harmonics{i,1} ) / sqrt( length(sorted_first_harmonics{1}) );
    ccw_first_harmonic(i) = mean(sorted_first_harmonics{i,2} );
    ccw_harmonic_std_err(i) = std(sorted_first_harmonics{i,2} ) / sqrt( length(sorted_first_harmonics{1}) );    
end

subplot(2, 2, 3);
plot(sorted_cond(:,1), cw_spike_rates, 'b-', sorted_cond(:,1), ccw_spike_rates, 'r-');
hold on;
errorbar(sorted_cond(:,1),  cw_spike_rates, cw_std_err, cw_std_err, 'bo');
hold on;
errorbar(sorted_cond(:,1),  ccw_spike_rates, ccw_std_err, ccw_std_err, 'ro');
hold on;
plot([min(sorted_cond(:,1) ) max(sorted_cond(:,1) )], [null_rate null_rate], 'k--');

xlabel('Background Color');
ylabel('Firing Rate (spk/s)');
xlim([min(sorted_cond(:,1) ) max(sorted_cond(:,1) ) ]);


subplot(2, 2, 4);
plot(sorted_cond(:,1), cw_first_harmonic, 'b-', sorted_cond(:,1), ccw_first_harmonic, 'r-');
hold on;
errorbar(sorted_cond(:,1),  cw_first_harmonic, cw_harmonic_std_err, cw_harmonic_std_err, 'bo');
hold on;
errorbar(sorted_cond(:,1),  ccw_first_harmonic, ccw_harmonic_std_err, ccw_harmonic_std_err, 'ro');
hold on;
plot([min(sorted_cond(:,1) ) max(sorted_cond(:,1) )], [null_first_harmonic null_first_harmonic], 'k--');

xlabel('Background Color');
ylabel('First Harmonic Amplitude (spk/s)');
xlim([min(sorted_cond(:,1) ) max(sorted_cond(:,1) ) ]);

