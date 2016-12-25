%-----------------------------------------------------------------------------------------------------------------------
%-- TimeBlockAnalysis.m -- Analyze CPs, etc over blocks of time within a discrimination run
%--	GCD, 3/1/01
%-----------------------------------------------------------------------------------------------------------------------
function TimeBlockAnalysis(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

num_trials = EndTrial - BegTrial + 1;
BlockSize = num_trials/8;
block_starts = floor(BegTrial:BlockSize:(EndTrial-BlockSize+1))
block_ends = floor((BegTrial+BlockSize-1):BlockSize:EndTrial)

grandCP = []; grandPVal = [];
NeuroThresh = []; MonkeyThresh = [];
for i=1:length(block_starts)
    [gcp, Pval] = Compute_ChoiceProb(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, block_starts(i), block_ends(i), StartOffset, StopOffset, PATH, FILE);
    grandCP(i) = gcp;
    grandPVal(i) = Pval;
    close 2 3 4
    
    [Nth, Mth] = NeuroPsychoPlot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, block_starts(i), block_ends(i), StartOffset, StopOffset, PATH, FILE);
    NeuroThresh(i) = Nth;
    MonkeyThresh(i) = Mth;
    close 2
end

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 50 400 573], 'Name', 'Time Block Analysis');
subplot(2,1,1);
plot((block_starts+block_ends)/2, grandCP, 'ko-', 'MarkerFaceColor', 'k');
hold on;
chance = ones(length(block_starts))*0.5;
plot((block_starts+block_ends)/2, chance, 'k--');
text(block_starts, (grandCP-0.06), num2str(grandPVal'));
xlim([0 EndTrial]);
ylim([0 1]);
hold off;
titl = [PATH FILE];
title(titl);
xlabel('Trial #');
ylabel('Grand Choice Prob.');

subplot(2,1,2);
hold on;
Handl(1) = plot((block_starts+block_ends)/2, NeuroThresh, 'ro-', 'MarkerFaceColor', 'r');
Handl(2) = plot((block_starts+block_ends)/2, MonkeyThresh, 'go-', 'MarkerFaceColor', 'g');
xlim([0 EndTrial]);
xlabel('Trial #');
ylabel('Threshold (% dots)');

legend(Handl, 'Neuron', 'Monkey', 2);

return;