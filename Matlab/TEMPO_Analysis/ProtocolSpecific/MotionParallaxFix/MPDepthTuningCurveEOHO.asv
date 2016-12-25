%-----------------------------------------------------------------------------------------------------------------------
%-- MPDepthTuningCurveEOHO.m -- Plots tuning curves for depth from motion parallax.
%-- Kluge for EOHO conditions separates this from MPDepthTuningCurve.m
%-- No longer can calculate PDIs.
%-- Started by JWN, 12/16/05
%-- Last by JWN, 04/10/06
%-----------------------------------------------------------------------------------------------------------------------
function MPDepthTuningCurveEOHO(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

%TEMPO_Defs;
ProtocolDefs;
Path_Defs;

symbols = {'bo' 'rs' 'gd' 'kv' 'm<' 'c>' 'bv' 'rv'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'm<-' 'c>-' 'yo-' 'bs-'};
NULL_VALUE = -9999;

disp(sprintf('(MPDepthTuningCurveEOHO) Started at %s.',datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
% Get the mean firing rates for all the trials

area = 'MT';  % Kluge! 80 for MT and 80 for MST (see Kruse et al 2002), +80 for transfer function delay
if(strcmp(area,'MT'))  % Don't change this one!
    latency = 160;  % MT guess
else
    latency = 160;  % MST guess
end 
begin_time = find(data.event_data(1,:,1)==StartCode) + latency;
% end_time = find(data.event_data(1,:,1)==StopCode) + latency;
end_time = begin_time + 1999; % 2s trial
if(max(max(max(data.spike_data))) > 1)
    data.spike_data = cast(data.spike_data>0,'double');
end
raw_spikes = data.spike_data(1,begin_time:end_time,:);
spike_rates = 1000*squeeze(mean(raw_spikes))';  % The hard way

% Calculate number of trials
trials = size(MPphase,2);

set(1, 'HandleVisibility', 'off');
figure;  set(gcf, 'HandleVisibility', 'on');
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573],'Name', 'Depth Tuning Curves');  % Better for printing?
final_data = zeros(10,28); % To be written out to files
final_data(:,1)=1:10;

%%%  Calculate PDI and PDImod  %%%
% Now also calculate iPDIs (PDIs for the four individual conditions)
% mis = MPGetMI(latency, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, FILE);
if(size(unique(MPtrial_types),2)==2) % Only for EOHO sets
%     for i = 1:2 % no true PDI (MP-RM) calced for EOHO 
%         if i>2 % true PDI
%             for j = 1:num_depths-1
%                 mean_dataMP(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == 0)';
%                 mean_dataRM(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == 2)'; 
%             end
%             mean_data = mean_dataMP - mean_dataRM;
%             misMP = squeeze(mis(1,:,:));
%             misRM = squeeze(mis(3,:,:));
%             mis_data = misMP-misRM;
%         else  % iPDIs
%             for j = 1:num_depths-1
%                 mean_data(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == i+4-1)'; %+4 to shift for EOHO
%             end
%             mis_data = squeeze(mis(i,:,:));
%         end
%         meanmean_data = mean(mean_data);  % mean of the trials from the means measure
%         meanmis = mean(mis_data,2);  % mean of the trials from the mis measure
%         stdmean_data = std(mean_data);
%         stdmis = std(mis_data,0,2);
%         for k = 1:num_depths/2-1
%             nearm = meanmean_data(k);
%             farm = meanmean_data(num_depths-k);
%             nearstd = stdmean_data(k);
%             farstd = stdmean_data(num_depths-k);
%             PDI(k) = (farm-nearm)/(abs(farm-nearm)+sqrt((nearstd^2+farstd^2)/2));
%             nearm = meanmis(k+1);  % Remember to shift by one to lose spont data
%             farm = meanmis(num_depths-(k-1));  % and by six to get to the far data
%             nearstd = stdmis(k+1);
%             farstd = stdmis(num_depths-(k-1));
%             PDImod(k) = (farm-nearm)/(abs(farm-nearm)+sqrt((nearstd^2+farstd^2)/2));
%         end
%         final_data(i,22:25) = PDI;
%         final_data(i+7,22:25) = PDImod;
%         mPDI(i) = mean(PDI);
%         mPDImod(i) = mean(PDImod);
%         final_data(i,26) = mPDI(i);
%         final_data(i+7,26) = mPDImod(i);
%         % Send data off to MPBootstrap for significance testing.
%         pmPDI(i) = MPBootstrap(mean_data, mPDI(i));
%         pmPDImod(i) = MPBootstrap(mis_data(2:end,:)', mPDImod(i));
%         final_data(i,27:28) = [pmPDI(i) (pmPDI(i)>0.025)+1];
%         final_data(i+7,27:28) = [pmPDImod(i) (pmPDImod(i)>0.025)+1];
%     end
else  % if not full set, just zero things out so you don't bug the figure
    PDI = zeros(1,4);
    PDImod = zeros(1,4);
    mPDI = zeros(1,7);
    mPDImod = zeros(1,7);
end    

clf(gcf);
% Plot tuning curves with error bars
titles = { 'Motion Parallax', 'Binocular Disparity', 'Retinal Motion', 'Congruent', 'Eye Movement Only', 'Head Movement Only' };
MItitles = {};  % Fill in once we see what conditions we have
for i=1:6  % Six blocks, only one spike channel
    subplot(2, 1, 1);	% Use a different subplot for means vs. MIs.
    hold on;
    indices = logical((MPtrial_types == (i-1)) & (MPdepths ~= NULL_VALUE));  % Both phases in this version
    if(sum(indices)>0)
        MItitles = [MItitles titles(i)];
        [unique_x, mean_rate, std_err, spk_max, spk_min] = PlotTuningCurve(MPdepths(indices)', spike_rates(indices)', symbols{i}, line_types4{i},1,1);
        indices = logical((MPtrial_types == (i-1)) & (MPdepths == NULL_VALUE));
        final_data(i+0,2) = mean(spike_rates(indices)); % 2 is null, i+0 for means
        final_data(i+0,3:11) = mean_rate'; % Tuning curve data points
        final_data(i+0,12) = std(spike_rates(indices))/sqrt(size(spike_rates(indices),2)); % 12 is null err, i+0 for means
        final_data(i+0,13:21) = std_err'; % Tuning curve error bars
    end
end
% Do MP-RM (i=7) final_data stuffing here to keep near other final_data; not plotted
for j = 1:num_depths
    tempMP(:,j) = spike_rates(MPdepths == uMPdepths(j) & MPtrial_types == 0)';
    tempRM(:,j) = spike_rates(MPdepths == uMPdepths(j) & MPtrial_types == 2)'; 
end
% subtracted_rates = tempMP - tempRM;
% final_data(7+0,2:11) = mean(subtracted_rates); % Tuning curve data points including null
% final_data(7+0,12:21) = std(subtracted_rates)/sqrt(size(subtracted_rates,2)); % Tuning curve error bars including null

h = axis;
h(1:2) = [-2.2 2.2];
axis(h);
YLabel('Response (spikes/s)');
XLabel('Depth (deg)');
RFinfo = sprintf('(%2.1f,%2.1f) %2.0f', data.neuron_params(RF_XCTR), data.neuron_params(RF_YCTR), data.neuron_params(RF_DIAMETER));
PDIinfo = '';
iPDIinfo = '';
%PDIinfo = sprintf('%1.2f ',PDI);
%iPDIinfo = sprintf('%1.2f ',mPDI(1:2)); % Only 2 for EOHO
Title(sprintf('%s  %s  %s\n%s\niPDIs=[ %s]', area, RFinfo, FILE,'Means',iPDIinfo));  % Simplified for EOHO

% % Get MIs (we already have from above PDI calcs) and plot them
% % mis = MPGetMI(latency, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin);
% meanmis = mean(mis,3);
% errmis = std(mis,0,3)/sqrt(size(mis,3));
% % if(size(meanmis)==[6 10]) % Only for full sets, pulled out for EO/HO
% %     subtracted_mis = (misMP-misRM)';  % misMP and misRM are calculated way above
% %     final_data(6:9,2:11) = meanmis; % Tuning curve data points including null
% %     final_data(6:9,12:21) = errmis; % Tuning curve error bars including null
% %     final_data(7+5,2:11) = mean(subtracted_mis); % Tuning curve data points including null
% %     final_data(7+5,12:21) = std(subtracted_mis)/sqrt(size(subtracted_mis,2)); % Tuning curve error bars including null
% % end
% subplot(2, 1, 2);  % Get subplot in main figure window
% hold on;
% for i=1:size(mis,1)
%     switch(char(MItitles(i)));  % Legend kluge
%         case 'Motion Parallax', color = 1;
%         case 'Binocular Disparity', color = 2;
%         case 'Retinal Motion', color = 3;
%         case 'Congruent', color = 4;
%         case 'Eye Movement Only', color = 5;
%         case 'Head Movement Only', color = 6;
%     end
%     errorbar(uMPdepths(find(uMPdepths~=NULL_VALUE)),meanmis(i,2:num_depths),errmis(i,2:num_depths),line_types5{color});
%     anova_data = squeeze(mis(i,:,:))';
% end
% children = get(gca,'children');
% Legend(children(end:-1:1),MItitles);
% Legend(gca,'boxoff')
% hold off;
% XLabel('Depth (deg)');
% YLabel('Amplitude (spikes/s)');
% % PDIinfo = sprintf('%1.2f ',PDImod);
% % iPDIinfo = sprintf('%1.2f ',mPDImod(1:2)); % Only 2 for EOHO
% Title(sprintf('%s\niPDIs=[ %s]', 'Modulation Indices',iPDIinfo));% Simplified for EOHO
% h = axis;
% h(1:2) = [-2.2 2.2];
% axis(h);
%print(2); % Uncomment for printing.
%close(gcf); % Uncomment for printing.

%Write results for this cell to 13 files (MP,BD,RM,C,EO,HO,MP-RM + means,MIs == 1,2,3,4,5,6,7 + 0,7 and then an ALL)
PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
filenames = {'MP_meansEOHO','BD_meansEOHO','RM_meansEOHO','C_meansEOHO','EO_meansEOHO','HO_meansEOHO','MP-RM_meansEOHO','MP_MIsEOHO','BD_MIsEOHO','RM_MIsEOHO','C_MIsEOHO','MP-RM_MIsEOHO','ALLEOHO'};
for i = [1 2 3 4 5 6]
    outfile = cell2mat(strcat(PATHOUT,area,'_',filenames(i),'.txt'));
    headerflag = 0;
    if (exist(outfile) == 0) % File does not yet exist, so print a header
        headerflag = 1;
    end
    fid = fopen(outfile, 'a');  % Open text file.
    if(i<15)
        if (headerflag)
            fprintf(fid, 'FILE Condition null neg2.0 neg1.5 neg1.0 neg0.5 0.0 0.5 1.0 1.5 2.0 errnull errneg2.0 errneg1.5 errneg1.0 errneg0.5 err0.0 err0.5 err1.0 err1.5 err2.0 PDI2.0 PDI1.5 PDI1.0 PDI0.5 mPDI p sigp');
            fprintf(fid, '\r\n');
        end
        fprintf(fid,'%10s', strtok(FILE,'.'));
        fprintf(fid,' %+2.5f', final_data(i,:));
    else
        if (headerflag)
            fprintf(fid, 'FILE MPiPDI BDiPDI RMiPDI CiPDI PDI absMPiPDI absBDiPDI absRMiPDI absCiPDI absPDI sigpMPiPDI sigpBDiPDI sigpRMiPDI sigpCiPDI sigpPDI MPiPDImod BDiPDImod RMiPDImod CiPDImod PDImod absMPiPDImod absBDiPDImod absRMiPDImod absCiPDImod absPDImod sigpMPiPDImod sigpBDiPDImod sigpRMiPDImod sigpCiPDImod sigpPDImod');
            fprintf(fid, '\r\n');
        end
        fprintf(fid,'%10s', strtok(FILE,'.'));
        all_output = [ final_data(1:5,26)' abs(final_data(1:5,26))' final_data(1:5,28)' final_data(6:10,26)' abs(final_data(6:10,26))' final_data(6:10,28)'];
        fprintf(fid,' %+2.4f', all_output);
    end    
    fprintf(fid,'\r\n');
    fclose(fid);
end
disp('(MPDepthTuningCurveEOHO) Done.');

% MP_PSTH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

return;