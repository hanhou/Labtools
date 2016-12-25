%-----------------------------------------------------------------------------------------------------------------------
%-- MP_PSTH -- Plots PSTHs. 
%-- Started by JWN 02/11/05
%-- Last by JWN 10/19/07
%-----------------------------------------------------------------------------------------------------------------------

function MP_PSTH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

% TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

line_types = {'bo-' 'ro-' 'go-' 'ko-' 'mo-' 'co-' 'yo-' 'bs-'};
symbols = {'bo' 'ro' 'go' 'ko' 'mo' 'co' 'yo' 'bs'};
colors = {'b' 'r' 'g' 'k' 'm' 'm' 'y' 'b'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types3 = {'k:' 'r:' 'g:' 'b:' 'g^-' 'b^-' 'r-^' 'k-^'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
NULL_VALUE = -9999;

disp(sprintf('(MP_PSTH) Started at %s.',datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);

if(max(max(max(data.spike_data))) > 1)
    disp(sprintf('(MP_PSTH) WARNING: %d corrupt values in data.spike_data.',sum(sum(sum(data.spike_data>1)))));
    data.spike_data = cast(data.spike_data>0,'double');
end

% Extract the bin indices corresponding to the codes
area = 'MST';  % Kluge! 80 for MT and 80 for MST (see Kruse et al 2002), +80 for transfer function delay
if(strcmp(area,'MT'))  % Don't change this one!
    latency = 80;  % MT guess  only 80 not 160 because not including neural latency
else
    latency = 80;  % MST guess
end
begin_time = find(data.event_data(1,:,1)==StartCode) + latency;
end_time = begin_time + 1999; % 2s trial
pre_time = 100;  post_time = 100;
begin_time = begin_time - pre_time;
end_time = end_time + post_time;
total_spike_bins = end_time - begin_time;

num_reduced_bins = 43;
bin_width = total_spike_bins/(num_reduced_bins+1);  % ~2200ms/(43+1) = ~50ms;
align_time = begin_time + pre_time;
marker_time = end_time - post_time;

uMPtrial_types = unique(MPtrial_types);
uMPdepths = unique(MPdepths);
uMPphase = unique(MPphase);
num_trial_types = size(uMPtrial_types,2);
num_depths = size(uMPdepths,2);
num_phase = size(uMPphase,2);
count_data = zeros(num_trial_types,num_depths,num_phase,num_reduced_bins+1);

% Loop through to extract histogram data and get upper bound for Y-axis
for i=1:num_trial_types  % Four MPtrial_type blocks (0=MP,1=BD,2=RM,3=C)
    for j=1:num_depths  % Ten MPdepths (which includes null)
        for k=1:num_phase    % Two phases (0 and 180)
            indices = logical((MPtrial_types == uMPtrial_types(i)) & (MPdepths == uMPdepths(j)) & (MPphase == uMPphase(k)));
            raw_spikes = data.spike_data(1,begin_time:end_time,indices);
            hist_data = sum(raw_spikes,3);
            [bins, counts] = SpikeBinner(hist_data, 1, bin_width, pre_time);
            count_data(i,j,k,:) = counts*(1000/bin_width)/sum(indices);  % Convert to instantaneous firing rates
        end
    end
end       
ymax = max(max(max(max(count_data))));
xx = [0 0];
xx2 = [marker_time-align_time marker_time-align_time];
yy = [0 ymax];
    
% Loop through for plotting
titles = { 'Motion Parallax', 'Binocular Disparity', 'Retinal Motion', 'Congruent', 'Eye Movement Only', 'Head Movement Only' };
% for i=1:num_trial_types  % Six MPtrial_type blocks (0=MP,1=BD,2=RM,3=C,4=EO,5=HO)
for i=[1 3 5 6]  % Just MP and RM and EOHO
%     MPhandle = figure(i+2);
    MPhandle = figure;
    clf(MPhandle);
    set(MPhandle,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', titles{uMPtrial_types(i)+1});
    MPsubplot = 0;
    for j=1:num_depths  % Ten MPdepths (which includes null)
        for k=1:num_phase    % Two phases (0 and 180)
            MPsubplot = MPsubplot+1;
            subplot(num_depths,num_phase,MPsubplot);
            hold on;
            % Plot counts as instantaneous firing rates
            bar(bins,squeeze(count_data(i,j,k,:)),1,colors{k});
            plot(xx, yy, 'k--');
            plot(xx2, yy,'k--');
            axis([0-pre_time,2000+post_time,0,ymax]);
            set(gca,'XTick',[0 500 1000 1500 2000])
            if(MPsubplot==1)
                YLabel('Null');
                Title(sprintf('%s\n %s ',FILE,titles{uMPtrial_types(i)+1}));
            end
        end
    end
%    print(MPhandle); % Uncomment for printing.
%    close(MPhandle); % Uncomment for printing.
end          

% Write results for this cell to a file
if(num_trial_types==6) % Only for full sets
    PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\PSTH\';
    filenames = {'PSTH_MP', 'PSTH_BD', 'PSTH_RM', 'PSTH_C', 'PSTH_EO', 'PSTH_HO'};
    for i = 1:6
        outfile = cell2mat(strcat(PATHOUT,strtok(FILE,'.'),'_',filenames(i),'.txt'));
        fid = fopen(outfile, 'w');  % Open text file.
        fprintf(fid, 'Condition Bins Null aNull neg2.0 aneg2.0 neg1.5 aneg1.5 neg1.0 aneg1.0 neg0.5 aneg0.5 0.0 a0.0 0.5 a0.5 1.0 a1.0 1.5 a1.5 2.0 a2.0');
        fprintf(fid, '\r\n');
        final_data = [zeros(size(count_data,4),1)+i bins squeeze(count_data(i,1,:,:))' squeeze(count_data(i,2,:,:))' squeeze(count_data(i,3,:,:))' squeeze(count_data(i,4,:,:))' squeeze(count_data(i,5,:,:))' squeeze(count_data(i,6,:,:))' squeeze(count_data(i,7,:,:))' squeeze(count_data(i,8,:,:))' squeeze(count_data(i,9,:,:))' squeeze(count_data(i,10,:,:))'];
        for j = 1:size(count_data,4)
            fprintf(fid,' %+2.4f', final_data(j,:));
            fprintf(fid,'\r\n');
        end
        fclose(fid);
    end
end

disp('(MP_PSTH) Done.');
return;