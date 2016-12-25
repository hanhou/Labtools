%----------------------------------------------------------------------------------------------------------
%-- LoadTEMPOData.m: This function opens the TEMPO HTB file and reads in all of the data and header information
%--	contained therein.  All data is returned in the large structure arrays, good_data and bad_data.  By doing it this way,
%-- 	I am keeping every scrap of data about every trial all together in these big structures.  Note that various
%--	defines in TEMPO_Defs.m are needed to access the data from good/bad_data.  GCD, starting 12/25/99
%----------------------------------------------------------------------------------------------------------
function	[return_value, good_data, bad_data, newChans] = LoadTEMPOData(PATH,FILE)

global listtext num_recorded_spike_channels

TEMPO_Defs;		%reads in some definitions that we need
Path_Defs;		%some path definitions

good_data = [];
bad_data = [];
newChans = [];

l = length(FILE);
if (FILE(l-3:l) == '.htb')	% .htb extension already there
    filename = [PATH FILE];   %the HTB data file
    logfile = [PATH FILE(1:l-4) '.log'];   %the TEMPO log file
else	%no extension in FILE, add extensions
    filename = [PATH FILE '.htb'];   %the HTB data file
    logfile = [PATH FILE '.log'];   %the TEMPO log file
    FILE = [FILE '.htb'];
end

%% Rescue available trial information from CED in cases where .htb/.log files are lost.
% Note that this needs to be careful and checked manually for each cell. HH20150624
% Also should add the entries in "BATCH_GUI_Tempo_Analysis.m" for overiding in Batch analysis
rescue_from_CED = {   % FileName   Protocol   Conditions   Headings  Coherence
    'Z:\Data\MOOG\Polo\raw\m5c77r2.htb','HD',[1 2 3],[-8 -4 -2 -1 0 1 2 4 8], 35;
    'Z:\Data\MOOG\Polo\raw\m5c90r2.htb','HD',[1 2 3],[-6 -3 -1.5 -0.75 0 0.75 1.5 3 6], 35;
    'Z:\Data\MOOG\Polo\raw\m5c91r1.htb','MemSac', nan, 0:45:350, nan;
    'Z:\Data\MOOG\Polo\raw\m5c118r3.htb','HD',[1 2 3],[-6 -3 -1.5 -0.75 0 0.75 1.5 3 6], 30;
    'Z:\Data\MOOG\Messi\raw\m10c50r3.htb','HD',[1 2 3],[-8 -4 -2 -1 0 1 2 4 8], 15;    
    'Z:\Data\MOOG\Messi\raw\m10c70r3.htb','HD',[1 2 3],[-8 -4 -2 -1 0 1 2 4 8], 16;
    'Z:\Data\MOOG\Messi\raw\m10c102r5.htb','DelSac',nan, 0:45:350, nan;
    'Z:\Data\MOOG\Messi\raw\m10c104r9.htb','HD',[1 2 3], [-8 -4 -2 -1 0 1 2 4 8], 14;
    'Z:\Data\MOOG\Messi\raw\m10c168r3.htb','HD',[1 2 3], [-8 -4 -2 -1 0 1 2 4 8], 8;
    };

ff = 1;
while ff <= size(rescue_from_CED,1) && isempty(strfind(rescue_from_CED{ff,1},[PATH FILE]))
    ff = ff + 1;
end

if ff <= size(rescue_from_CED,1) % Cell that needs to be rescued
    beep; pause(0.1); beep;
    fprintf('=== Rescuing data from CED ... \n');
    rescue_info = rescue_from_CED(ff,2:end);
    good_data = RescueFromCED(PATH,FILE,rescue_info);
    fprintf('=========  Success !! =========\n');

else
    %% Read htb START (normal condition)
    
    %disp('opening  HTB file...');
    fid = HTBOPEN(filename);            % Open the HTB file
    if (fid == -1)		%file could not be opened
        return_value = -1;
        return;
    end
    
    %disp('counting databases...');
    ndbs = HTBCOUNT(fid);               % Get number of databases in it
    
    %read in data and headers from the TEMPO .htb file
    %identify type of each database and store it in appropriate array
    temp1 = []; temp2 = []; temp3 = []; temp4 = [];
    for database = 1:ndbs               % Loop though each one...
        %store the database headers in good_data and in bad_data
        temp_hd = HTBGETHD(fid, database);
        
        if strcmp(temp_hd.title, 'Eye Traces')
            good_data.htb_header{EYE_DB} = temp_hd;
            bad_data.htb_header{EYE_DB} = temp_hd;
            %disp('reading eye data...');
            temp1 = HTBGETDA(good_data.htb_header{EYE_DB});  %gets data for *all* epochs into a single matrix with dimentions: (sample#*trial#, channel#)
            temp1 = reshape(temp1', [good_data.htb_header{EYE_DB}.nchannels, good_data.htb_header{EYE_DB}.period, good_data.htb_header{EYE_DB}.sweep]);
            %*NOTE*:after this manipulation, eye_data has the following indices: (channel #, sample #, trial #).  I am doing this because it is 10x faster
            %than reading in the data for each trial using htbGetEp()
        end
        
        if strcmp(temp_hd.title, 'Spikes')
            good_data.htb_header{SPIKE_DB} = temp_hd;
            bad_data.htb_header{SPIKE_DB} = temp_hd;
            %disp('reading spike data...');
            temp2 = HTBGETDA(good_data.htb_header{SPIKE_DB});
            %                 temp2 = zeros(1000, 5);
            if size(temp2) == [good_data.htb_header{SPIKE_DB}.period * good_data.htb_header{SPIKE_DB}.sweep good_data.htb_header{SPIKE_DB}.nchannels]
                temp2 = reshape(temp2', [good_data.htb_header{SPIKE_DB}.nchannels, good_data.htb_header{SPIKE_DB}.period, good_data.htb_header{SPIKE_DB}.sweep]);
            else              % *********************************************************************************************************
                temp2 = [];   % CRF -- 3/8/06 -- TEMPORARY -- error in Tempo screwed up the spike DB in some of my data files (discrim_2I)
            end               % *********************************************************************************************************
            num_recorded_spike_channels = size(temp2, 1);
        end
        
        if strcmp(temp_hd.title, 'Events')
            good_data.htb_header{EVENT_DB} = temp_hd;
            bad_data.htb_header{EVENT_DB} = temp_hd;
            %disp('reading event data...');
            temp3 = HTBGETDA(good_data.htb_header{EVENT_DB});
            temp3 = reshape(temp3', [good_data.htb_header{EVENT_DB}.nchannels, good_data.htb_header{EVENT_DB}.period, good_data.htb_header{EVENT_DB}.sweep]);
        end
        
        if ( strcmp(temp_hd.title, 'Local Field Potentials') | strcmp(temp_hd.title, 'LFP') )
            good_data.htb_header{LFP_DB} = temp_hd;
            bad_data.htb_header{LFP_DB} = temp_hd;
            disp('reading LFP data...');
            temp4 = HTBGETDA(good_data.htb_header{LFP_DB});
            temp4 = reshape(temp4', [good_data.htb_header{LFP_DB}.nchannels, good_data.htb_header{LFP_DB}.period, good_data.htb_header{LFP_DB}.sweep]);
        end
        
        
    end
    
    err = HTBCLOSE(fid);                % Close HTB file
    
    %determine number of trials so that we can declare arrays in ReadTEMPOLog
    num_trials = good_data.htb_header{EVENT_DB}.sweep;
    %now, read in parameters from the TEMPO log file, which contains all experimental parameters for each trial
    [revcorr_params, dots_params, moog_params, targ_params, cue_params, misc_params, one_time_params, eye_calib_params, obj_params, bar_params, bkgnd_params, neuron_params] = ReadTEMPOLog(logfile, num_trials);
    %dots_params[] has indices (param#, trial#, patch#)
    %moog_params[] has indices (param#, trial#, item#)
    %targ_params[] has indices (param#, trial#, targ#)
    %cue_params[] has indices (param#, trial#, cue#)
    %misc_params[] has indices (param#, trial#)
    %one_time_params[] has index (param#)
    %eye_calib_params[] has indices (param#, value#)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % added by BJP 1/3/01 to incorporate binding data
    %obj_params[] has indices (param#, trial#, obj#)
    %bar_params[] has indices (param#, trial#, obj#, bar#)
    %bkgnd_params[] has indices (param#, trial#)
    %neuron_params[] has indices (param#, #spike channels)
    
    %now, remove any trials that were not completed successfully as we store the data in good_data
    %and store data for trials not completed in bad_data
    if isempty(temp2)
        all_trials = 1:size(temp3,3);	%list of indices for all trials; CRF, 2/21/06 -- added IF statement because some of my behavior-only (discrim) datasets have no spike DB (and thus temp2 cannot be used to index all_trials)
    else
        all_trials = 1:size(temp2,3);	%list of indices for all trials; GCD changed this 5/18/04 to deal with a strange data file that had extra trial in event channel; should be OK this way
    end
    good_trials = find(temp3(:,:,:) == SUCCESS_CD);
    fprintf('Success/Total Trials = %g / %g\n',length(good_trials),num_trials);  % HH20130825
    
    good_trials = ceil(good_trials/good_data.htb_header{EVENT_DB}.period);  %these are now trial indices
    %NOTE: don't use htb_header.sweep hereafter for the number of trials, since this includes trials that were
    %not completed successfully (e.g. breaks of fixation).  GCD, 12/28/99
    all_trials(good_trials) = NaN;	%mark the good trials with NaNs
    bad_trials = all_trials(~isnan(all_trials));  %and then the bad trials are the members of all_trials ~= NaN
    
    % Override the default eye channel settings. HH20150722
    if ~isnan(one_time_params(LEFT_EYE_X_CHANNEL))
        LEYE_H = one_time_params(LEFT_EYE_X_CHANNEL);
    end
    
    if ~isnan(one_time_params(LEFT_EYE_Y_CHANNEL))
        LEYE_V = one_time_params(LEFT_EYE_Y_CHANNEL);
    end
    
    if ~isnan(one_time_params(RIGHT_EYE_X_CHANNEL))
        REYE_H = one_time_params(RIGHT_EYE_X_CHANNEL);
    end
    
    if ~isnan(one_time_params(RIGHT_EYE_Y_CHANNEL))
        REYE_V = one_time_params(RIGHT_EYE_Y_CHANNEL);
    end
    
    %now, I am going to stash the data into the appropriate
    %parts of my big data structures, called good_data and bad_data
    if isempty(temp1)
        good_data.eye_data = [];
        bad_data.eye_data = [];
    else
        %Convert the eye_data here from A/D units to degrees of visual angle
        %This is a change in the approach implemented to make software calibration easier.  GCD 12/28/00
        temp1(LEYE_H, :, :) = temp1(LEYE_H, :, :) .* (one_time_params(X_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
        temp1(LEYE_V, :, :) = temp1(LEYE_V, :, :) .* (one_time_params(Y_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
        temp1(REYE_H, :, :) = temp1(REYE_H, :, :) .* (one_time_params(X_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
        temp1(REYE_V, :, :) = temp1(REYE_V, :, :) .* (one_time_params(Y_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
        
        % if there are D/A channels in addition to eye channels (not overlapped, HH20150722)
        if (size(temp1,1) > 4) && (isempty(intersect([LEYE_H,LEYE_V,REYE_H,REYE_V],[DA_H DA_V])))  
            temp1(DA_H, :, :) = temp1(DA_H, :, :) .* (one_time_params(X_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
            temp1(DA_V, :, :) = temp1(DA_V, :, :) .* (one_time_params(Y_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
        end
        
        %Do the software calibration of eye position signals here, if software calibration was used in collection
        %Here, I operate directly on the temp1 array, replacing the previous values
        %Added by GCD, 12/22/00
        if (one_time_params(SOFTWARE_CALIB_STATUS) == 1)
            junk = temp1;
            if isnan(eye_calib_params(SOFT_CAL_LEYE_HORIZ,4)) %if only 3 calibration parameters
                temp1(LEYE_H, :, :) = eye_calib_params(SOFT_CAL_LEYE_HORIZ,1) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,2).*junk(LEYE_H,:,:) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,3).*junk(LEYE_V,:,:);
                temp1(LEYE_V, :, :) = eye_calib_params(SOFT_CAL_LEYE_VERT,1) + eye_calib_params(SOFT_CAL_LEYE_VERT,2).*junk(LEYE_V,:,:) + eye_calib_params(SOFT_CAL_LEYE_VERT,3).*junk(LEYE_H,:,:);
                temp1(REYE_H, :, :) = eye_calib_params(SOFT_CAL_REYE_HORIZ,1) + eye_calib_params(SOFT_CAL_REYE_HORIZ,2).*junk(REYE_H,:,:) + eye_calib_params(SOFT_CAL_REYE_HORIZ,3).*junk(REYE_V,:,:);
                temp1(REYE_V, :, :) = eye_calib_params(SOFT_CAL_REYE_VERT,1) + eye_calib_params(SOFT_CAL_REYE_VERT,2).*junk(REYE_V,:,:) + eye_calib_params(SOFT_CAL_REYE_VERT,3).*junk(REYE_H,:,:);
            else %if all 4 parameters are used
                temp1(LEYE_H, :, :) = eye_calib_params(SOFT_CAL_LEYE_HORIZ,1) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,2).*junk(LEYE_H,:,:) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,3).*junk(LEYE_V,:,:) + eye_calib_params(SOFT_CAL_LEYE_HORIZ,4).*junk(LEYE_V,:,:).*junk(LEYE_H,:,:);
                temp1(LEYE_V, :, :) = eye_calib_params(SOFT_CAL_LEYE_VERT,1) + eye_calib_params(SOFT_CAL_LEYE_VERT,2).*junk(LEYE_V,:,:) + eye_calib_params(SOFT_CAL_LEYE_VERT,3).*junk(LEYE_H,:,:) + eye_calib_params(SOFT_CAL_LEYE_VERT,4).*junk(LEYE_H,:,:).*junk(LEYE_V,:,:);
                temp1(REYE_H, :, :) = eye_calib_params(SOFT_CAL_REYE_HORIZ,1) + eye_calib_params(SOFT_CAL_REYE_HORIZ,2).*junk(REYE_H,:,:) + eye_calib_params(SOFT_CAL_REYE_HORIZ,3).*junk(REYE_V,:,:) + eye_calib_params(SOFT_CAL_REYE_HORIZ,4).*junk(REYE_V,:,:).*junk(REYE_H,:,:);
                temp1(REYE_V, :, :) = eye_calib_params(SOFT_CAL_REYE_VERT,1) + eye_calib_params(SOFT_CAL_REYE_VERT,2).*junk(REYE_V,:,:) + eye_calib_params(SOFT_CAL_REYE_VERT,3).*junk(REYE_H,:,:) + eye_calib_params(SOFT_CAL_REYE_VERT,4).*junk(REYE_H,:,:).*junk(REYE_V,:,:);
            end
            ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
            listtext{length(listtext)+1} = 'Using software-calibrated eye data.';
            set(ListHandle, 'String', listtext);
        end
        good_data.eye_data = temp1(:,:,good_trials);
        bad_data.eye_data = temp1(:,:,bad_trials);
    end
    clear temp1;
    
    if isempty(temp2)
        good_data.spike_data = [];
        bad_data.spike_data = [];
    else
        %default
        good_data.spike_data =  temp2(:,:,good_trials);
        bad_data.spike_data = temp2(:,:,bad_trials);
        
        %%%Yong's preference
        %     good_data.spike_data =  uint8(temp2(:,:,good_trials));
        %     bad_data.spike_data = uint8(temp2(:,:,bad_trials));
    end
    clear temp2;
    
    if isempty(temp3)
        good_data.event_data = [];
        bad_data.event_data = [];
    else
        good_data.event_data = temp3(:,:,good_trials);
        bad_data.event_data = temp3(:,:,bad_trials);
    end
    clear temp3;
    
    if isempty(temp4)
        good_data.lfp_data = [];
        bad_data.lfp_data = [];
    else
        for i = 1:size(temp4,1)
            temp4(i, :, :) = temp4(i, :, :)/one_time_params(AD_RANGE);
        end
        good_data.lfp_data = temp4(:,:,good_trials);
        bad_data.lfp_data = temp4(:,:,bad_trials);
    end
    clear temp4;
    
    good_data.targ_params = targ_params(:,good_trials,:);
    good_data.misc_params = misc_params(:,good_trials,:);
    good_data.one_time_params = one_time_params;
    good_data.eye_calib_params = eye_calib_params;
    good_data.neuron_params = neuron_params;		%neuron params array is used for both binding and dots protocols
    
    if (~isempty(revcorr_params))
        good_data.revcorr_params = revcorr_params(:, good_trials, :);
    end
    if (~isempty(moog_params))
        good_data.moog_params = moog_params(:, good_trials, :);
    end
    if ~isempty(dots_params)
        good_data.dots_params = dots_params(:,good_trials,:);
    end
    if ~isempty(cue_params)
        good_data.cue_params = cue_params(:,good_trials,:);
    end
    if (~isempty(obj_params))
        good_data.obj_params = obj_params(:,good_trials,:);
    end
    if (~isempty(bar_params))
        good_data.bar_params = bar_params(:,good_trials,:,:);
    end
    if (~isempty(bkgnd_params))
        good_data.bkgnd_params = bkgnd_params(:,good_trials,:);
    end
    %else%if ~isempty(moog_params)
    %good_data.moog_params = moog_params(:,good_trials,:);
    %else
    %binding protocol values
    %    good_data.obj_params = obj_params(:,good_trials,:);
    %    good_data.bar_params = bar_params(:,good_trials,:,:);
    %   good_data.bkgnd_params = bkgnd_params(:,good_trials,:);
    %end
    
    bad_data.targ_params = targ_params(:,bad_trials,:);
    bad_data.misc_params = misc_params(:,bad_trials,:);
    bad_data.one_time_params = one_time_params;
    bad_data.eye_calib_params = eye_calib_params;
    bad_data.neuron_params = neuron_params;
    
    if (~isempty(revcorr_params))
        bad_data.revcorr_params = revcorr_params(:, bad_trials, :);
    end
    if ~isempty(cue_params)
        bad_data.cue_params = cue_params(:,bad_trials,:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % added by BJP 1/3/01 for binding stuff
    if ~isempty(dots_params)
        bad_data.dots_params = dots_params(:,bad_trials,:);
    elseif ~isempty(moog_params)
        bad_data.moog_params = moog_params(:,bad_trials,:);
    else
        bad_data.obj_params = obj_params(:,bad_trials,:);
        bad_data.bar_params = bar_params(:,bad_trials,:,:);
        bad_data.bkgnd_params = bkgnd_params(:,bad_trials,:);
    end
    
    % Check to see if the htb file was in a standard location.
    PATH2 = upper(PATH);
    wasinstdloc = findstr(PATH2, 'Z:\DATA\TEMPO');
    
    % This adds the good data from spikesort2.
    fa = find(PATH == '\');
    
    % Prevent and indexing error if the htb file was in a nonstandard location.
    if (~isempty(wasinstdloc))
        fName = ['Z:\Data\CED', PATH(fa(3):fa(4)), 'SortedSpikes\', FILE(1:length(FILE) - 4)];
    else
        if findstr(PATH, 'Z:\Data\MOOG')%Add by AHC
            fName = ['Z:\Data\MOOG', PATH(fa(3):fa(4)), 'CED\SortedSpikes\', FILE(1:length(FILE) - 4)]; %Add by AHC
            good_data.fName=fName;
        else
            fName = [];
        end
    end
    
    
    %%%% Read htb END %%%%%
    
    
    check_spont = 0;
    global LOAD_MU_DATA;
    if (LOAD_MU_DATA & exist([fName, '.mat'], 'file'))
        disp('Running Packaroni, Booyah!!!');
        [good_data winDiscrim] = Packaroni(good_data, fName, 1, 0);  %get window discriminator data
        [good_data shiftValue]= SpikeCorr(good_data, 1, 1, winDiscrim);
        if  ( (shiftValue == -50) & (sum( sum( good_data.spike_data(winDiscrim,:,:) ) == 0) ) )
            % remove window discrimination data if not being used
            good_data.spike_data = good_data.spike_data(1:size(good_data.spike_data,1) - 1,:,:);
            shiftValue = 0;
        end
        [good_data newChans] = Packaroni(good_data, fName, 0, shiftValue);  %get rest of data
        for i = 1:length(newChans)
            [good_data shiftValue]= SpikeCorr(good_data, 0, 1, newChans(i));
            if check_spont == 1
                num_temp = size(good_data.event_data);
                num_t = num_temp(3);
                spont = 0;
                for j=1:num_t
                    start_code_bin = find(good_data.event_data(1, :, j) == 3);
                    stop_code_bin = find(good_data.event_data(1, :, j) == 4);
                    spikes = sum(good_data.spike_data(newChans(i), start_code_bin:stop_code_bin, j));
                    time = stop_code_bin - start_code_bin;
                    spont = spikes/time + spont;
                end
                spont_avg = spont/j * 1000
            end
        end
    end
end

%% Load CED spike sorting data %%

return_value = 1;		%indicates completed OK

%add by AHC 02-22-06
fa = find(PATH == '\');
sName = ['Z:\Data\MOOG', PATH(fa(3):fa(4)), 'Analysis\SortedSpikes2\', FILE(1:length(FILE) - 4),'.mat'];
if (exist(sName,'file'))
    [good_data SortedChannel] = PackData(good_data,sName);%get the Spike2 sorted data
    if SortedChannel >0
        disp('Loading Spike2 sorted data!');
    elseif SortedChannel == -1 % Not match htb & spike2
        return_value = -999;  % Not match htb & spike2
        %       error('Htb & spike2 files not matched ...');  % I throw an error message only when we are in batch mode
    end
end

%  added2loadtempodata

% Check spike bins. HH20150909
abnormal_spike_bins = find(good_data.spike_data > 1);
if ~isempty(abnormal_spike_bins)
    [abnormal_chans,~,~] = ind2sub(size(good_data.spike_data),abnormal_spike_bins);
    unique_abnormal_chans = unique(abnormal_chans);
    fprintf('\n\n***** WARNING: Abnormal spike bins detected *****\n');
    for aa = 1:length(unique_abnormal_chans)
        if unique_abnormal_chans(aa) >= 5 && unique_abnormal_chans(aa) < 20
            beep; 
        end        
        abnormal_values = good_data.spike_data(abnormal_spike_bins(abnormal_chans == unique_abnormal_chans(aa)));
        unique_abnormal_values = unique(abnormal_values);
        fprintf('Chan %d:\n',unique_abnormal_chans(aa));
        for uu = 1:length(unique_abnormal_values)
            fprintf('\tSpike count = %d, %d time(s)\n',unique_abnormal_values(uu),sum(abnormal_values == unique_abnormal_values(uu)));    
        end
    end
    
    good_data.spike_data(abnormal_spike_bins) = 1;
    
    fprintf('******* All abnormal values were set to 1 *******\n');
end

return;