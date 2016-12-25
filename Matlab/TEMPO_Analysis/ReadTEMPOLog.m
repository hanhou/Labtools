%----------------------------------------------------------------------------------------------------------
%-- ReadTEMPOLog.m: This function reads in all the data from the TEMPO *.log file, which contains a list
%--	of experiment parameters for each trial.	GCD, 12/29/99
%----------------------------------------------------------------------------------------------------------
function [revcorr_params, dots_params, moog_params, targ_params, cue_params, misc_params, one_time_params, eye_calib_params, obj_params, bar_params, bkgnd_params, neuron_params] = ReadTEMPOLog(logfile, num_trials);

TEMPO_Defs;	%some defines that we'll need

%% main bottleneck is reading and parsing data out into these two cell arnrays - BJP 5/2/01
logfile
[keys, data] = textread(logfile, '%s %[^\n]', 'bufsize', 300000);
num_lines_read = size(keys, 1);

%now, until we reach the end of the file, do the following:
%	get a line, pull off the keyword, and act on the data according to what the keyword is
prev_trial = 0;

%initialize arrays and preallocate space to optimize speed - BJP 5/2/01
%preallocation speeds processing up significantly
dots_params = ones(NUM_DOTS_PARAMS, num_trials, NUM_PATCHES)*NaN; 
moog_params = ones(NUM_MOOG_PARAMS, num_trials, NUM_MOOG_ITEMS)*NaN;
%revcorr_params = zeros(NUM_REVCORR_PARAMS, num_trials, numRevcorrParams)*NaN;
misc_params = ones(NUM_MISC_PARAMS, num_trials, NUM_CUES)*NaN;
cue_params = ones(NUM_CUE_PARAMS, num_trials, NUM_CUES)*NaN;
targ_params = ones(NUM_TARG_PARAMS, num_trials, NUM_TARGETS)*NaN;
obj_params = ones(NUM_OBJ_PARAMS, num_trials, NUM_OBJECTS)*NaN;
bar_params = ones(NUM_BAR_PARAMS, num_trials, NUM_OBJECTS, NUM_BARS)*NaN;
bkgnd_params = ones(NUM_BKGND_PARAMS, num_trials)*NaN;
one_time_params = ones(NUM_ONE_TIME_PARAMS, 1)*NaN;
eye_calib_params = ones(4,4)*NaN; %now compatible with 4 calibration params. VR 7/6/06
dots_protocol = 0;
moog_protocol = 0;
revcorr_protocol = 0;

initRevCorrParams = 1;
revcorr_params = [];

for i=1:num_lines_read
    switch keys{i}
    case 'TRIAL#'
        curr_trial = sscanf(data{i},'%d');  
        if ( (curr_trial - prev_trial) > 1)
            disp('WARNING: There is a missing Trial # in the TEMPO LOG file!!');   
            line = sprintf('prev_trial = %d, curr_trial = %d', prev_trial, curr_trial);
            disp(line);
            beep;
        end
        prev_trial = curr_trial;   
    case dots_keywords
        index = eval(keys{i});
        patch_values = sscanf(data{i}, '%f');
        if(size(patch_values,1) ~= size(dots_params,3))  % JWN 01/31/05
            disp(sprintf('(ReadTEMPOLog) WARNING: Dimension mismatch.  Ignoring %s.',keys{i}));
        else
            dots_params(index, curr_trial, :) = patch_values;
        end
        dots_protocol = 1;		%set flag if dots protocol
    case moog_keywords
        index = eval(keys{i});
        moog_values = sscanf(data{i}, '%f');
        moog_params(index, curr_trial, 1:length(moog_values)) = moog_values;
        moog_protocol = 1;		%set flag if moog protocol
    case revcorr_keywords
        index = eval(keys{i});

        revcorr_values = sscanf(data{i}, '%f');
        
        % Create a big set of revcorr_params if we're loading a REVCORR
        % protocol.  Otherwise, make a small set.
        if initRevCorrParams
            if one_time_params(PROTOCOL) == DIRECTION_REVCORR
                revcorr_params = zeros(NUM_REVCORR_PARAMS, num_trials, 1500)*NaN;
            else
                revcorr_params = zeros(NUM_REVCORR_PARAMS, num_trials, 5)*NaN;
            end
            initRevCorrParams = 0;
        end
        
        % Put a dummy value in revcorr_values if the keyword has no data
        % associated with it.
        if (isempty(revcorr_values))
            revcorr_values = 0;
        end
        
        revcorr_params(index, curr_trial, 1:length(revcorr_values)) = revcorr_values;
        revcorr_protocol = 1;
    case targ_keywords
        index = eval(keys{i});
        targ_values = sscanf(data{i}, '%f');
        % I had to adapt the following lines to account for older data files that 
        % had three targets rather than five. - BJP 11/2/01
        targ_params(index, curr_trial, 1:length(targ_values) ) = targ_values;
    case cue_keywords
        index = eval(keys{i});
        cue_values = sscanf(data{i}, '%f');
        cue_params(index, curr_trial, 1:length(cue_values) ) = cue_values;
    case misc_keywords
        index = eval(keys{i});
        misc_values = sscanf(data{i}, '%f');
        misc_params(index, curr_trial, :) = misc_values;      
	case neuron_keywords
        index = eval(keys{i});
        neuron_values = sscanf(data{i}, '%f');
        neuron_params(index, :) = neuron_values';
        
        %stash copy of params for first neuron into one_time params db
        %clean up protocol specific analysis tools later
        one_time_params(index, :) = neuron_values(1);
    case one_time_keywords
        index = eval(keys{i});
        one_time_values = sscanf(data{i}, '%f');
        one_time_params(index, :) = one_time_values;
    case eye_calib_keywords
        index = eval(keys{i});
        eye_calib_values = sscanf(data{i}, '%f');
        eye_calib_params(index, 1:length(eye_calib_values)) = eye_calib_values'; %now compatible with either 3 or 4 
                                                                                 %calibration params. VR 7/8/06
    case obj_keywords
  	  	  index = eval(keys{i});
        obj_values = sscanf(data{i}, '%f');
        obj_params(index, curr_trial, :) = obj_values';
    case bar_keywords
        %curr_trial
        index = eval(keys{i});
        bar_values = sscanf(data{i}, '%f');
        % right now index obj 1 = 0 BJP 1/3/01
        obj_num = bar_values(1) + 1;
        bar_values = bar_values(2:end);
        bar_params(index, curr_trial, obj_num, :) = bar_values';
    case bkgnd_keywords
        index = eval(keys{i});
        bkgnd_values = sscanf(data{i}, '%f');
        bkgnd_params(index, curr_trial, :) = bkgnd_values'; 
    otherwise
        disp(sprintf('(ReadTEMPOLog) WARNING: Unknown keyword %s.',keys{i}));  % JWN 01/31/05
    end
end

%empty appropriate matrices depending on protocol
%empty matrices used in LoadTEMPOData - BJP 5/2/01
% Added &'s to the ifs for stuff coming off the downstairs MOOG, which is both dots and moog.  JWN 01/31/05
if (dots_protocol & ~moog_protocol) % dots protocol 
   moog_params = [];
   obj_params = [];
   bar_params = [];
   bkgnd_params = [];
   %revcorr_params = [];
elseif (moog_protocol & ~dots_protocol) % moog protocol
   dots_params = [];
   obj_params = [];
   bar_params = [];
   bkgnd_params = [];
   %revcorr_params = [];
   cue_params = [];
elseif (~moog_protocol & ~dots_protocol)	%binding protocol
   dots_params = [];
   moog_params = [];
   %revcorr_params = [];
   cue_params = [];
end

return;