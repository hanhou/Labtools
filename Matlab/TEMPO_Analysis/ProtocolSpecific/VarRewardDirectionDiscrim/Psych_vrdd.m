%-----------------------------------------------------------------------------------------------------------------------
%-- Psych_vrdd.m -- Plots psychometric curve sorted by cue validity but collapsed across direction
%--	VR, 9/23/08
%-----------------------------------------------------------------------------------------------------------------------
function Psych_vrdd(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%% Set up some basic variables related to the task conditions and outcomes
TEMPO_Defs;
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,BegTrial:EndTrial,PATCH1);
unique_direction = munique(direction');
Pref_direction = data.one_time_params(PREFERRED_DIRECTION);
if unique_direction(1)~=Pref_direction  %this ensures that the preferred direction (T1) is the first element of unique_direction array.
    unique_direction = [unique_direction(2) unique_direction(1)];
end

%get the motion coherences
coherence = data.dots_params(DOTS_COHER, BegTrial:EndTrial, PATCH1);  %this is the coherence irrespective of the direction of motion
unique_coherence = munique(coherence');
%unique_coherence(isnan(unique_coherence)) = [];
signed_coherence = coherence.*(-1+2.*(direction==Pref_direction)); %this assigns coherence in the preferred (T1) direction to positive numbers.
%null (T2) direction motion coherences are assigned negative values
unique_signed_coherence = [-unique_coherence' unique_coherence'];

%get the reward ratios, as coded in the MICROSTIM variable.
reward_ratio = data.misc_params(MICROSTIM,BegTrial:EndTrial,PATCH1);
unique_reward_ratio = munique(reward_ratio');

%get outcome for each trial: 0=incorrect, 1=correct
trials_outcomes = logical(data.misc_params(OUTCOME,BegTrial:EndTrial) == CORRECT);

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(coherence);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

OPTIONS = OPTIMSET('MaxIter', 1000000,'MaxFunEvals',200000);  %this will be useful for fitting

%% Do analysis here
% You'll need to select the trials for every combination of reward ratio and signed_coherence and compute
% the fraction T1-choices.  Make a plot of the fraction T1-choices vs
% signed_coherences (just like what the online code does) and also compute
% fits for the data.  (Use logistic_fit.m to fit the data to a logistic function and look at the
% code to figure out which output is the bias parameter.)
% Then go ahead and do the permutation analysis, and any other analyses you
% choose to.  For now, just do this on a single saved data set... We can
% worry about combining data across datasets later.
% You may find it useful to look at the some of my (probably inelegant and poorly commented!) matlab code in the
% Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim
% directory.  Also, feel free to ask me any questions.  This might be a bit
% of a steep learning curve, that to some extent you won't get to use again since
% you're only doing the psychophsyics in TEMPO and Matlab... Nonetheless,
% it'll probably be a useful experience.  Good luck!
% P.S.  Oh yeah, to run the code, you can first load the gui by typing
% TEMPO_GUI, then navigate to the datafile of interest in the red box by
% typing the data directory (i.e., Z:\Data\Tempo\Baskin\Raw\) and selecting
% the datafile (or you can navigate it by clicking the Browse button).
% Then click the Load Data button and that should, after about 20 seconds
% load the data.  There is only one available analysis option (ie., Plot
% Psychometric).  Don't worry about anything else and click the green
% Analyze button to execute this file.
% And remember: breakpoints are your friends!
% also, to accomplish this you won't need to to edit any other files (and
% you probably won't have write access to do so....).  But read whatever
% code you want to understand better or to get ideas of how to accomplish
% this.

%make direction positive
ics = unique_direction < 0;
unique_direction(ics) = unique_direction(ics) + 360;
ics = direction < 0;
direction(ics) = direction(ics) + 360;

%split into two files for the two directions:
odirection = direction;
ocoherence = coherence;
oreward_ratio = reward_ratio;
oselect_trials = select_trials;
osigned_coherence = signed_coherence;
otrials = trials;
otrials_outcomes = trials_outcomes;

for u = 1 : length(unique_reward_ratio),
    %determine whether moving to high reward or low reward
%    signtowardT1 = double(unique_direction(u) == Pref_direction) - 0.5; %positive if the shown motion was meant to be toward T1
 %   ics = ((odirection == unique_direction(u)) & (signtowardT1 * reward_ratio > 0)), %pick trials for given directions where we find high reward at this direction
    %find indices for one direction
    ics = oreward_ratio * unique_reward_ratio(u) > 0; %will first find the indices where T1 has smaller reward, on second iteration, where T1 has the larger reward
    highrwdir = unique_direction(ceil(-sign(unique_reward_ratio(u)) / 2 + 1)); %reports the direction of the large reward
    
    direction = odirection(ics);
    coherence = ocoherence(ics);
    reward_ratio = oreward_ratio(ics);
    select_trials = oselect_trials(ics);
    signed_coherence = osigned_coherence(ics);
    trials = otrials(ics);
    trials_outcomes = otrials_outcomes(ics);

    %psychoreward(signed_coherence, unique_signed_coherence, reward_ratio, unique_reward_ratio, direction, Pref_direction, trials_outcomes, FILE);

    %create record
    load record;

    rc = length(record) + 1;
    datenum = datevec(date());
    datestr = sprintf('%04d-%02d-%02d', datenum(1), datenum(2), datenum(3));
    record(rc).builddate = datestr;
    record(rc).buildrunno = rc;
    
    record(rc).highrwdir = highrwdir;
    record(rc).rw = unique_reward_ratio(u);
    record(rc).eyewnd = sprintf('[%s]', num2str(data.targ_params(1:2,1,1)'));
    [dct1] = data.targ_params(1:2,1,2)-data.targ_params(1:2,1,1); ecc = round(sqrt(dct1' * dct1));
    record(rc).ecc = ecc;
    record(rc).rewrat = sprintf('[%d %d]', unique_reward_ratio(1), unique_reward_ratio(2));
    record(rc).coh = sprintf('[%d %d]', unique_coherence(1), unique_coherence(end));
    record(rc).weight = NaN;
    %find the date of the file by looking at the file properties
    cdr = cd;
    cd(PATH);
    a=dir;
    cd(cdr);
    acell = {a(:).name};
    fix = strmatch(FILE, acell);
    datestr = a(fix(end)).date;
    datenum = datevec(datestr);
    datestr = sprintf('%04d-%02d-%02d', datenum(1), datenum(2), datenum(3));
    record(rc).date = datestr;

    %save T1 label
    if highrwdir == Pref_direction,
        highrwdirisT1 = true;
    else
        highrwdirisT1 = false;
    end
    record(rc).Pref_direction = Pref_direction;
    record(rc).highrwdirisT1 = highrwdirisT1;
    
    %compute run number
    try
        prevdate = record(rc - 1).date;
        prevorigfile = record(rc - 1).origfile;
        run = record(rc - 1).run;
    catch
        prevdate = 'initDate';
        prevorigfile = 'initFile';
        run = 1;
    end

    if ~strcmp(record(rc).date, prevdate)
        run = 1;
    else
        if ~strcmp(record(rc).origfile, prevorigfile)
            run = run + 1;
        end
    end
    record(rc).run = run;
    
    %save file
    [path,mfile,a2,a3]=fileparts(FILE);
    record(rc).origfile = mfile;
    record(rc).origpath = PATH;
   
    path = ['C:\BaskinData\mfiles\'];
    record(rc).path = path;
    mfile = sprintf('%s_highrwdir=%d_ecc=%d', mfile, highrwdir, ecc);
    record(rc).filename = mfile;
        
    save([path, mfile]);
    
    %save record
    save record record;    
end