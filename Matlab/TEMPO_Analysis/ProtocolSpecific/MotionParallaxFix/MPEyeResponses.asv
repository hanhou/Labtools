%-----------------------------------------------------------------------------------------------------------------------
%-- MPEyeReponses.m -- Looks at the effects of eye position/velocity on responses.  Looks at the temporal correlations 
% between eye velocity differences and response differences using a range of possible response latencies.  Has not 
% been run in a long time, and is no longer up to date.   
%-- Started by JWN, 5/17/05
%-- Last by JWN, 05/03/07
%-----------------------------------------------------------------------------------------------------------------------
function MPEyeResponses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

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

disp(sprintf('(MPEyeResponses) Started at %s.',datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
uMPtrial_types = unique(MPtrial_types);  % Conditions present
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);

% % In data.eye_data, Channels 3&4 are eye (x&y), 5&6 are Moog (x&y).
% % Only analyze during visual stim + 200
pre_time = 0;  post_time = 200;
begin_time = find(data.event_data(1,:,1)==StartCode) - pre_time;
end_time = find(data.event_data(1,:,1)==StopCode) + post_time;
total_spike_bins = end_time - begin_time;
% eye_xy = data.eye_data(3:4,200:642,:);
% Moog_xy = data.eye_data(5:6,200:642,:);
% % Rescale Moog analog signals into degrees
% % MoogX = ((amoogx/nSets) * WIN_XRANGE/high_res)/AVAL_RANG
% % MoogZ = ((amoogz/nSets) * WIN_YRANGE/high_res)/AVAL_RANGE;
% % writef("c:\\TEMPO\\MoogProtocol\\params.log AD_RANGE %8d\n", AVAL_RANGE);	//number of levels in the A/D conversion (e.g., 12 bits -> 4096)
% % float high_res = 10.0;	//multiplies the _zoom variables to have high-res animated graphs
% % hide float constant WIN_XRANGE = 100.0*high_res;	//full scale will be this many degrees
% % hide float constant WIN_YRANGE = 100.0*high_res;
% % nSets = 10?
% 
% % Moog_xy = Moog_xy / 10 * 1000 / 65536;  We don't need to do this now for some reason?  JWN 100305
% 
% 
% % Realign axes to match preferred direction
% pref = data.neuron_params(PREFERRED_DIRECTION);
% opp = tan(pref/(180/pi));
% u = [1 opp] / sqrt(1+opp^2);
% v = [-u(2) u(1)];
% for i=1:size(eye_xy,3)
%     eye_uv(1,:,i) = u*eye_xy(:,:,i);
%     eye_uv(2,:,i) = v*eye_xy(:,:,i);
%     Moog_uv(1,:,i) = u*Moog_xy(:,:,i);
%     Moog_uv(2,:,i) = v*Moog_xy(:,:,i);
% end
% 
% % Velocity in deg/s
% veye_uv = diff(eye_uv(1,:,:))*200;
% vMoog_uv = diff(Moog_uv(1,:,:))*200;

eye_xyl = data.eye_data(1:2,200:642,:);
eye_xyr = data.eye_data(3:4,200:642,:);
Moog_xy = data.eye_data(5:6,200:642,:);
% Realign axes to match preferred direction
pref = data.neuron_params(PREFERRED_DIRECTION);
opp = tan(pref/(180/pi));
u = [1 opp] / sqrt(1+opp^2);
v = [-u(2) u(1)];
for i=1:size(eye_xyl,3)
    eye_uvr(1,:,i) = u*eye_xyr(:,:,i);
    eye_uvr(2,:,i) = v*eye_xyr(:,:,i);
    eye_uvl(1,:,i) = u*eye_xyl(:,:,i);
    eye_uvl(2,:,i) = v*eye_xyl(:,:,i);
    Moog_uv(1,:,i) = u*Moog_xy(:,:,i);
    Moog_uv(2,:,i) = v*Moog_xy(:,:,i);
end
eye_uv = (eye_uvl+eye_uvr)/2;  % Average of the two eyes
% Eye check based on file name (set up for Barracuda and Ovid 050107)
[monkid, cellid, runstr]=strread(FILE,'m%dc%dr%s.htb');
switch monkid
    case 9, % Barracuda
        if cellid==155, eye_uv = eye_uvl; end
        if cellid==156, eye_uv = eye_uvl; end
    case 15, % Ovid
        if cellid<35, eye_uv = eye_uvl; end
        if cellid>46, eye_uv = eye_uvr; end
    otherwise
        disp('(MPEyeResponses) WARNING: Unknown monkid');
end
% Velocity in deg/s
%veye_uvr = diff(eye_uvr(1,:,:))*200;
%veye_uvl = diff(eye_uvl(1,:,:))*200;
vMoog_uv = diff(Moog_uv(1,:,:))*200;
veye_uv = diff(eye_uv(1,:,:))*200;

% Defines for use in eye response binning (5ms bins to match velocity data)
num_reduced_bins = 441;
bin_width = total_spike_bins/(num_reduced_bins+1);  % ~2200ms/(442+1) = ~5ms;
boxcarwidth = 50;  % Used for later boxcar filtering with boxcarfilter(vector, width)
saccade_threshold = 100;  % deg/s to be considered a saccade-like noise

% For shifted correlation coefficients
shifts = 20:30:200; % 7 shifts are in ms
shifts = 20:10:200; % 19 shifts are in ms

% set(1, 'HandleVisibility', 'off');
% figure(2);  set(2, 'HandleVisibility', 'on');  clf(2);
% hold on;
% set(2,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573],'Name', 'Cross Correlograms');  % Better for printing?

% Generate a header
% note = MPGetNote(FILE,'z:\Users\Jacob\BarracudaNotes.txt');
% note = sprintf('(%d,%d) %d  %s', data.neuron_params(RF_XCTR), data.neuron_params(RF_YCTR), data.neuron_params(RF_DIAMETER), FILE);

% For the motion parallax condition, compare a 0 phase near with a 180
% phase far (and vice versa).
if(isempty(find(uMPtrial_types==0))) disp('(MPEyeResponses) Breakout: No MP');  return;  end;
i = 1;  % Only look at MP condition
for j = 1:num_depths  % Look at all 10 depths...
    if(uMPdepths(j)==0 | uMPdepths(j)==-9999) continue; end  % but throw away two of them (0 and null)
    indices_0 = logical((MPtrial_types == (i-1)) & MPphase == 0 & MPdepths == uMPdepths(j));
    trials_0 = find(indices_0==1);  % Get zero phase trials for this depth
    indices_180 = logical((MPtrial_types == (i-1)) & MPphase == 180 & MPdepths == -uMPdepths(j));
    trials_180 = find(indices_180==1);  % Get 180 phase trials for the -matched depth
%     % FOR TRIAL PAIR ANALYSIS (not currently used)
%     for k = 1:size(trials_0,2)  % Should have 3-5 trial pairs
%         % For each trial pair get eye data (velocity)
%         % delta_eyev is a subtraction after a sign flip
%         delta_veye = boxcarfilter(veye_uv(1,:,trials_0(k)),boxcarwidth) + boxcarfilter(veye_uv(1,:,trials_180(k)),boxcarwidth);
%         % For each trial pair get response data
%         % Responses, binned to match velocity data (taken at 200Hz)
%         hist_data = data.spike_data(1,begin_time:end_time,trials_0(k));
%         % hist_data = sum(raw_spikes,3);  % Not needed for single trial binning
%         [bins, counts0] = SpikeBinner(hist_data, 1, bin_width, pre_time);
%         counts0 = counts0*(1000/bin_width);  % Convert to instantaneous firing rates
%         % Again, for 180 this time
%         hist_data = data.spike_data(1,begin_time:end_time,trials_180(k));
%         % hist_data = sum(raw_spikes,3);  % Not needed for single trial binning
%         [bins, counts180] = SpikeBinner(hist_data, 1, bin_width, pre_time);
%         counts180 = counts180*(1000/bin_width);  % Convert to instantaneous firing rates
%         % Compute the differences in response, and smooth with a boxcar
%         delta_resp = boxcarfilter((counts0-counts180)',boxcarwidth);
%         
%         % Compute the cross-correlation function estimates for this trial pair with xcorr()
%         cc = xcorr(delta_veye,delta_resp,40);  % Limit the lag shifts to max of 200ms
%         cc = cc(41:81);  % Clip off negative lags
%         max_cc(j,k) = max(cc);  % Max coeffiecient
%         max_ct(j,k) = (find(cc==max(cc))-1)*5;  % Time of max in ms
%     end
    % FOR MEAN ANALYSIS
    m_veye_uv_0 = mean(veye_uv(1,:,indices_0),3); % Get zero phase mean for this depth
    m_veye_uv_180 = mean(veye_uv(1,:,indices_180),3);  % Get 180 phase mean for the -matched depth
    delta_veye = boxcarfilter(m_veye_uv_0,boxcarwidth) + boxcarfilter(m_veye_uv_180,boxcarwidth);
    raw_spikes = data.spike_data(1,begin_time:end_time,indices_0);
    hist_data = sum(raw_spikes,3);
    [bins, counts0] = SpikeBinner(hist_data, 1, bin_width, pre_time);
    counts0 = counts0*(1000/bin_width)/sum(indices_0);  % Convert to instantaneous firing rates
    raw_spikes = data.spike_data(1,begin_time:end_time,indices_180);
    hist_data = sum(raw_spikes,3);
    [bins, counts180] = SpikeBinner(hist_data, 1, bin_width, pre_time);
    counts180 = counts180*(1000/bin_width)/sum(indices_180);  % Convert to instantaneous firing rates
    % Compute the differences in response, and smooth with a boxcar
    delta_resp = boxcarfilter((counts0-counts180)',boxcarwidth);
   
    % Compute the cross-correlation function estimates for this trial pair with xcorr()
    ccraw = xcorr(delta_veye,delta_resp,40);  % Limit the lag shifts to max of 200ms
    cc(j,:) = ccraw(41:81);  % Clip off negative lags
    max_cc(j) = max(cc(j,:));  % Max coeffiecient
    max_ct(j) = (find(cc(j,:)==max(cc(j,:)))-1)*5;  % Time of max in ms
    
    % Compute single correlation coefficients for a few shifted times
    for k = 1:size(shifts,2)
      lag = shifts(k)/5;  % Shift is in ms, divide by 5 to put it in index units 
      tmp = corrcoef(delta_veye(1+lag:end),delta_resp(1:end-lag));
      corrcoefs(j,k) = tmp(2);
    end    
end

% Kill the empty rows (depths 0 and null)
%max_cc = [max_cc(2:5,:); max_cc(7:10,:)];
%max_ct = [max_ct(2:5,:); max_ct(7:10,:)];
max_cc = [max_cc(2:5) max_cc(7:10)];
max_ct = [max_ct(2:5) max_ct(7:10)];
cc = [cc(2:5,:); cc(7:10,:)];
corrcoefs = [corrcoefs(2:5,:); corrcoefs(7:10,:)];

% Plot peak xcorr for velocity differences vs. response differences
% Title(note);
% for k=1:8
%     subplot(4,2,k);
%     plot(0:5:200,cc(k,:));
% end

% figure(3);
% hold on;
% plot(shifts,corrcoefs);
% plot(shifts,mean(corrcoefs),'ko');
% legend('-2.0','-1.5','1.0','-0.5','0.5','1.0','1.5','2.0','Mean');
% XLabel('Presumed latency (ms)');
% YLabel('Correlation coefficient');
% Title('Eye velocity vs. response, shifted by X ms');

% Write results for this cell to a file
PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\MPEyeResponses.txt';
outfile = PATHOUT;
headerflag = 0;
if (exist(outfile) == 0) % File does not yet exist, so print a header
    headerflag = 1;
end
fid = fopen(outfile, 'a');  % Open text file.
if (headerflag)
    fprintf(fid, 'FILE monkid cellid cc20 cc30 cc40 cc50 cc60 cc70 cc80 cc90 cc100 cc110 cc120 cc130 cc140 cc150 cc160 cc170 cc180 cc190 cc200');
    fprintf(fid, '\r\n');
end
fprintf(fid,'%10s', strtok(FILE,'.'));
fprintf(fid,' %+2.4f', monkid, cellid, mean(corrcoefs));
fprintf(fid,'\r\n');
fclose(fid);

disp('(MPEyeResponses) Done.');
return;