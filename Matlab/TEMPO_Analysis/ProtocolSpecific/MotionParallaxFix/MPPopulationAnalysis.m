%-----------------------------------------------------------------------------------------------------------------------
%-- MPPopulationAnalyses.m -- Do everything, and save the information to a file.
%-- Started by JWN, 6/17/05
%-- Last by JWN, 6/17/05
%-----------------------------------------------------------------------------------------------------------------------
function MPPopulationAnalyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

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

area = 'MT';
disp(sprintf('(MPPopulationAnalysis) Started at %s.  Analyzing %s %s',datestr(now,14),area,FILE));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);

% Extract the bin indices corresponding to the codes
% I am assuming all trials have same begin and end times!
begin_time = find(data.event_data(1,:,1)==StartCode);
end_time = find(data.event_data(1,:,1)==StopCode);
total_spike_bins = end_time - begin_time;

uMPtrial_types = unique(MPtrial_types);
uMPdepths = unique(MPdepths);
uMPphase = unique(MPphase);
num_trial_types = size(uMPtrial_types,2);
num_depths = size(uMPdepths,2);
num_phase = size(uMPphase,2);

% Use BD condition to compute latency for this cell
indices = logical(MPtrial_types == 1 & MPphase == 0 & MPdepths ~= -9999);
try
    latencies = Latency(data.spike_data(1,begin_time-200:end_time,indices),200);
catch
    if(strcmp(area,'MT'))
        latencies = [110 110 110];  % MT guess
    else
        latencies = [130 130 130];  % MST guess
    end
end
pseudo_latency = latencies(2);
delay_modifier = 50;  % A guess (in ms) about delay between stim on command and actual stim on
estimated_true_latency = pseudo_latency - delay_modifier;  % Take into account that delay

%%%  Calculate PDI and PDImod  %%%
begint = begin_time + pseudo_latency;  % Modified begin and end times
endt = end_time + pseudo_latency;
raw_spikes = data.spike_data(1,begint:endt,:);
spike_rates = 1000*squeeze(mean(raw_spikes))';  % Mean spike rates
for j = 1:num_depths-1
    mean_data(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == 0)';    
end
mis = MPGetMI(pseudo_latency, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin);
meanmean_data = mean(mean_data);
meanmis = mean(mis,3);
stdmean_data = std(mean_data);
stdmis = std(mis,0,3);
for k = 1:num_depths/2-1
    nearm = meanmean_data(k);
    farm = meanmean_data(num_depths-k);
    nearstd = stdmean_data(k);
    farstd = stdmean_data(num_depths-k);
    PDI(k) = (farm-nearm)/(abs(farm-nearm)+(nearstd+farstd)/2);
    nearm = meanmis(1,k+1);  % Remember to shift by one to lose spont data
    farm = meanmis(1,k+6);  % and by six to get to the far data
    nearstd = stdmis(1,k+1);
    farstd = stdmis(1,k+6);
    PDImod(k) = (farm-nearm)/(abs(farm-nearm)+(nearstd+farstd)/2);
end
mPDI = mean(PDI);
mPDImod = mean(PDImod); 

%%%  Compute Eye-Response correlations  %%%

% Defines for use in eye response binning (5ms bins to match velocity data)
num_reduced_bins = 399;
bin_width = total_spike_bins/(num_reduced_bins+1);  % ~2000ms/(399+1) = ~5ms;

% Get eye and Moog data
eye_xy = data.eye_data(3:4,round((begin_time+delay_modifier)/5):round((begin_time+delay_modifier)/5)+num_reduced_bins+1,:);
Moog_xy = data.eye_data(5:6,round((begin_time+delay_modifier)/5):round((begin_time+delay_modifier)/5)+num_reduced_bins+1,:);
Moog_xy = Moog_xy / 10 * 1000 / 65536;

% Realign axes to match preferred direction
pref = data.neuron_params(PREFERRED_DIRECTION);
opp = tan(pref/(180/pi));
u = [1 opp] / sqrt(1+opp^2);
v = [-u(2) u(1)];
for i=1:size(eye_xy,3)
    eye_uv(1,:,i) = u*eye_xy(:,:,i);
    eye_uv(2,:,i) = v*eye_xy(:,:,i);
    Moog_uv(1,:,i) = u*Moog_xy(:,:,i);
    Moog_uv(2,:,i) = v*Moog_xy(:,:,i);
end

% Velocity in deg/s
veye_uv = diff(eye_uv(1,:,:))*200;
vMoog_uv = diff(Moog_uv(1,:,:))*200;
boxcarwidth = 50;  % Used for later boxcar filtering with boxcarfilter(vector, width)
saccade_threshold = 100;  % deg/s to be considered a saccade-like noise

i = 1;  % Only look at MP condition
for j = 1:num_depths  % Look at all 10 depths...
    if(uMPdepths(j)==0 | uMPdepths(j)==-9999) continue; end  % but throw away two of them (0 and null)
    indices_0 = logical((MPtrial_types == (i-1)) & MPphase == 0 & MPdepths == uMPdepths(j));
    trials_0 = find(indices_0==1);  % Get zero phase trials for this depth
    indices_180 = logical((MPtrial_types == (i-1)) & MPphase == 180 & MPdepths == -uMPdepths(j));
    trials_180 = find(indices_180==1);  % Get 180 phase trials for the -matched depth
    m_veye_uv_0 = mean(veye_uv(1,:,indices_0),3); % Get zero phase mean for this depth
    m_veye_uv_180 = mean(veye_uv(1,:,indices_180),3);  % Get 180 phase mean for the -matched depth
    delta_veye(j,:) = boxcarfilter(m_veye_uv_0,boxcarwidth) + boxcarfilter(m_veye_uv_180,boxcarwidth);
    
    raw_spikes = data.spike_data(1,begint:endt,indices_0);
    hist_data = sum(raw_spikes,3);
    [bins, counts0] = SpikeBinner(hist_data, 1, bin_width, 0);
    counts0 = counts0*(1000/bin_width)/sum(indices_0);  % Convert to instantaneous firing rates
    raw_spikes = data.spike_data(1,begint:endt,indices_180);
    hist_data = sum(raw_spikes,3);
    [bins, counts180] = SpikeBinner(hist_data, 1, bin_width, 0);
    counts180 = counts180*(1000/bin_width)/sum(indices_180);  % Convert to instantaneous firing rates
    % Compute the differences in response, and smooth with a boxcar
    delta_resp(j,:) = boxcarfilter((counts0-counts180)',boxcarwidth);
    
    % Compute single Spearman rank correlation coefficients for a latency shifted time
    lag = round(0/5);  % Replace 0 to put in lags that are more or less than estimated_true_latency
    [r,p] = corr(delta_veye(j,1+lag:end)',delta_resp(j,1:end-lag)','type','Spearman');
    ccrs(j) = r;
    ccps(j) = p;
end
% Compute a giant correlation for all depths & phases concatenated together
alldelta_veye = reshape([delta_veye(2:5,:) delta_veye(7:10,:)],3200,1);
alldelta_resp = reshape([delta_resp(2:5,:) delta_resp(7:10,:)],3200,1);
[r,p] = corr(alldelta_veye,alldelta_resp,'type','Spearman');
allr = r;
allp = p;

% Now write it all to a file.
fid = fopen('Z:\Data\MOOG\Barracuda\Analysis\PopulationAnalysis.txt', 'a');  % Open text file.
fprintf(fid,'%10s %4d %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f %+2.4f', strtok(FILE,'.'), estimated_true_latency, PDI, PDImod, mPDI, mPDImod, [ccrs(2:5) ccrs(7:10)], [ccps(2:5) ccps(7:10)], allr, allp);
fprintf(fid,'\r\n');
fclose(fid);

% Plot out Eye-response figures
set(1, 'HandleVisibility', 'off');
figure(2);  set(2, 'HandleVisibility', 'on');  clf(2);
hold on;
set(2,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573],'Name', 'Eye-Response plots');  % Better for printing?
for i = 1:4
    subplot(4,3,(i-1)*3+1);
    hold on;
    plot(delta_veye(i+1,:),'b');
    plot(delta_veye((num_depths-i)+1,:),'r');
    if(i==1)
        Title('delta EYE');
    end
    if(i==4)
        Xlabel('Time (5ms bins)');
    end
    subplot(4,3,(i-1)*3+2);
    hold on;
    plot(delta_resp(i+1,:),'b');
    plot(delta_resp((num_depths-i)+1,:),'r');
    if(i==1)
        note = sprintf('%s %s PDI=%2.4f PDImod=%2.4f (R=%2.4f P=%2.4f)',area,FILE,mPDI,mPDImod,allr,allp);
        Title(sprintf('%s\n%s',note,'delta RESP'));
    end
    if(i==4)
        Xlabel('Time (5ms bins)');
    end
    subplot(4,3,(i-1)*3+3);
    hold on;
    plot(delta_resp(i+1,:),delta_veye(i+1,:),'b.');
    plot(delta_resp((num_depths-i)+1,:),delta_veye((num_depths-i)+1,:),'r.');
    if(i==1)
        Title('dRESP(x) vs dEYE(y)');
    end
    Legend(sprintf('R=%+2.4f P=%+2.4f',ccrs(i+1),ccps(i+1)),sprintf('R=%+2.4f P=%+2.4f',ccrs((num_depths-i)+1),ccps((num_depths-i)+1)));
    Legend(gca,'boxoff');
end

% Print this figure
% print(2);

disp('(MPPopulationAnalysis) Done.');
return;