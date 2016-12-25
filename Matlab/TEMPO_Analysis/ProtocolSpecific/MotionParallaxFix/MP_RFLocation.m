%-----------------------------------------------------------------------------------------------------------------------
%-- MP_RFLocation.m -- Used in analyses based on RF locations.
%-- Started by JWN, 10/30/05
%-- Last by JWN, 10/30/05
%-----------------------------------------------------------------------------------------------------------------------
function MP_RFLocation(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

TEMPO_Defs;
ProtocolDefs;
Path_Defs;

symbols = {'bo' 'rs' 'gd' 'kv' 'bd' 'rd' 'bv' 'rv'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'mo-' 'co-' 'yo-' 'bs-'};
NULL_VALUE = -9999;

disp(sprintf('(MP_RFLocation) Started at %s.',datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
% Get the mean firing rates for all the trials

area = 'MT';  % Kluge! 80 for MT, 100 for MST, +50 for transfer function delay 

if(strcmp(area,'MT'))  % Don't change this one!
    latency = 130;  % MT guess
else
    latency = 150;  % MST guess
end 
begin_time = find(data.event_data(1,:,1)==StartCode) + latency;
end_time = find(data.event_data(1,:,1)==StopCode) + latency;
if(max(max(max(data.spike_data))) > 1)
    disp(sprintf('(MPDepthTuningCurve) WARNING: %d corrupt values in data.spike_data.',sum(sum(sum(data.spike_data>1)))));
    data.spike_data = cast(data.spike_data>0,'double');
end
raw_spikes = data.spike_data(1,begin_time:end_time,:);
spike_rates = 1000*squeeze(mean(raw_spikes))';  % The hard way

% Calculate number of trials
trials = size(MPphase,2);

%%%  Calculate PDI and PDImod  %%%
% Now also calculate iPDIs (PDIs for the four individual conditions)
mis = MPGetMI(latency, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin);
if(size(unique(MPtrial_types),2)==4) % Only for full sets
    for i = 1:5 % 5 is the true PDI (MP-RM)
        if i>4 % true PDI
            for j = 1:num_depths-1
                mean_dataMP(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == 0)';
                mean_dataRM(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == 2)'; 
            end
            mean_data = mean_dataMP - mean_dataRM;
            misMP = squeeze(mis(1,:,:));
            misRM = squeeze(mis(3,:,:));
            mis_data = misMP-misRM;
        else  % iPDIs
            for j = 1:num_depths-1
                mean_data(:,j) = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == i-1)';
            end
            mis_data = squeeze(mis(i,:,:));
        end
        meanmean_data = mean(mean_data);  % mean of the trials from the means measure
        meanmis = mean(mis_data,2);  % mean of the trials from the mis measure
        stdmean_data = std(mean_data);
        stdmis = std(mis_data,0,2);
        for k = 1:num_depths/2-1
            nearm = meanmean_data(k);
            farm = meanmean_data(num_depths-k);
            nearstd = stdmean_data(k);
            farstd = stdmean_data(num_depths-k);
            PDI(k) = (farm-nearm)/(abs(farm-nearm)+(nearstd+farstd)/2);
            nearm = meanmis(k+1);  % Remember to shift by one to lose spont data
            farm = meanmis(num_depths-(k-1));  % and by six to get to the far data
            nearstd = stdmis(k+1);
            farstd = stdmis(num_depths-(k-1));
            PDImod(k) = (farm-nearm)/(abs(farm-nearm)+(nearstd+farstd)/2);
        end
        mPDI(i) = mean(PDI);
        mPDImod(i) = mean(PDImod);
        % Send data off to MPBootstrap for significance testing.
        pmPDI(i) = MPBootstrap(mean_data, mPDI(i));
        pmPDImod(i) = MPBootstrap(mis_data(2:end,:)', mPDImod(i));
    end
else  % if not full set, just zero things out so you don't bug the figure
    PDI = zeros(1,4);
    PDImod = zeros(1,4);
    mPDI = zeros(1,5);
    mPDImod = zeros(1,5);
end    

% Get RF information
RFx = data.neuron_params(RF_XCTR);
RFy = data.neuron_params(RF_YCTR);
RFd = data.neuron_params(RF_DIAMETER);
pref = data.neuron_params(PREFERRED_DIRECTION);
% Calculate eccentricity, polar angle, and realigned polar angle (based on pref)
RFecc = sqrt(RFx^2 + RFy^2);
RFang = atan2(RFy,RFx)*180/pi;  % between -180 and 180
if(RFang<0)
    RFang = 360+RFang; % make it between 0 and 360 like pref is
end
rRFang = abs(pref-RFang);  % larger angle minus smaller angle
if(rRFang>=180)
    rRFang = 360-rRFang; % get the smaller of the two angle differences (the one less than 180)
end
if(rRFang>=90)
    rRFang = 180-rRFang; % convert to angle between 0 and 90
end
if(rRFang-45 < -15)
    sigrRFang = 1;
else if(rRFang-45 >15)
    sigrRFang = 3;
else
    sigrRFang = 2;
    end
end     

% Write results for this cell to 1 file
PATHOUT = 'Z:\Data\MOOG\Barracuda\Analysis\';
filenames = {'RFLocations'};
for i = 1:1
    outfile = cell2mat(strcat(PATHOUT,area,'_',filenames(i),'.txt'));
    headerflag = 0;
    if (exist(outfile) == 0) % File does not yet exist, so print a header
        headerflag = 1;
    end
    fid = fopen(outfile, 'a');  % Open text file.
    if (headerflag)
        fprintf(fid, 'FILE RFx RFy RFd RFecc Pref RFang rRFang sigrRFang MPiPDI BDiPDI RMiPDI CiPDI PDI absMPiPDI absBDiPDI absRMiPDI absCiPDI absPDI sigpMPiPDI sigpBDiPDI sigpRMiPDI sigpCiPDI sigpPDI MPiPDImod BDiPDImod RMiPDImod CiPDImod PDImod absMPiPDImod absBDiPDImod absRMiPDImod absCiPDImod absPDImod sigpMPiPDImod sigpBDiPDImod sigpRMiPDImod sigpCiPDImod sigpPDImod');
        fprintf(fid, '\r\n');
    end
    fprintf(fid,'%10s', strtok(FILE,'.'));
    fprintf(fid,' %+2.3f', RFx, RFy, RFd, RFecc, pref, RFang, rRFang, sigrRFang);
    fprintf(fid,' %+2.5f', mPDI, abs(mPDI), (pmPDI>0.025)+1, mPDImod, abs(mPDI), (pmPDI>0.025)+1);    
    fprintf(fid,'\r\n');
    fclose(fid);
end

disp('(MP_RFLocation) Done.');
return;