%-----------------------------------------------------------------------------------------------------------------------
%-- MPTrainingGains.m -- Comes from SelectiveAnalaysis.m
%-- Started by JWN, 10/22/06
%-----------------------------------------------------------------------------------------------------------------------
function MPTrainingGains(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

ver = '1.0';
TEMPO_Defs;
ProtocolDefs;
Path_Defs;
symbols = {'bo' 'rs' 'gd' 'kv' 'm<' 'c>' 'bv' 'rv'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'm<-' 'c>-' 'yo-' 'bs-'};
NULL_VALUE = -9999;

disp(sprintf('(MPTrainingGains v%s) Started at %s.',ver,datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
uMPtrial_types = unique(MPtrial_types);  % Conditions present
%Place breakouts here.  This is what SelectiveAnalysis could be all about!
%if(isempty(find(uMPtrial_types==###))) return;  end;

%if(isempty(find(uMPtrial_types==0))) disp('(MPTrainingGains) Breakout: No MP');  return;  end;  % BREAKOUT ENABLED!

num_trial_types = length(uMPtrial_types);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
uMPphase = unique(MPphase);
num_phase = size(uMPphase,2);
if(num_phase ~= 2)
    disp('(MPTrainingGains) Fatal Error: Two phases required to calculate modulation indices.');
    return;
end
trials = size(MPphase,2);
reps = floor(trials/(num_trial_types*num_depths*num_phase));

% In data.eye_data, Channels 1,2,3&4 are eye (x&y), 5&6 are Moog (x&y).
% Only analyze stimulus time 200:642.
eye_xyr = data.eye_data(1:2,200:642,:);
eye_xyl = data.eye_data(3:4,200:642,:);
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
eye_uv = eye_uvr; % Only use left eye because Ovid only has left eye, r is left eye
% Velocity in deg/s
veye_uvr = diff(eye_uvr(1,:,:))*200;
veye_uvl = diff(eye_uvl(1,:,:))*200;
vMoog_uv = diff(Moog_uv(1,:,:))*200;
veye_uv = diff(eye_uv(1,:,:))*200;
%indices = logical((MPtrial_types == 4) | (MPtrial_types == 0)); % only care about MP/EO, lumping phases together
indices = logical((MPtrial_types == 4)); % only care about MP or EO, lumping phases together
% Do fft on Moog first for baseline, throwing out v component
fft_vMoog = squeeze(abs(fft(vMoog_uv(1,:,indices)))); %All ffts
Moog_amplitude = mean(fft_vMoog(2,:));
median_Moog_amplitude = median(fft_vMoog(2,:));
% Do fft on average eye
fft_veye = squeeze(abs(fft(veye_uv(1,:,indices)))); %All ffts
eye_amplitude = mean(fft_veye(2,:));
median_eye_amplitude = median(fft_veye(2,:));
pursuit_gain = eye_amplitude/Moog_amplitude
median_pursuit_gain = median_eye_amplitude/median_Moog_amplitude
% Do near and far too
indices = logical((MPtrial_types == 0) & (MPdepths == -2) & MPphase == 0);
fft_veye = squeeze(abs(fft(veye_uv(1,:,indices)))); %All ffts
eye_amplitude = mean(fft_veye(2,:));
near_gain = eye_amplitude/Moog_amplitude
indices = logical((MPtrial_types == 0) & (MPdepths == 2) & MPphase == 0);
fft_veye = squeeze(abs(fft(veye_uv(1,:,indices)))); %All ffts
eye_amplitude = mean(fft_veye(2,:));
far_gain = eye_amplitude/Moog_amplitude
        
% Get properties of movement (magnitude, direction, eyewindow size) if possible
move_magnitude = mean(data.moog_params(MOVE_MAGNITUDE,:,MOOG))
pref = data.neuron_params(PREFERRED_DIRECTION);
eyewinx = data.targ_params(TARG_WIN_FWIDTH);
eyewiny = data.targ_params(TARG_WIN_FHEIGHT);
if(pref < 30 | pref > 330 | pref>150&pref<210)
    pref_cat = 1;
    relevanteyewin = eyewinx;
elseif(pref>60&pref<120 | pref>240&pref<300)
    pref_cat = 3;
    relevanteyewin = eyewiny;
else
    pref_cat = 2;
    relevanteyewin = sqrt(eyewinx^2+eyewiny^2);
end

% Write results for this cell to 1 file
% PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
% filenames = {'TrainingGains'};
% for i = 1:1
%     outfile = cell2mat(strcat(PATHOUT,filenames(i),'.txt'));
%     headerflag = 0;
%     if (exist(outfile) == 0) % File does not yet exist, so print a header
%         headerflag = 1;
%     end
%     fid = fopen(outfile, 'a');  % Open text file.
%     if (headerflag)
%         fprintf(fid, 'FILE ');
%         fprintf(fid, 'MPgain MPmediangain MPneargain MPfargain ');
%         fprintf(fid, 'magnitude pref prefcat eyewinX eyewinY reyewin');
%         fprintf(fid, '\r\n');
%     end
%     fprintf(fid,'%10s', strtok(FILE,'.'));
%     fprintf(fid,' %+2.5f', pursuit_gain, median_pursuit_gain, near_gain, far_gain, move_magnitude, pref, pref_cat, eyewinx, eyewiny, relevanteyewin);
%     fprintf(fid,'\r\n');
%     fclose(fid);
% end

disp('(MPTrainingGains) Done.');
return;