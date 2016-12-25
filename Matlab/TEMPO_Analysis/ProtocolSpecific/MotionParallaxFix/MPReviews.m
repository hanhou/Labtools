%-----------------------------------------------------------------------------------------------------------------------
%-- MPEyeStdDevs.m - Gets the std deviations of eye position (along prefdir axis aka axis of motion) for all
%-- conditions to see if there is more wiggle in the HO case than RM.  Named MPREviews because no access to ProtocolDefs
%-- Started by JWN, 10/25/08
%-- Last by JWN, 10/25/08
%-----------------------------------------------------------------------------------------------------------------------
function MPReviews(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

ver = '1.0';
TEMPO_Defs;
Path_Defs;
symbols = {'bo' 'rs' 'gd' 'kv' 'm<' 'c>' 'bv' 'rv'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'm<-' 'c>-' 'yo-' 'bs-'};
NULL_VALUE = -9999;
area = 'MT'
disp(sprintf('(MPEyeStdDevs v%s) Started at %s.',ver,datestr(now,14)));

[monkid, cellid, runstr]=strread(FILE,'m%dc%dr%s.htb');
% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
uMPtrial_types = unique(MPtrial_types);  % Conditions present
%Place breakouts here.  This is what SelectiveAnalysis could be all about!
%if(isempty(find(uMPtrial_types==###))) return;  end;
%if(isempty(find(uMPtrial_types==0))) disp('(MPSelectiveAnalysis) Breakout: No MP');  return;  end;  % BREAKOUT ENABLED!

num_trial_types = length(uMPtrial_types);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
uMPphase = unique(MPphase);
num_phase = size(uMPphase,2);
if(num_phase ~= 2)
    disp('(MPEyeStdDevs) Fatal Error: Two phases required to calculate modulation indices.');
    return;
end
trials = size(MPphase,2);

pref = data.neuron_params(PREFERRED_DIRECTION);

% Take a break from firing to look at eye movements
% In data.eye_data, Channels 1,2,3&4 are eye (x&y), 5&6 are Moog (x&y).
% Only analyze stimulus time 214:614 (2s long).
eye_xyl = data.eye_data(1:2,215:614,:);
eye_xyr = data.eye_data(3:4,215:614,:);
Moog_xy = data.eye_data(5:6,215:614,:);
% Realign axes to match preferred direction
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
switch monkid
    case 9, % Barracuda
        if cellid==155, eye_uv = eye_uvl; end
        if cellid==156, eye_uv = eye_uvl; end
        gain_constant = .18;  % Based on viewing distance and interocular distance
    case 15, % Ovid
        if cellid<35, eye_uv = eye_uvl; end
        if cellid>46, eye_uv = eye_uvr; end
        gain_constant = .15;  % Based on viewing distance and interocular distance
    otherwise
        disp('(MPEyeStdDevs) WARNING: Unknown monkid');
end
%Get stddev now for each condition
eyeposstd = zeros(1,6)+888;
for i = 1:6
   if(isempty(find(uMPtrial_types==i-1))) continue;  end;  % Break out if none from that condition
   indices = (MPtrial_types == i-1);
   eyeposstd(i) = mean(std(eye_uv(1,:,indices)-Moog_uv(1,:,indices)));
end

% Write results for this cell to 1 file
PATHOUT = 'Z:\Users\Jacob\';
filenames = {'MPEyeStdDevsJ'};
for i = 1:1
    outfile = cell2mat(strcat(PATHOUT,area,'_',filenames(i),'.txt'));
    headerflag = 0;
    if (exist(outfile) == 0) % File does not yet exist, so print a header
        headerflag = 1;
    end
    fid = fopen(outfile, 'a');  % Open text file.
    if (headerflag)
        fprintf(fid, 'FILE ');
        fprintf(fid, 'monkid cellid ');
        fprintf(fid, 'eyeposstdMP eyeposstdBD eyeposstdRM eyeposstdC eyeposstdEO eyeposstdHO ');
        fprintf(fid, '\r\n');
    end
    fprintf(fid,'%10s', strtok(FILE,'.'));
    fprintf(fid,' %+2.5f', monkid, cellid, eyeposstd(1:6));
    fprintf(fid,'\r\n');
    fclose(fid);
end

disp('(MPEyeStdDevs) Done.');
return;