%-----------------------------------------------------------------------------------------------------------------------
%-- MPEyeTracesWrite.m -- Writes eye movement data to a text file for origin plotting. 
%-- Started by JWN, 07/24/07
%-- Last by JWN, 07/24/07
%-----------------------------------------------------------------------------------------------------------------------
function MPEyeTracesWrite(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

ProtocolDefs;
Path_Defs;

disp(sprintf('(MPEyeTracesWrite) Started at %s.',datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);

% In data.eye_data, Channels 1,2,3&4 are eye (x&y), 5&6 are Moog (x&y).
% Only analyze stimulus time 214:614 (2s long).
expansion = 100;
startt = 215 - expansion/2;
endt = 614 + expansion/2;
boxcarwidthpos = 10;
boxcarwidthvel = 40;
eye_xyl = data.eye_data(1:2,startt:endt,:);
eye_xyr = data.eye_data(3:4,startt:endt,:);
Moog_xy = data.eye_data(5:6,startt:endt,:);
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
        disp('(MPSelectiveAnalysis) WARNING: Unknown monkid');
end
% Post-hoc eye calibration
indices = logical(MPtrial_types == 0);
initial_fixation = mean(eye_uv(1,1:expansion/2,indices),3);
offset = mean(initial_fixation);
eye_uv = eye_uv - offset;
% Convert from dva to cm
Moog_uv = Moog_uv * .56;  % ~ Barracuda's viewing distance (31.9) * tan(1/(180/pi))
eye_uv = eye_uv * .56;  % ~ Barracuda's viewing distance (31.9) * tan(1/(180/pi))
% Velocity in deg/s
vMoog_uv = diff(Moog_uv(1,:,:))*200;
veye_uv = diff(eye_uv(1,:,:))*200;
% Compute 16 things
%0 Phase
indices = logical(MPtrial_types == 0 & MPphase == 0); % only care about MP
mpMoog0 = boxcarfilter(squeeze(mean(Moog_uv(1,:,indices),3)),boxcarwidthpos);
vmpMoog0 = [boxcarfilter(squeeze(mean(vMoog_uv(1,:,indices),3)),boxcarwidthvel) 0];
avg0 = boxcarfilter(squeeze(mean(eye_uv(1,:,indices),3)),boxcarwidthpos);
vavg0 = [boxcarfilter(squeeze(mean(veye_uv(1,:,indices),3)),boxcarwidthvel) 0];
indices = logical(MPtrial_types == 0 & MPphase == 0 & MPdepths == -2);  %Near
nears = find(indices);
near0 = boxcarfilter(squeeze(mean(eye_uv(1,:,indices),3)),boxcarwidthpos);
vnear0 = [boxcarfilter(squeeze(mean(veye_uv(1,:,indices),3)),boxcarwidthvel) 0];
indices = logical(MPtrial_types == 0 & MPphase == 0 & MPdepths == 2);  %Far
fars = find(indices);
far0 = boxcarfilter(squeeze(mean(eye_uv(1,:,indices),3)),boxcarwidthpos);
vfar0 = [boxcarfilter(squeeze(mean(veye_uv(1,:,indices),3)),boxcarwidthvel) 0];
%180 Phase
indices = logical(MPtrial_types == 0 & MPphase == 180);
mpMoog180 = boxcarfilter(squeeze(mean(Moog_uv(1,:,indices),3)),boxcarwidthpos);
vmpMoog180 = [boxcarfilter(squeeze(mean(vMoog_uv(1,:,indices),3)),boxcarwidthvel) 0];
avg180 = boxcarfilter(squeeze(mean(eye_uv(1,:,indices),3)),boxcarwidthpos);
vavg180 = [boxcarfilter(squeeze(mean(veye_uv(1,:,indices),3)),boxcarwidthvel) 0];
indices = logical(MPtrial_types == 0 & MPphase == 180 & MPdepths == -2);  %Near
nears180 = find(indices);
near180 = boxcarfilter(squeeze(mean(eye_uv(1,:,indices),3)),boxcarwidthpos);
vnear180 = [boxcarfilter(squeeze(mean(veye_uv(1,:,indices),3)),boxcarwidthvel) 0];
indices = logical(MPtrial_types == 0 & MPphase == 180 & MPdepths == 2);  %Far
fars180 = find(indices);
far180 = boxcarfilter(squeeze(mean(eye_uv(1,:,indices),3)),boxcarwidthpos);
vfar180 = [boxcarfilter(squeeze(mean(veye_uv(1,:,indices),3)),boxcarwidthvel) 0];

%Get individual traces, smoothed
for i = 1:size(nears,2)
    inear0(i,:) = boxcarfilter(squeeze(eye_uv(1,:,nears(i))),boxcarwidthpos);
    ifar0(i,:) = boxcarfilter(squeeze(eye_uv(1,:,fars(i))),boxcarwidthpos);
    vinear0(i,:) = [boxcarfilter(squeeze(veye_uv(1,:,nears(i))),boxcarwidthvel) 0];
    vifar0(i,:) = [boxcarfilter(squeeze(veye_uv(1,:,fars(i))),boxcarwidthvel) 0];
    inear180(i,:) = boxcarfilter(squeeze(eye_uv(1,:,nears180(i))),boxcarwidthpos);
    ifar180(i,:) = boxcarfilter(squeeze(eye_uv(1,:,fars180(i))),boxcarwidthpos);
    vinear180(i,:) = [boxcarfilter(squeeze(veye_uv(1,:,nears180(i))),boxcarwidthvel) 0];
    vifar180(i,:) = [boxcarfilter(squeeze(veye_uv(1,:,fars180(i))),boxcarwidthvel) 0];
end

figure
hold on
plot(mpMoog0,'k')
% plot(avg0,'.')
% plot(near0,'r')
% plot(far0,'b')
plot(mpMoog180,'k')
% plot(avg180,'.')
% plot(near180,'r')
% plot(far180,'b')
plot(inear0','r')
plot(ifar0','g')

figure
hold on
plot(vmpMoog0,'k')
% plot(vavg0,'.')
% plot(vnear0,'r')
% plot(vfar0,'b')
plot(vmpMoog180,'k')
% plot(vavg180,'.')
% plot(vnear180,'r')
% plot(vfar180,'b')
plot(vinear0','r')
plot(vifar0','g')

% Write results for this cell to a file
PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
filenames = {'Supp3'};
for i = 1:1
    outfile = cell2mat(strcat(PATHOUT,strtok(FILE,'.'),'_',filenames(i),'.txt'));
    fid = fopen(outfile, 'w');  % Open text file.
    fprintf(fid, 'Time mpMoog0 near0 far0 avg0 mpMoog180 near180 far180 avg180 vmpMoog0 vnear0 vfar0 vavg0 vmpMoog180 vnear180 vfar180 vavg180 near1 near2 near3 near4 near5 far1 far2 far3 far4 far5 vnear1 vnear2 vnear3 vnear4 vnear5 vfar1 vfar2 vfar3 vfar4 vfar5 180near1 180near2 180near3 180near4 180near5 180far1 180far2 180far3 180far4 180far5 180vnear1 180vnear2 180vnear3 180vnear4 180vnear5 180vfar1 180vfar2 180vfar3 180vfar4 180vfar5');
    fprintf(fid, '\r\n');
    final_data = [(-.295:.005:2.200)' mpMoog0' near0' far0' avg0' mpMoog180' near180' far180' avg180' vmpMoog0' vnear0' vfar0' vavg0' vmpMoog180' vnear180' vfar180' vavg180' inear0' ifar0' vinear0' vifar0' inear180' ifar180' vinear180' vifar180']; % kluged time
    for j = 1:400+expansion
        fprintf(fid,' %+2.4f', final_data(j,:));
        fprintf(fid,'\r\n');
    end
    fclose(fid);
end

disp('(MPEyeTracesWrite) Done.');
return;
