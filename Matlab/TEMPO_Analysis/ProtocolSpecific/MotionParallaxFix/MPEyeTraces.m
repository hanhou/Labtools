%-----------------------------------------------------------------------------------------------------------------------
%-- MPEyeTraces.m -- Plots eye traces. 
%-- Started by JWN, 4/14/04
%-- Last by JWN, 10/03/05
%-----------------------------------------------------------------------------------------------------------------------
function MPEyeTraces(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

%TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
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

disp(sprintf('(MPEyeTraces) Started at %s.',datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);

% In data.eye_data, Channels 1,2,3&4 are eye (x&y), 5&6 are Moog (x&y).
% Only analyze 0.5s to 3.5s (100 to 700 data points at 200Hz)
eye_xyr = data.eye_data(1:2,100:700,:);
eye_xyl = data.eye_data(3:4,100:700,:);
Moog_xy = data.eye_data(5:6,100:700,:);
% Rescale Moog analog signals into degrees
% MoogX = ((amoogx/nSets) * WIN_XRANGE/high_res)/AVAL_RANG
% MoogZ = ((amoogz/nSets) * WIN_YRANGE/high_res)/AVAL_RANGE;
% writef("c:\\TEMPO\\MoogProtocol\\params.log AD_RANGE %8d\n", AVAL_RANGE);	//number of levels in the A/D conversion (e.g., 12 bits -> 4096)
% float high_res = 10.0;	//multiplies the _zoom variables to have high-res animated graphs
% hide float constant WIN_XRANGE = 100.0*high_res;	//full scale will be this many degrees
% hide float constant WIN_YRANGE = 100.0*high_res;
% nSets = 10?

% Moog_xy = Moog_xy / 10 * 1000 / 65536;  We don't need to do this now for some reason?  JWN 100305

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

% Velocity in deg/s
veye_uvr = diff(eye_uvr(1,:,:))*200;
veye_uvl = diff(eye_uvl(1,:,:))*200;
vMoog_uv = diff(Moog_uv(1,:,:))*200;
boxcarwidth = 50;  % Used for later boxcar filtering with boxcarfilter(vector, width)
saccade_threshold = 100;  % deg/s to be considered a saccade-like noise

set(1, 'HandleVisibility', 'off');
figure(2);  set(2, 'HandleVisibility', 'on');  clf(2);
hold on;
set(2,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573],'Name', 'Eye Traces');  % Better for printing?

% Print a header
% note = MPGetNote(FILE,'z:\Users\Jacob\BarracudaNotes.txt');
% note = sprintf('(%d,%d) %d  %s', data.neuron_params(RF_XCTR), data.neuron_params(RF_YCTR), data.neuron_params(RF_DIAMETER), note);

% Plot eye traces
i = 1;  % Only look at MP condition (1=MP, 2=BD, 3=RM, 4=C, 5=EO, 6=HO)
hold on;
indices = logical((MPtrial_types == (i-1)) & MPphase == 0);
mpMoog = mean(vMoog_uv(1,:,indices),3);
mpMoog = boxcarfilter(mpMoog,boxcarwidth);  % Clean Moog movement
plot(.505:.005:3.5,mpMoog,'k');
indices = logical((MPdepths == 2) & (MPtrial_types == (i-1)) & MPphase == 0);  
if(sum(indices)>0)
    farr = mean(veye_uvr(1,:,indices),3);
    farl = mean(veye_uvl(1,:,indices),3);
    farr = boxcarfilter(farr,boxcarwidth);
    farl = boxcarfilter(farl,boxcarwidth);
    far = mean(eye_uv(1,:,indices),3);
    ifarr = squeeze(eye_uvr(1,:,indices));
    ifarl = squeeze(eye_uvl(1,:,indices));
    ifar = squeeze(eye_uv(1,:,indices));
    plot(.505:.005:3.5, farr, 'b');
    plot(.505:.005:3.5, farl, 'b--');
    %plot(.5:.005:3.5, farr+std(eye_uvr(1,:,indices),0,3), 'b--');
    %plot(.5:.005:3.5, farr-std(eye_uvr(1,:,indices),0,3), 'b--');
end
indices = logical((MPdepths == -2) & (MPtrial_types == (i-1)) & MPphase == 0); 
if(sum(indices)>0)
    nearr = mean(veye_uvr(1,:,indices),3);
    nearl = mean(veye_uvl(1,:,indices),3);
    nearr = boxcarfilter(nearr,boxcarwidth);
    nearl = boxcarfilter(nearl,boxcarwidth);
    near = mean(eye_uv(1,:,indices),3);
    inearr = squeeze(eye_uvr(1,:,indices));
    inearl = squeeze(eye_uvl(1,:,indices));
    inear = squeeze(eye_uv(1,:,indices));
    plot(.505:.005:3.5, nearr, 'r');
    plot(.505:.005:3.5, nearl, 'r--');
    %plot(.5:.005:3.5, nearr+std(eye_uvr(1,:,indices),0,3), 'r--');
    %plot(.5:.005:3.5, nearr-std(eye_uvr(1,:,indices),0,3), 'r--');
end
%Title(note);
legend('Moog', 'Far R', 'Far L', 'Near R', 'Near L');
Legend(gca,'boxoff');
YLabel('Eye position (deg)');
xlim([0.7 3.3]);
ylim([-12 12]);
% ylim([-4 4]);
% diffs(1,:) = farr-nearr;
% verg(1,:) = farl-farr;
% verg(2,:) = nearl-nearr;

% Lets do a 2-way anova (time and depth)
% reye_anova_vals = zeros(size(MPdepths,2)/8, size(nearr,2));  % 50 rows, 601 columns
% reye_anova_vals180 = zeros(size(MPdepths,2)/8, size(nearr,2));  % 50 rows, 601 columns
% leye_anova_vals = zeros(size(MPdepths,2)/8, size(nearr,2));  % 50 rows, 601 columns
% leye_anova_vals180 = zeros(size(MPdepths,2)/8, size(nearr,2));  % 50 rows, 601 columns
% reps = size(reye_anova_vals,1)/num_depths;  % 5 reps, usually
% for j=1:num_depths
%     indices = logical((MPdepths == uMPdepths(j)) & (MPtrial_types == (i-1)) & MPphase == 0);
%     reye_anova_vals(5*(j-1)+1:5*(j-1)+5,:) = squeeze(eye_uvr(1,:,indices))';
%     leye_anova_vals(5*(j-1)+1:5*(j-1)+5,:) = squeeze(eye_uvl(1,:,indices))';
%     indices = logical((MPdepths == uMPdepths(j)) & (MPtrial_types == (i-1)) & MPphase == 180);
%     reye_anova_vals180(5*(j-1)+1:5*(j-1)+5,:) = squeeze(eye_uvr(1,:,indices))';
%     leye_anova_vals180(5*(j-1)+1:5*(j-1)+5,:) = squeeze(eye_uvl(1,:,indices))';
% end
% pr = anova2(reye_anova_vals,reps,'off');
% pl = anova2(leye_anova_vals,reps,'off');
% pr180 = anova2(reye_anova_vals180,reps,'off');
% pl180 = anova2(leye_anova_vals180,reps,'off');

% subplot(5,1,2);
% hold on;
% indices = logical((MPtrial_types == (i-1)) & MPphase == 180);
% plot(.5:.005:3.5,mean(Moog_uv(1,:,indices),3),'k');
% indices = logical((MPdepths == 1.5 | MPdepths == 2) & (MPtrial_types == (i-1)) & MPphase == 180);  
% if(sum(indices)>0)
%     farr = mean(eye_uvr(1,:,indices),3);
%     farl = mean(eye_uvl(1,:,indices),3);
%     plot(.5:.005:3.5, farr, 'b');
%     plot(.5:.005:3.5, farl, 'b--');    
% end
% indices = logical((MPdepths == -1.5 | MPdepths == -2) & (MPtrial_types == (i-1)) & MPphase == 180); 
% if(sum(indices)>0)
%     nearr = mean(eye_uvr(1,:,indices),3);
%     nearl = mean(eye_uvl(1,:,indices),3);
%     plot(.5:.005:3.5, nearr, 'r');
%     plot(.5:.005:3.5, nearl, 'r--');    
% end
% YLabel('Eye position (deg)');
% diffs(2,:) = farr-nearr;
% verg(3,:) = farl-farr;
% verg(4,:) = nearl-nearr;
% 
% subplot(5,1,3);
% hold on;
% indices = logical((MPtrial_types == (i-1)) & MPphase == 0);
% vcorrect = boxcarfilter(mean(vMoog_uv(1,:,indices),3),boxcarwidth);  % Clean Moog movement
% plot(.5:.005:3.495,vcorrect,'k');  % Plotted Moog movement
% indices = logical((MPdepths == 1.5 | MPdepths == 2) & (MPtrial_types == (i-1)) & MPphase == 0);  
% if(sum(indices)>0)
%     % Removing saccade-like noise (and maybe an occasional saccade) from right eye
%     raw = veye_uvr(1,:,indices);
%     sacs = find(raw>saccade_threshold|raw<-saccade_threshold);  size(sacs)
%     sacscorrect = mod(sacs,600);
%     if(size(sacs,1)>0)
%         raw(sacs) = vcorrect(sacscorrect);
%     end
%     farr = boxcarfilter(mean(raw,3),boxcarwidth);
%     % Removing saccade-like noise (and maybe an occasional saccade) from left eye
%     raw = veye_uvl(1,:,indices);
%     sacs = find(raw>saccade_threshold|raw<-saccade_threshold);  size(sacs)
%     sacscorrect = mod(sacs,600);
%     if(size(sacs,1)>0)
%         raw(sacs) = vcorrect(sacscorrect);
%     end
%     farl = boxcarfilter(mean(raw,3),boxcarwidth);
%     % Plot
%     plot(.5:.005:3.495, farr, 'b');
%     plot(.5:.005:3.495, farl, 'b--');
% end
% indices = logical((MPdepths == -1.5 | MPdepths == -2) & (MPtrial_types == (i-1)) & MPphase == 0);  
% if(sum(indices)>0)
%     % Removing saccade-like noise (and maybe an occasional saccade) from right eye
%     raw = veye_uvr(1,:,indices);
%     sacs = find(raw>saccade_threshold|raw<-saccade_threshold);  size(sacs)
%     sacscorrect = mod(sacs,600);
%     if(size(sacs,1)>0)
%         raw(sacs) = vcorrect(sacscorrect);
%     end
%     nearr = boxcarfilter(mean(raw,3),boxcarwidth);
%     % Removing saccade-like noise (and maybe an occasional saccade) from left eye
%     raw = veye_uvl(1,:,indices);
%     sacs = find(raw>saccade_threshold|raw<-saccade_threshold);  size(sacs)
%     sacscorrect = mod(sacs,600);
%     if(size(sacs,1)>0)
%         raw(sacs) = vcorrect(sacscorrect);
%     end
%     nearl = boxcarfilter(mean(raw,3),boxcarwidth);
%     % Plot
%     plot(.5:.005:3.495, nearr, 'r');
%     plot(.5:.005:3.495, nearl, 'r--');
% end
% YLabel('Eye velocity (deg/s)');
% diffs(3,:) = [farr 0]-[nearr 0];
% 
% subplot(5,1,4);
% hold on;
% indices = logical((MPtrial_types == (i-1)) & MPphase == 180);
% vcorrect = boxcarfilter(mean(vMoog_uv(1,:,indices),3),boxcarwidth);  % Clean Moog movement
% plot(.5:.005:3.495,vcorrect,'k');  % Plotted Moog movement
% indices = logical((MPdepths == 1.5 | MPdepths == 2) & (MPtrial_types == (i-1)) & MPphase == 180); 
% if(sum(indices)>0)
%     % Removing saccade-like noise (and maybe an occasional saccade) from right eye
%     raw = veye_uvr(1,:,indices);
%     sacs = find(raw>saccade_threshold|raw<-saccade_threshold);  size(sacs)
%     sacscorrect = mod(sacs,600);
%     if(size(sacs,1)>0)
%         raw(sacs) = vcorrect(sacscorrect);
%     end
%     farr = boxcarfilter(mean(raw,3),boxcarwidth);
%     % Removing saccade-like noise (and maybe an occasional saccade) from left eye
%     raw = veye_uvl(1,:,indices);
%     sacs = find(raw>saccade_threshold|raw<-saccade_threshold);  size(sacs)
%     sacscorrect = mod(sacs,600);
%     if(size(sacs,1)>0)
%         raw(sacs) = vcorrect(sacscorrect);
%     end
%     farl = boxcarfilter(mean(raw,3),boxcarwidth);
%     % Plot
%     plot(.5:.005:3.495, farr, 'b');
%     plot(.5:.005:3.495, farl, 'b--');
% end
% indices = logical((MPdepths == -1.5 | MPdepths == -2) & (MPtrial_types == (i-1)) & MPphase == 180);  
% if(sum(indices)>0)
%     % Removing saccade-like noise (and maybe an occasional saccade) from right eye
%     raw = veye_uvr(1,:,indices);
%     sacs = find(raw>saccade_threshold|raw<-saccade_threshold);  size(sacs)
%     sacscorrect = mod(sacs,600);
%     if(size(sacs,1)>0)
%         raw(sacs) = vcorrect(sacscorrect);
%     end
%     nearr = boxcarfilter(mean(raw,3),boxcarwidth);
%     % Removing saccade-like noise (and maybe an occasional saccade) from left eye
%     raw = veye_uvl(1,:,indices);
%     sacs = find(raw>saccade_threshold|raw<-saccade_threshold);  size(sacs)
%     sacscorrect = mod(sacs,600);
%     if(size(sacs,1)>0)
%         raw(sacs) = vcorrect(sacscorrect);
%     end
%     nearl = boxcarfilter(mean(raw,3),boxcarwidth);
%     % Plot
%     plot(.5:.005:3.495, nearr, 'r');
%     plot(.5:.005:3.495, nearl, 'r--');
% end
% YLabel('Eye velocity (deg/s)');
% diffs(4,:) = [farr 0]-[nearr 0];
% 
% subplot(5,1,5);
% hold on;
% 
% % For velocity differences
% % for j = 1:4
% %     plot(.5:.005:3.5,diffs(j,:),colors{j});
% % end
% % legend('Position 0', 'Position 180', 'Velocity 0', 'Velocity 180');
% % Legend(gca,'boxoff');
% % YLabel('Far-Near (deg or deg/s)');
% % h = axis;
% % h(3:4) = [-2 2];
% % axis(h);
% 
% % For vergence differences
% for j = 1:4
%     plot(.5:.005:3.5,verg(j,:),colors{j});
% end
% legend('Far 0', 'Near 180', 'Far 0', 'Near 180');
% Legend(gca,'boxoff');
% h = axis;
% h(3:4) = [-0.4 0.4];
% axis(h);
% YLabel('L Eye - R Eye (deg)');

XLabel('Time (s)');

%Loop back through and slap on begin/end times
begin_time = 1;  % Onset
end_time = 3;  % Offset
%for j = 1:5;
    %subplot(5,1,j);
    hold on;
    h = axis;
    plot([begin_time begin_time], [h(3) h(4)], 'k--');
    plot([end_time end_time], [h(3) h(4)],'k--');
%end

titles = { 'Motion Parallax', 'Binocular Disparity', 'Retinal Motion', 'Congruent', 'Eye Movement Only', 'Head Movement Only' };
Title(sprintf('%s\n %s ',FILE,titles{i}));

% Write results for this cell to a file
% PATHOUT = 'Z:\Data\MOOG\Barracuda\Analysis\EyeTraces\';
% filenames = {'EyePos'};
% for i = 1:1
%     outfile = cell2mat(strcat(PATHOUT,strtok(FILE,'.'),'_',filenames(i),'.txt'));
%     fid = fopen(outfile, 'w');  % Open text file.
%     fprintf(fid, 'Time AvgMoog AvgNear AvgFar Near1 Near2 Near3 Near4 Near5 Far1 Far2 Far3 Far4 Far5');
%     fprintf(fid, '\r\n');
%     final_data = [(-565:5:2435)' mpMoog' near' far' inear ifar]; % kluged time
%     for j = 1:601
%         fprintf(fid,' %+2.4f', final_data(j,:));
%         fprintf(fid,'\r\n');
%     end
%     fclose(fid);
% end

disp('(MPEyeTraces) Done.');
return;