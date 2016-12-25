% for Tilt/Trans task   modified DirectionTuning2D_pursuit
%--	KT, 2/14/06, MODIFIED 3/23/06
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuning2D_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_elevation = data.moog_params(ROT_AZIMUTH,:,MOOG);%That's tric! only change ELE to AZI
temp_fp_rotate = data.moog_params(TT_MODE,:,MOOG);%Same Tric! only change FP-ROTATE to TT_MODE
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             
%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_elevation);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
null_trials = logical( (temp_elevation == data.one_time_params(NULL_VALUE)) );
elevation = temp_elevation(~null_trials & select_trials);
fp_rotate = temp_fp_rotate(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials)

unique_elevation = munique(elevation');
unique_fp_rotate = munique(fp_rotate');

h_title{1}='Tilt+Trans';
h_title{2}='Tilt-Trans';
h_title{3}='Tilt only';
h_title{4}='Trans only';

% calculate spontaneous firing rate 
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));
repetition = floor( length(spike_rates) / (length(unique_elevation)*length(unique_fp_rotate)) ); % take minimum repetition

% creat basic matrix represents each response vector
resp = [];
for i=1:length(unique_elevation)
    for k=1:length(unique_fp_rotate)
        select = logical( (elevation==unique_elevation(i))  & (fp_rotate==unique_fp_rotate(k)) );
        for j = 1 : repetition; 
            spike_temp = spike_rates(select);   
            resp_trial{k}(j, i) = spike_temp( j );           
        end
        if (sum(select) > 0)
            resp(i, k) = mean(spike_rates(select));
            resp_err(i,k) = std(spike_rates(select)) / sqrt(repetition);
        else
            resp(i, k) = resp(1, k) ;
            resp_err(i,k) =  resp_err(1,k);
        end
    end
end

% vectorsum and calculate preferred direction
unique_azimuth_s(1:length(unique_elevation)) = 0;
for k = 1: length(unique_fp_rotate)
    [az(k), el(k), amp(k)] = vectorsumAngle(resp(:,k), unique_elevation, unique_azimuth_s );
end

%------------------------------------------------------------------
% % Define figure  X axis starts with 0-45-90-----, so change to same as 3D-tuning graph ---2/24/06 by Katsu
% xoffset=0;
% yoffset=0;
% figure(2);
% set(2,'Position', [5,15 980,650], 'Name', 'Tilt/Trans direction tuning');
% orient landscape;
% %set(0, 'DefaultAxesXTickMode', 'manual', 'DefaultAxesYTickMode', 'manual', 'DefaultAxesZTickMode', 'manual');
% 
% spon_elevation = min(unique_elevation) : 1 : max(unique_elevation);
% for k=1: length(unique_fp_rotate)     
%     if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
%         yoffset = yoffset-0.4;
%         xoffset = 0;
%     end
%     axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
%     errorbar(unique_elevation, resp(:,k), resp_err(:,k), 'bo-' );
%     hold on;
%     plot(spon_elevation, spon_resp, 'r-');
%     ylabel('spikes/s');
%     xlabel('Rot.Azimuth');
%     xlim( [min(unique_elevation), max(unique_elevation)] );    
%     title(h_title{k});
%     set(gca, 'xtick',[unique_elevation]);
%     hold on;
%     plot([az(k),az(k)],[min(resp(:,k)),max(resp(:,k))],'r-');
%     xoffset=xoffset+0.48;    
% end
% % show file name and some values in text
% axes('position',[0.05,0.85, 0.9,0.1] );
% xlim( [0,100] );
% ylim( [0,length(unique_fp_rotate)] );
% text(0, length(unique_fp_rotate), FILE);
% text(13,length(unique_fp_rotate),'prefer Rot.Azimuth       p-ANOVA             min               max');
% for k=1:length(unique_fp_rotate) 
%     % check significance by one-way-anova
%     p_1D{k} = anova1(resp_trial{k},'','off');
% % Now write
%     text(0,length(unique_fp_rotate)-k,h_title{k});
%     text(15,length(unique_fp_rotate)-k, num2str(az(k)) );
%     text(28,length(unique_fp_rotate)-k, num2str(p_1D{k}) );
%     text(40,length(unique_fp_rotate)-k, num2str( min(resp(:,k)) ) );
%     text(50,length(unique_fp_rotate)-k, num2str( max(resp(:,k)) ) );
% end
% axis off;
%------------------------------------------------------------------------------------------------
% Rotation Azimuth (270 225 180 135 90 45 0 315 270) is X axis
xoffset=0;
yoffset=0;
figure(3);
set(3,'Position', [5,15 980,650], 'Name', 'Tilt/Trans direction tuning Rot. Azimuth');
orient landscape;
%set(0, 'DefaultAxesXTickMode', 'manual', 'DefaultAxesYTickMode', 'manual', 'DefaultAxesZTickMode', 'manual');

unique_elevation_270=[270;225;180;15;90;45;0;315;270];%by KT

% spon_elevation = min(unique_elevation) : 1 : max(unique_elevation);never plot by dots I want to change to line

%re-arrangement of x-axis, start at 270;225;180;15;90;45;0;315;270
for i=1:(length(unique_elevation)+1)
    for k=1:length(unique_fp_rotate)
        if (i<8)
            resp_270(i,k)=resp(8-i,k);resp_err_270(i,k)=resp_err(8-i,k);
        elseif (i==8)
            resp_270(i,k)=resp(8,k);resp_err_270(i,k)=resp_err(8,k);
        else
            resp_270(i,k)=resp(7,k);resp_err_270(i,k)=resp_err(7,k);
        end
    end
end
%re-arrangement of az (preferreed direction)
for k=1:length(unique_fp_rotate)
    if az(k)<270
        az_270(k)=abs(az(k)-270);
        az_270(k)=((az_270(k).*8)./360)+1;% x-axis convert to 1----9
    else
        az_270(k)=360+270-az(k);
        az_270(k)=((az_270(k).*8)./360)+1;% x-axis convert to 1----9
    end
end
    
%NOW! plot errorbar and plot
for k=1: length(unique_fp_rotate)     
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.4;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
    x=[1,2,3,4,5,6,7,8,9];%temporaly x axis by number, later, change to 270-----0---270
    errorbar(x, resp_270(:,k), resp_err_270(:,k), 'bo-' );
    %errorbar(unique_elevation_270, resp_270(:,k), resp_err_270(:,k), 'bo-' );% x?? It's OK. later fix xlim
    hold on;
    plot([az_270(k),az_270(k)],[min(resp(:,k)),max(resp(:,k))],'r-');
    hold on;
    %plot(spon_elevation, spon_resp, 'r-');
    plot([1,9],[spon_resp,spon_resp],'r-');
    hold on;
    ylabel('spikes/s');
    xlabel('Rot.Azimuth');
    xlim( [1, length(unique_elevation)+1] ); 
    set(gca,'XTickMode','manual');
    set(gca,'xtick',[1,2,3,4,5,6,7,8,9]);
    set(gca,'xticklabel','270|225|180|135|90|45|0|315|270');
    title(h_title{k});
    %set(gca, 'xtick',[unique_elevation270]);
    hold on;
    xoffset=xoffset+0.48;    
end
% show file name and some values in text
axes('position',[0.05,0.85, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_fp_rotate)] );
text(0, length(unique_fp_rotate), FILE);
text(13,length(unique_fp_rotate),'spont                  min             max                        prefer Rot.Azimuth                  p-ANOVA');
for k=1:length(unique_fp_rotate) 
    % check significance by one-way-anova
    p_1D{k} = anova1(resp_trial{k},'','off');
% Now write
    text(0,length(unique_fp_rotate)-k,h_title{k});
    text(15,length(unique_fp_rotate)-k, num2str( spon_resp));
    text(27,length(unique_fp_rotate)-k, num2str( min(resp(:,k)) ) );
    text(37,length(unique_fp_rotate)-k, num2str( max(resp(:,k)) ) );
    text(52,length(unique_fp_rotate)-k, num2str(az(k)) );
    text(65,length(unique_fp_rotate)-k, num2str(p_1D{k}) );
end
axis off;
%Translation direction
axes('position',[0.59,0.05,0.4,0.05] );
xlim( [0,100] );
ylim( [0,10] );
text(0,5,'0    315     270     225     180     135      90      45      0');
text(30,2,'Translation Direction');
axis off;
%----------------------------------------------------------------------------------------------------------------
%Translation direction, so Rotation Azimuth is 180 135 90 45 0 315 270 225 180
xoffset=0;
yoffset=0;
figure(4);
set(4,'Position', [5,15 980,650], 'Name', 'Tilt/Trans direction tuning Trans. Direction');
orient landscape;
%set(0, 'DefaultAxesXTickMode', 'manual', 'DefaultAxesYTickMode', 'manual', 'DefaultAxesZTickMode', 'manual');

% unique_elevation_270=[270;225;180;15;90;45;0;315;270];%by KT

% spon_elevation = min(unique_elevation) : 1 : max(unique_elevation);never plot by dots I want to change to line

%re-arrangement of x-axis, start at 270;225;180;15;90;45;0;315;270
%As for rotation Azi. 180;135;90;45;0;315;270;225;180
for i=1:(length(unique_elevation)+1)
    for k=1:length(unique_fp_rotate)
        if (i<6)
            resp_180(i,k)=resp(6-i,k);resp_err_180(i,k)=resp_err(6-i,k);
%         elseif (i==8)
%             resp_270(i,k)=resp(8,k);resp_err_270(i,k)=resp_err(8,k);
        else
            resp_180(i,k)=resp(14-i,k);resp_err_180(i,k)=resp_err(14-i,k);
        end
    end
end
%re-arrangement of az (preferreed direction)
for k=1:length(unique_fp_rotate)
    if az(k)<180
        az_180(k)=180-az(k);
        az_180(k)=((az_180(k).*8)./360)+1;% x-axis convert to 1----9 in order to X lim 1--9
    else
        az_180(k)=360+180-az(k);
        az_180(k)=((az_180(k).*8)./360)+1;% x-axis convert to 1----9
    end
end
    
%NOW! plot errorbar and plot
for k=1: length(unique_fp_rotate)     
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.4;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
    x=[1,2,3,4,5,6,7,8,9];%temporaly x axis by number, later, change to 270-----0---270
    errorbar(x, resp_180(:,k), resp_err_180(:,k), 'bo-' );
    %errorbar(unique_elevation_270, resp_270(:,k), resp_err_270(:,k), 'bo-' );% x?? It's OK. later fix xlim
    hold on;
    plot([az_180(k),az_180(k)],[min(resp(:,k)),max(resp(:,k))],'r-');
    hold on;
    %plot(spon_elevation, spon_resp, 'r-');
    plot([1,9],[spon_resp,spon_resp],'r-');
    hold on;
    ylabel('spikes/s');
    xlabel('Translation Direction');
    xlim( [1, length(unique_elevation)+1] ); 
    set(gca,'XTickMode','manual');
    set(gca,'xtick',[1,2,3,4,5,6,7,8,9]);
    set(gca,'xticklabel','270|225|180|135|90|45|0|315|270');
    title(h_title{k});
    %set(gca, 'xtick',[unique_elevation270]);
    hold on;
    xoffset=xoffset+0.48;    
end
% show file name and some values in text
axes('position',[0.05,0.85, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_fp_rotate)] );
text(0, length(unique_fp_rotate), FILE);
text(13,length(unique_fp_rotate),'spont                  min             max                        prefer Rot.Azimuth                  p-ANOVA');
for k=1:length(unique_fp_rotate) 
    % check significance by one-way-anova
    p_1D{k} = anova1(resp_trial{k},'','off');
% Now write
    text(0,length(unique_fp_rotate)-k,h_title{k});
    text(15,length(unique_fp_rotate)-k, num2str( spon_resp));
    text(27,length(unique_fp_rotate)-k, num2str( min(resp(:,k)) ) );
    text(37,length(unique_fp_rotate)-k, num2str( max(resp(:,k)) ) );
    text(52,length(unique_fp_rotate)-k, num2str(az(k)) );
    text(65,length(unique_fp_rotate)-k, num2str(p_1D{k}) );
end
axis off;
%Rotation Azimuth
axes('position',[0.09,0.05,0.4,0.05] );
xlim( [0,100] );
ylim( [0,10] );
text(0,5,'180    135      90      45       0      315      270      225      180');
text(30,2,'Rotation Azimuth');
axis off;
%---------------------------------------------------------------------------------------
%prepare for min and max
for k=1:length(unique_fp_rotate)
    min_resp{k}=min(resp(:,k));
    max_resp{k}=max(resp(:,k));
end
%Also, write out some summary data to a cumulative summary file
buff = sprintf('%s\t %4.2f\t    %4.2f\t   %4.2f\t   %4.2f\t   %4.2f\t   %4.2f\t    %4.2f\t   %4.2f\t    %4.2f\t    %4.2f\t    %4.2f\t   %4.2f\t    %4.2f\t    %4.4f\t    %4.4f\t   %4.4f\t    %4.4f\t', ...
     FILE, spon_resp, min_resp{:}, max_resp{:}, az(:), p_1D{:} );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\TILT_TRANSLATION\Tuning2D_TT.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         Spon\t MinT+T\t MinT-T\t MinTilt\t MinTrans\t MaxT+T\t MaxT-T\t MaxTilt\t MaxTrans\t AzT+T\t AzT-T\t AzTilt\t AZTrans\t pT+T\t pT-T\t pTilt\t pTrans\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
return;