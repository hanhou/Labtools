% for RVOR/Puisuit task
%--	GY, 10/23/05   Modulated by KT 02/24/06 08/03/06
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuning2D_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
temp_fp_rotate = data.moog_params(FP_ROTATE,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             
%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_elevation);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
null_trials = logical( (temp_elevation == data.one_time_params(NULL_VALUE)) );

elevation = temp_elevation(~null_trials & select_trials);
fp_rotate = temp_fp_rotate(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_elevation = munique(elevation');
unique_fp_rotate = munique(fp_rotate');
unique_stim_type = munique(stim_type');

h_title{1,1}='Head-fixed';
h_title{2,1}='World-fixed';
h_title{3,1}='Pursuit only';
h_title{1,2}='Visual with fixation';
h_title{2,2}='Back Move pursuit';
h_title{3,2}='Back stable pursuit';
s_title{1}='Vestibular (without background dots)';
s_title{2}='Visual (without RVOR, with background dots)';

% calculate spontaneous firing rate 
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));
repetition = floor( length(spike_rates) / (length(unique_elevation)*length(unique_fp_rotate)*length(unique_stim_type)) ); % take minimum repetition

% creat basic matrix represents each response vector
resp = [];
for n=1:length(unique_stim_type)
for i=1:length(unique_elevation)
    for k=1:length(unique_fp_rotate)
        select = logical( (elevation==unique_elevation(i))  & (fp_rotate==unique_fp_rotate(k)) & (stim_type==unique_stim_type(n)));
        for j = 1 : repetition; 
            spike_temp = spike_rates(select);   
            resp_trial{k, n}(j, i) = spike_temp( j );           
        end
        if (sum(select) > 0)
            resp{n}(i, k) = mean(spike_rates(select));
            resp_err{n}(i,k) = std(spike_rates(select)) / sqrt(repetition);
        else
            resp{n}(i, k) = resp{n}(1, k) ;
            resp_err{n}(i,k) =  resp_err{n}(1,k);
        end
    end
end
end

% vectorsum and calculate preferred direction
unique_azimuth_s(1:length(unique_elevation)) = 0;
for n=1:length(unique_stim_type)
for k = 1: length(unique_fp_rotate)
    [az(k,n), el(k,n), amp(k,n)] = vectorsumAngle(resp{n}(:,k), unique_elevation, unique_azimuth_s );
end
end
%------------------------------------------------------------------
% Define figure  modified by KT new figure (3)
% xoffset=0;
% yoffset=0;
% figure(2);
% set(2,'Position', [5,15 980,650], 'Name', 'pursuit direction tuning');
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
%     xlabel('elevation');
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
% text(15,length(unique_fp_rotate),'prefer direction          p    min    max');
% for k=1:length(unique_fp_rotate) 
%     % check significance by one-way-anova
%     p_1D{k} = anova1(resp_trial{k},'','off');
%     text(0,length(unique_fp_rotate)-k,h_title{k});
%     text(15,length(unique_fp_rotate)-k, num2str(az(k)) );
%     text(25,length(unique_fp_rotate)-k, num2str(p_1D{k}) );
%     text(35,length(unique_fp_rotate)-k, num2str( min(resp(:,k)) ) );
%     text(45,length(unique_fp_rotate)-k, num2str( max(resp(:,k)) ) );
% end
% axis off;
%----------------------------------------------------------------------------------------------------
% X-axis 0-45-90-----315-0 again! Figure another graph!!
for n=1:length(unique_stim_type)
xoffset=0;
yoffset=0;
figure(n+1);
set(n+1,'Position', [5,15 980,650], 'Name', 'pursuit direction tuning');
orient landscape;
%spon_elevation = min(unique_elevation) : 1 : max(unique_elevation);never use this dot-plot method
% add 0 degrees end of x axis
for k=1:length(unique_fp_rotate)
    resp{n}(9,k)=resp{n}(1,k);% just add 0 degrees twice, so 0-45-90------315-0
    resp_err{n}(9,k)=resp_err{n}(1,k);
end
%
%To plot  %'plot([az(k,n),az(k,n)],[min(resp{n}(:,k)),max(resp{n}(:,k))],'r-');'
% az(k,n) ----> azplot(k,n)
%
%
%az(preferred direction) should be change to 1----9 x axis
% az(:,n)=((az(:,n).*8)./360)+1; %%%% az ---> azplot
azplot(:,n)=((az(:,n).*8)./360)+1;
%
%NOW! plot ; Caution X-axis is 1 2 3 4 5 6 7 8 9
for k=1: length(unique_fp_rotate)     
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.4;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
    x=[1,2,3,4,5,6,7,8,9];%temporaly x axis by number, later, change to 0-45-90----315-0
    errorbar(x, resp{n}(:,k), resp_err{n}(:,k), 'bo-' );
    hold on;
%     plot([az(k,n),az(k,n)],[min(resp{n}(:,k)),max(resp{n}(:,k))],'r-');
    plot([azplot(k,n),azplot(k,n)],[min(resp{n}(:,k)),max(resp{n}(:,k))],'r-');
    hold on;
    plot([1,9],[spon_resp,spon_resp],'r-');
    hold on;
    ylabel('Spikes/s');
    xlabel('Elevation');
    xlim( [1, length(unique_elevation)+1] ); 
    set(gca,'XTickMode','manual');
    set(gca,'xtick',[1,2,3,4,5,6,7,8,9]);
    set(gca,'xticklabel','0|45|90|135|180|225|270|315|0');
    title(h_title{k,n});
    hold on;
    xoffset=xoffset+0.48;    
end
% show file name and some values in text
axes('position',[0.05,0.85, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_fp_rotate)+1] );
text(0, length(unique_fp_rotate)+1, FILE);
text(15, length(unique_fp_rotate)+1, s_title{n});
text(11,length(unique_fp_rotate),'Prefer Direct.    p-ANOVA       spont.       min         max ');
% 
for k=1:length(unique_fp_rotate) 
    % check significance by one-way-anova
    p_1D{k,n} = anova1(resp_trial{k,n},'','off');%ANOVA p-value caluculation
    
    
    text(0,length(unique_fp_rotate)-k,h_title{k,n});
    text(14,length(unique_fp_rotate)-k, num2str(az(k,n)) );
    text(21,length(unique_fp_rotate)-k, num2str(p_1D{k,n}) );
    text(31,length(unique_fp_rotate)-k, num2str(spon_resp) );
%     text(38,length(unique_fp_rotate)-k, num2str( min(resp{n}(:,k)) ) );
%     text(45,length(unique_fp_rotate)-k, num2str( max(resp{n}(:,k)) ) );
    min_resp(k,n)=min(resp{n}(:,k));
    max_resp(k,n)=max(resp{n}(:,k));
     text(38,length(unique_fp_rotate)-k, num2str( min_resp(k,n) ) );
     text(45,length(unique_fp_rotate)-k, num2str( max_resp(k,n) ) );
end
axis off;
end
%--------------------------------------------------------------------------------------
%Also, write out some summary data to a cumulative summary file
buff = sprintf('%s\t %4.3f\t    %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t    %4.3f\t   %4.3f\t    %4.3f\t    %4.3f\t    %4.3f\t   %4.3f\t    %4.3f\t    %4.3f\t    %4.3f\t   %4.3f\t    %4.3f\t   %4.3f\t    %4.3f\t   %4.3f\t    %4.23f\t    %4.3f\t    %4.3f\t   %4.3f\t    %4.3f\t    %4.3f\t    %4.3f\t   %4.3f\t    %4.3f\t', ...
     FILE, spon_resp, min_resp(:), max_resp(:), az(:), p_1D{:} );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit\DirectionTuning2D_pursuit.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit\DirectionTuning2D_Que_pursuit.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         Spon\t Head-min\t World-min\t Pursuit-min\t Visual-min\t BackMove-min\t BackStable-min\t Head-max\t World-max\t Pursuit-max\t Visual-max\t BackMove-max\t BackStable-max\t HeadPreDir\t WorldPreDir\t PursuitPreDir\t VisualPreDir\t BackMovePreDir\t BackStablePreDir\t HeadANOVA\t WorldANOVA\t PursuitANOVA\t VisualANOVA\t BackMoveANOVA\t BackStableANOVA\t');
%     fprintf(fid, 'FILE\t         Spon\t Head-min\t World-min\t Pursuit-min\t Head-max\t World-max\t Pursuit-max\t HeadPreDir\t WorldPreDir\t PursuitPreDir\t HeadANOVA\t WorldANOVA\t PursuitANOVA\t');
fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
return;