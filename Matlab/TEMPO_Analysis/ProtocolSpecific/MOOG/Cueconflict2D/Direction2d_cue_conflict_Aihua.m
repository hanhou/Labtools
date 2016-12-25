%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	GCD, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------
function Direction2d_cue_conflict_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth_moog = data.moog_params(HEADING,:,MOOG);
temp_azimuth_cam = data.moog_params(HEADING,:,CAMERAS);
temp_preferred_azimuth = data.moog_params(PREFERRED_AZIMUTH,:,MOOG);
temp_preferred_elevation = data.moog_params(PREFERRED_ELEVATION,:,CAMERAS);
preferred_azimuth = data.one_time_params(PREFERRED_AZIMUTH);
preferred_elevation = data.one_time_params(PREFERRED_ELEVATION);

temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth_moog == data.one_time_params(NULL_VALUE)) & (temp_azimuth_cam == data.one_time_params(NULL_VALUE)) );

% %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth_moog);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );

azimuth_moog = temp_azimuth_moog(~null_trials & select_trials);
azimuth_cam = temp_azimuth_cam(~null_trials & select_trials);
% elevation = temp_elevation(~null_trials & select_trials);
motion_coherence = temp_motion_coherence(~null_trials & select_trials);

spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth_moog = munique(azimuth_moog');
unique_azimuth_cam = munique(azimuth_cam');
% unique_elevation = munique(elevation');
unique_motion_coherence = munique(motion_coherence');

%% ADD CODE HERE FOR PLOTTING
% create basic matrix represents each response vector    
resp = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
resp_std = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
resp_ste = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
resp_trial = zeros( ...
    sum(select_trials) / ( length(unique_azimuth_moog)*length(unique_azimuth_cam) ), ...
    length(unique_azimuth_moog) , length(unique_azimuth_cam) );
for i=1:length(unique_azimuth_moog)
    for j=1:length(unique_azimuth_cam)
        select = logical( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(j)) );
        if (sum(select) > 0) 
            resp_trial(1:sum(select),i,j) = spike_rates(select);
            resp(i, j) = mean(spike_rates(select));
            resp_std(i,j) = std(spike_rates(select));
%             resp_ste(i,j) = resp_std(i,j) / sqrt(length(find( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(j)) )) );
            resp_ste(i,j) = resp_std(i,j) / sqrt( sum(select) );
        end
    end
end 

% the first one is the conflict dataset
resp_conflict = resp(2:length(unique_azimuth_moog), 2:length(unique_azimuth_cam) );
resp_trial_conflict = resp_trial( :, 2:length(unique_azimuth_moog), 2:length(unique_azimuth_cam) );
% Rearrange resp_trial_conflict for use with the anova2 function.
nrep_anova2 = size( resp_trial_conflict,1 );
% resp_trial_conflict_anova2 = zeros( nrep_anova2 * size( resp_trial_conflict, 3 ) , size( resp_trial_conflict, 2) );
% for i=1:nrep_anova2 % Number of reps
%     for j=1:size( resp_trial_conflict, 2) % MOOG
%         for k=1:size( resp_trial_conflict, 3 ) % cam
%             resp_trial_conflict_anova2( nrep_anova2*(k-1) + i , j ) = ...
%                 resp_trial_conflict(i,j,k);
%         end
%     end
% end
resp_trial_conflict_anova2 = transpose(reshape( permute(resp_trial_conflict,[2 1 3]), ...
    [ size(resp_trial_conflict,2) size(resp_trial_conflict,1)*size(resp_trial_conflict,3)]));

% the second is the moog control
resp_cam = resp( 1, 2 : length(unique_azimuth_cam) );
resp_cam_ste = resp_ste( 1, 2 : length(unique_azimuth_cam) );
resp_trial_cam = squeeze( resp_trial( :, 1, 2 : length(unique_azimuth_cam) ) );
% the third is the camera control
resp_ves = resp( 2:length(unique_azimuth_moog) , 1 );
resp_ves_ste = resp_ste( 2:length(unique_azimuth_moog) , 1 );
resp_trial_ves = squeeze( resp_trial( :, 2:length(unique_azimuth_moog) , 1 ) );


% Compute the p values using ANOVA.
p_anova1_ves = anova1(resp_trial_ves,'','off');
p_anova1_cam = anova1(resp_trial_cam,'','off');
p_anova1_conflict = anova1( ...
    reshape(resp_trial_conflict,...
    [size(resp_trial_conflict,1) size(resp_trial_conflict,2)*size(resp_trial_conflict,3)]) ,...
    '','off');
p_anova2_conflict = anova2( resp_trial_conflict_anova2 , nrep_anova2 , 'off' );

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));

% vector_num = length(unique_azimuth) * (length(unique_elevation)-2) + 2;
% %vector_num = 6 ;

%------------------------------------------------------------------
% Define figure

figure(2);
clf(2,'reset');
set(2,'Position', [5,15 1200,900], 'Name', '3D Direction Tuning');

axes('position',[0.05,0.3,0.6,0.55]);
%contourf( unique_azimuth_moog(2:end) , unique_azimuth_cam(2:end) , resp_conflict );
contourf(unique_azimuth_cam(2:end), ...
    unique_azimuth_moog(2:end), ...
    resp_conflict(:,:)' );
h_cont=gca;
% The transpose is because x varies with column number and y varies with
% row number. I also plotted axes x=camera and y=moog.

% xlim([0 315]);
% ylim([0 315]);
ylabel('visual')
set(gca,'xlim', [ min(unique_azimuth_moog(2:end)) max(unique_azimuth_moog(2:end)) ], ...
    'ylim', [ min(unique_azimuth_cam(2:end)) max(unique_azimuth_cam(2:end)) ],...
    'xtick', unique_azimuth_cam(2:end),...
    'xticklabel',round(unique_azimuth_cam(2:end)),...
    'ytick', unique_azimuth_moog(2:end),...
    'yticklabel',round(unique_azimuth_moog(2:end)) );
text( max(get(gca,'xlim')), max(get(gca,'ylim')),...
    sprintf('ANOVA1 p=%0.5g  ANOVA2 p_{vestibular}=%0.4g p_{visual}=%0.4g p_{interaction}=%0.4g',...
    p_anova1_conflict,p_anova2_conflict(1),p_anova2_conflict(2),p_anova2_conflict(3)),...
    'HorizontalAlignment','right','VerticalAlignment','bottom');

% set(gca,'xtick',unique_azimuth_moog(2:end),'ytick',unique_azimuth_cam(2:end));
colorbar;
%view(90,270);

% x_axis = 0:45:315;
% Along the bottom (x-axis) of the contour plot, plot the moog-only
% control.
axes('position',[0.05,0.05,0.6,0.2]);
errorbar( unique_azimuth_moog(2:end) , resp_ves(:,:), resp_ves_ste(:,:) , 'o-');
set(gca,'xlim',[ min(unique_azimuth_moog(2:end)) max(unique_azimuth_moog(2:end)) ],...
    'xtick',unique_azimuth_moog(2:end) );
xlabel('vestibular');
set(gca,'xticklabel',round(unique_azimuth_moog(2:end)));
% Set width to match contour map.
pos_ves=get(gca,'position');
pos_cont=get(h_cont,'position');
set(gca,'position', [ pos_cont(1) pos_ves(2) pos_cont(3) pos_ves(4) ]);
% errorbar(x_axis,resp_cam , resp_cam_ste,  'o-');
% set(gca, 'xtick',[0,45,90,135,180,225,270,315]);
% xlabel('visual');
% xlim([0,315]);
text( max(get(gca,'xlim')), max(get(gca,'ylim')),...
    sprintf('ANOVA p=%0.5g',p_anova1_ves),...
    'HorizontalAlignment','right','VerticalAlignment','top');

% To the right (y-axis) of the contour plot, plot the visual-only control.
axes('position',[0.7,0.3,0.2,0.55]);
errorbar( unique_azimuth_cam(2:end), resp_cam(:,:) , resp_cam_ste(:,:) , 'o-');
set(gca,'xlim',[ min(unique_azimuth_cam(2:end)) max(unique_azimuth_cam(2:end)) ], ...
    'xtick',unique_azimuth_cam(2:end),...
    'xticklabel',round(unique_azimuth_cam(2:end)),...
    'XAxisLocation','top');
xlabel('visual');
% Set position to align with contour map.
pos_cam=get(gca,'position');
pos_cont=get(h_cont,'position');
set(gca,'position', [ pos_cam(1) pos_cont(2) pos_cam(3) pos_cont(4) ]);
% errorbar(x_axis, resp_ves, resp_ves_ste, 'o-');
% xlim([0,315]);
%plot(resp_ves(:,:) , 'o-');
% set(gca, 'xtick',[0,45,90,135,180,225,270,315]);
view(90,270);
% xlabel('vestibular');
text( min(get(gca,'xlim')), max(get(gca,'ylim')),...
    sprintf('ANOVA p=%0.5g',p_anova1_cam),...
    'HorizontalAlignment','right','VerticalAlignment','top','rotation',270)

% % calculate min and max firing rate, standard deviation, DSI, Vectorsum
% Min_resp(k) = min( min( resp_mat_tran(k,:,:)) );
% Max_resp(k) = max( max( resp_mat_tran(k,:,:)) );
% resp_std(k) = sum( sum(resp_mat_std(k,:,:)) ) / vector_num;  % notice that do not use mean here, its 26 vectors intead of 40
% M=squeeze(resp_mat(k,:,:));     % notice that here DSI should use resp_temp without 0 value set manually
% DSI_temp(k) = DSI(M,spon_resp,resp_std(k));
% N=squeeze(resp_mat(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually 
% [Azi, Ele, Amp] =vectorsum(N);
% Vec_sum{k}=[Azi, Ele, Amp];

% %-------------------------------------------------------------------
% %check significance of DSI and calculate p value
% perm_num=1000;
% bin = 0.005;
% spike_rates_perm = [];
% for n=1: perm_num 
%     for k=1:length(unique_condition_num)   
%         spike_rates_pe{k} = spike_rates( find( condition_num==unique_condition_num(k) ) );
%         spike_rates_pe{k} = spike_rates_pe{k}( randperm(length(spike_rates_pe{k})) );
%         spike_rates_perm=[spike_rates_perm,spike_rates_pe{k}];            % get the permuted spikerate to re-calculate DSI for each condition
%     end
% 
%     % re-creat a matrix similar as resp_mat              
%     resp_vector_perm = [];
%     for i=1:length(unique_azimuth)
%         for j=1:length(unique_elevation)
%             for k=1:length(unique_condition_num)
%                 select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
%                 if (sum(select) > 0)
%                     resp_mat_perm(k,j,i) = mean(spike_rates_perm(select));
%                     resp_mat_perm_std(k,j,i) = std(spike_rates_perm(select));
%                 else
%                     resp_mat_perm(k,j,i) = 0;
%                     resp_mat_perm_std(k,j,i) = 0;
%                 end
%             end        
%         end
%     end
%     % re-calculate DSI now
%     for k=1: length(unique_condition_num)
%         resp_perm_std(k) = sum( sum(resp_mat_perm_std(k,:,:)) ) / vector_num; 
%         M_perm=squeeze(resp_mat_perm(k,:,:));
%         DSI_perm(k,n) = DSI(M_perm, spon_resp, resp_perm_std(k) );
%     end
%     
% end
% x_bin = 0 : bin : 1;
% for k = 1 : length(unique_condition_num)
%     histo(k,:) = hist( DSI_perm(k,:), x_bin );
%     bin_sum = 0;
%     n = 0;
%     while ( n < (DSI_temp(k)/bin) )
%           n = n+1;
%           bin_sum = bin_sum + histo(k, n);
%           p{k} = (perm_num - bin_sum)/ perm_num;    % calculate p value
%     end 
% end

%------------------------------------------------------------------

% Now show vectorsum, DSI, p and spontaneous at the top of figure
h_title{1}='conflict';
h_title{2}='vestibular';
h_title{3}='visual';
axes('position',[0.05,0.9, 0.9,0.05] );
xlim( [0,100] );
ylim( [0,3] );
h_spon = num2str(spon_resp);
text(0, 3, FILE);
text(15,3,'Coherence');

text(30,3,'Spon');
text(40,3,'pre-azi');
text(50,3,'pre-ele');

for k=1:3
    h_text{k}=num2str( [spon_resp ] );
    text(0,3-k,h_title{k});
    text(30,3-k, h_text{k} );    
end
text(15, 2, num2str(unique_motion_coherence));
text(40, 0, num2str( preferred_azimuth ));
text(50, 0, num2str( preferred_elevation ));

axis off;

%---------------------------------------------------------------------------------------
%Also, write out some summary data to a cumulative summary file

buff = sprintf('%s\t %4.2f\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t', ...
     FILE,  unique_motion_coherence, spon_resp, preferred_azimuth , preferred_elevation, p_anova1_conflict, p_anova2_conflict);
 
outfile = ['Z:\Users\Aihua\z_tempOutputs\DirectionTuningSum.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t    Coherence\t     SPon\t PrefAzi\t PrefEle\t pAnova\t pVeb\t pVis\t pInteract\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
%Save the figures
OutputPath=['Z:\Users\Aihua\z_tempOutputs\'];
FigureIndex=2;
 figure(FigureIndex); 
 set(gcf, 'PaperOrientation', 'portrait');
 saveas(gcf,[OutputPath FILE(1:end-4) '_Ch' num2str(SpikeChan) '.png'],'png');
 close(FigureIndex);
%---------------------------------------------------------------------------------------

return;