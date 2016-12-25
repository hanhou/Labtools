%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	GCD, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------
function Direction2d_cue_conflict_fit_compare(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth_moog = data.moog_params(HEADING,:,MOOG);
temp_azimuth_cam = data.moog_params(HEADING,:,CAMERAS);
temp_preferred_azimuth = data.moog_params(PREFERRED_AZIMUTH,:,MOOG);
temp_preferred_elevation = data.moog_params(PREFERRED_ELEVATION,:,CAMERAS);
preferred_azimuth = data.one_time_params(PREFERRED_AZIMUTH);
preferred_elevation = data.one_time_params(PREFERRED_ELEVATION);

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

spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth_moog = munique(azimuth_moog');
unique_azimuth_cam = munique(azimuth_cam');
% unique_elevation = munique(elevation');


%% ADD CODE HERE FOR PLOTTING
% create basic matrix represents each response vector    
resp = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
resp_std = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
resp_ste = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
for i=1:length(unique_azimuth_moog)
    for j=1:length(unique_azimuth_cam)
        select = logical( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(j)) );
        if (sum(select) > 0) 
            resp(i, j) = mean(spike_rates(select));
            resp_std(i,j) = std(spike_rates(select));
%             resp_ste(i,j) = resp_std(i,j) / sqrt(length(find( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(j)) )) );
            resp_ste(i,j) = resp_std(i,j) / sqrt( sum(select) );
        end
    end
end 

% Grab the vestibular-only responses for ANOVA.
resp_trial_ves=zeros( sum((azimuth_moog==unique_azimuth_moog(2)) & (azimuth_cam==unique_azimuth_cam(2)) ) , ...
    length(unique_azimuth_moog)-1 );
% Loop starting from 2 to avoid the control.
for i=2:length(unique_azimuth_moog)
    % unique_azimuth_cam(1) should be the null camera case.
    select = logical( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(1)) );
    % Because looping 2:length, index with i-1.
    resp_trial_ves(:,i-1)=spike_rates(select);
end
% Compute the p value using ANOVA.
p_anova1_ves = anova1(resp_trial_ves(:,:),'','off');

% Grab the visual-only responses for ANOVA.
resp_trial_cam=zeros( sum((azimuth_moog==unique_azimuth_moog(2)) & (azimuth_cam==unique_azimuth_cam(2)) ) , ...
    length(unique_azimuth_cam)-1 );
% Loop starting from 2 to avoid the control.
for i=2:length(unique_azimuth_cam)
    % unique_azimuth_moog(1) should be the null camera case.
    select = logical( (azimuth_moog==unique_azimuth_moog(1)) & (azimuth_cam==unique_azimuth_cam(i)) );
    % Because looping 2:length, index with i-1.
    resp_trial_cam(:,i-1)=spike_rates(select);
end
% Compute the p value using ANOVA.
p_anova1_cam = anova1(resp_trial_cam(:,:),'','off');

% the first one is the conflict dataset
resp_conflict = resp(2:length(unique_azimuth_moog), 2:length(unique_azimuth_cam) );
% the second is the moog control
resp_cam = resp( 1, 2 : length(unique_azimuth_cam) );
resp_cam_ste = resp_ste( 1, 2 : length(unique_azimuth_cam) );
% the third is the camera control
resp_ves = resp( 2:length(unique_azimuth_moog) , 1 );
resp_ves_ste = resp_ste( 2:length(unique_azimuth_moog) , 1 );

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));

% vector_num = length(unique_azimuth) * (length(unique_elevation)-2) + 2;
% %vector_num = 6 ;

%%
% Fit the conflict data based on the two single cue responses. Compare
% different fitting methods.

% Subtract spontaneous rate for fitting.
conflict = resp_conflict - spon_resp;
visual = resp_cam - spon_resp;
vestibular = resp_ves - spon_resp;

% Turn the single cue data into grids for fitting the 2D conflict array.
[VES,VIS]=ndgrid( vestibular, visual );
%[VIS,VES]=ndgrid( visual , vestibular ); % This line is wrong. I'm just testing.

% Set a few fitting parameters.
A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
OPTIONS = optimset('fmincon');
OPTIONS = optimset('LargeScale', 'off', 'LevenbergMarquardt', 'on','MaxIter', 10000, 'Display', 'off');

% % One parameter fit using scaled product of visual and vestibular responses.
% fit1 = @(x) fullfit([ 0 0 0 0 x(1)]);
% error1 = @(x) sum( sum( ( fit1(x) - conflict ).^2 ) );
% es1 = [0.5]; % Initial weights
% LB1 = [-10]; % Weight lower bounds.
% UB1 = [10]; % Weight upper bounds.
% 
% weights1 = fmincon(error1,es1,A,B,Aeq,Beq,LB1,UB1, NONLCON, OPTIONS); % fminsearch
% error1_SSE = error1( weights1 );
% RMS1 = sqrt(error1_SSE) / sqrt( sum(sum( conflict.^2)) );
% prediction1 = fit1(weights1);
% 
% % Two parameter fit using weighted sum of visual and vestibular responses.
% fit2 = @(x) fullfit([ x(1) x(2) 0 0 0 ]);
% error2 = @(x) sum( sum( ( fit2(x) - conflict ).^2 ) );
% es2 = [0.5,0.5]; % Initial weights
% LB2 = [-10,-10]; % Weight lower bounds.
% UB2 = [10,10]; % Weight upper bounds.
% 
% weights2 = fmincon(error2,es2,A,B,Aeq,Beq,LB2,UB2, NONLCON, OPTIONS); % fminsearch
% error2_SSE = error2( weights2 );
% RMS2 = sqrt(error2_SSE) / sqrt( sum(sum( conflict.^2)) );
% prediction2 = fit2(weights2);
% 
% % Three parameter fit using weighted sum of visual and vestibular plus the
% % weighted product of them.
% fit3 = @(x) fullfit([ x(1) x(2) 0 0 x(3) ]);
% error3 = @(x) sum(sum( ( fit3(x) - conflict ).^2 ) );
% es3 = [0.5,0.5,0.5];  
% LB3 = [-10,-10,-10];
% UB3 = [10,10,10];
% 
% weights3 = fmincon(error3,es3,A,B,Aeq,Beq,LB3,UB3, NONLCON, OPTIONS); % fminsearch
% error3_SSE = error3( weights3 );
% RMS3 = sqrt(error3_SSE) / sqrt( sum(sum( conflict.^2)) );
% prediction3 = fit3(weights3);
% 
% % Four parameter fit using squares of single cue responses.
% fit4 = @(x) fullfit([ x(1) x(2) x(3) x(4) 0]);
% error4 = @(x)sum(sum( ( fit4(x) - conflict ).^2 ) );
% es4 = [0.5,0.5,0.5,0.5];  
% LB4 = [-10,-10,-10,-10];
% UB4 = [10,10,10,10];
% 
% weights4 = fmincon(error4,es4,A,B,Aeq,Beq,LB4,UB4, NONLCON, OPTIONS); % fminsearch
% error4_SSE = error4( weights4 );
% RMS4 = sqrt(error4_SSE) / sqrt( sum(sum( conflict.^2)) );
% prediction4 = fit4(weights4);
% 
% % Five parameter fit using squares and product of single cue responses.
% fit5 = @(x) fullfit(x);
% error5 = @(x)sum(sum( ( fit5(x) - conflict ).^2 ) );
% es5 = [0.5,0.5,0.5,0.5,0.5];  
% LB5 = [-10,-10,-10,-10,-10];
% UB5 = [10,10,10,10,10];
% 
% weights5 = fmincon(error5,es5,A,B,Aeq,Beq,LB5,UB5, NONLCON, OPTIONS); % fminsearch
% error5_SSE = error5( weights5 );
% RMS5 = sqrt(error5_SSE) / sqrt( sum(sum( conflict.^2)) );
% prediction5 = fit5(weights5);

% Define the most number of parameters to feed to the model.
maxparams = 5;

% Define the full fit function with linear, second order and interaction terms.
fullfit = @(x) x(1)*VES + x(2)*VIS + x(3)*VES.^2 + x(4)*VIS.^2 + x(5)*(VES.*VIS);

% Create an array used to test versions of this model with some
% coefficients set to zero.
fitselect = zeros(maxparams,maxparams,maxparams);
% For the 1 parameter fit, select only the cross term.
fitselect(1,5,1) = 1;
% For the 2 parameter fit, select only the two linear terms.
fitselect(1,1,2) = 1;
fitselect(2,2,2) = 1;
% For the 3 parameter fit, select the two linear terms and the cross term.
fitselect(1,1,3) = 1;
fitselect(2,2,3) = 1;
fitselect(3,5,3) = 1;
% For the 4 parameter fit, select the two linear terms and the squared terms.
fitselect(1,1,4) = 1;
fitselect(2,2,4) = 1;
fitselect(3,3,4) = 1;
fitselect(4,4,4) = 1;
% For the 5 parameter fit, select all terms.
fitselect(:,:,maxparams) = eye(maxparams);

clear weights error_SSE RMS prediction
for i=1:maxparams
    
    fit = @(x) fullfit( x * fitselect(1:i, : , i ) );
    sse = @(x) sum(sum( ( fit(x) - conflict ).^2 ) );
    es = 0.5 * ones(1,i);
    LB = -10 * ones(1,i);
    UB = 10 * ones(1,i);
    
    weights{i} = fmincon(sse,es,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS); % fminsearch
    error_SSE{i} = sse( weights{i} );
    RMS{i} = sqrt( error_SSE{i} ) / sqrt( sum(sum ( conflict.^2 ) ) );
    prediction{i} = fit( weights{i} );
    
end

% Do sequential F tests to compare fits with more parameters to the two
% parameter linear combination fit.

% I think that 64 is correct as the number of heading pairs, but I should
% double check with the bosses.

% % 3 weights versus 2 weights
% F_3_2 = ( (error2_SSE - error3_SSE) / ( length(weights3) - length(weights2) ) ) / ...
%     ( error3_SSE / ( prod(size(conflict)) - length(weights3) ) );
% P_3_2 = 1 - fcdf( F_3_2 , ( length(weights3) - length(weights2) ) , ( prod(size(conflict)) - length(weights3) ) );

% Compare each model to models with fewer coefficients.
for j=2:maxparams
    for i=1:(j-1)
                
        F{j,i} = ( ( error_SSE{i} - error_SSE{j} ) / ...
            ( ( j - i ) ) ) / ...
            ( error_SSE{j} / ( prod(size(conflict)) - j ) );
        P{j,i} = 1 - fcdf( F{j,i} , ...
            ( j - i ) , ...
            ( prod(size(conflict)) - j ) );

%         % Convert the number of parameters to strings.
%         is=num2str(i);
%         js=num2str(j);
%         % Compute the F statistic.
%         eval([ 'F_' js '_' is '=  ' ...
%             '( (error' is '_SSE - error' js '_SSE) / ' ...
%             '( length(weights' js ') - length(weights' is ') ) ) /' ...
%             '( error' js '_SSE / ( prod(size(conflict)) - length(weights' js ') ) );' ]);
%         % Get the corresponding P value for the F statistic.
%         eval([ 'P_' js '_' is '= 1 - fcdf( F_' js '_' is ' , '...
%             '( length(weights' js ') - length(weights' is ') ) , ' ...
%             '( prod(size(conflict)) - length(weights' js ') ) );' ]);
        
    end
end


%%
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
    'xtick',unique_azimuth_cam(2:end));
set(gca,'xticklabel',round(unique_azimuth_cam(2:end)));
xlabel('visual');
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
text(15,3,'Spon');
text(30,3,'pre-azi');
text(35,3,'pre-ele');
for k=1:3
    h_text{k}=num2str( [spon_resp ] );
    text(0,3-k,h_title{k});
    text(15,3-k, h_text{k} );    
end
text(30, 0, num2str( preferred_azimuth ));
text(35, 0, num2str( preferred_elevation ));

axis off;

figure(3);
subplot(2,2,1);
contourf(unique_azimuth_cam(2:end), ...
    unique_azimuth_moog(2:end), ...
    conflict' );
set(gca,'xlim', [ min(unique_azimuth_moog(2:end)) max(unique_azimuth_moog(2:end)) ], ...
    'ylim', [ min(unique_azimuth_cam(2:end)) max(unique_azimuth_cam(2:end)) ],...
    'xtick', unique_azimuth_cam(2:end),...
    'xticklabel',round(unique_azimuth_cam(2:end)),...
    'ytick', unique_azimuth_moog(2:end),...
    'yticklabel',round(unique_azimuth_moog(2:end)) );
colorbar;
title('Mean responses');

for k=1:3
    subplot(2,2,k+1);
    contourf(unique_azimuth_cam(2:end), ...
        unique_azimuth_moog(2:end), ...
        prediction{k}' );
    set(gca,'xlim', [ min(unique_azimuth_moog(2:end)) max(unique_azimuth_moog(2:end)) ], ...
        'ylim', [ min(unique_azimuth_cam(2:end)) max(unique_azimuth_cam(2:end)) ],...
        'xtick', unique_azimuth_cam(2:end),...
        'xticklabel',round(unique_azimuth_cam(2:end)),...
        'ytick', unique_azimuth_moog(2:end),...
        'yticklabel',round(unique_azimuth_moog(2:end)) );
    colorbar;
    title([num2str(k) ' parameter(s) weight(s)=' num2str(weights{k}) ' P=' num2str(cat(2,P{k,:})) ]);

end



%---------------------------------------------------------------------------------------
%Also, write out some summary data to a cumulative summary file

% buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
%      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, DSI_temp, p{:} , resp_std );
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\DirectionTuningSum.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_DSI\t Vis_DSI\t Comb_DSI\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% 
%---------------------------------------------------------------------------------------

return;