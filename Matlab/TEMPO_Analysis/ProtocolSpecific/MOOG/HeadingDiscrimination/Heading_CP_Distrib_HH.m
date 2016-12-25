%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04
%-----------------------------------------------------------------------------------------------------------------------

function Heading_CP_Distrib_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_total_trials = data.misc_params(OUTCOME, :);

trials = 1:length(temp_heading);		% a vector of trial indices
%bad_tri = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
% if ( bad_tri ~= NaN)
%    select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_tri) );
% else 
%    select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
% end
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
spike_rates = temp_spike_rates( select_trials);
total_trials = temp_total_trials( select_trials);
unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
%disc_heading = unique_heading( (length(unique_heading)/2+1) : end );

% divide data into two groups according to monkey's choice
LEFT = 1;
RIGHT = 2;
for i=1 : length(spike_rates)
    temp = data.event_data(1,:,i+BegTrial-1);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = RIGHT;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = LEFT;
    else
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end
end
%choice(889) =2; % for cell m2c384r2

repetition = floor( length(spike_rates)/length(unique_heading)/length(unique_stim_type) ); % take minimum repetition
resp_heading = [];
Z_Spikes = spike_rates;
% z-score data for later cp analysis across headings
for k = 1:length(unique_stim_type)
    for i = 1:length(unique_heading)
        select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
        z_dist = spike_rates(select);
        z_dist = (z_dist - mean(z_dist))/std(z_dist);
        Z_Spikes(select) = z_dist;
    end
end
% now group neuronal data into two groups according to monkey's choice
for k = 1:length(unique_stim_type)    % notice that the condition is double than disc_heading     
    for i = 1:length(unique_heading)
        select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
        resp{k,i} = spike_rates(select); 
        % calculate firing rate of each trial
        for j = 1 : repetition; 
            spike_temp = spike_rates(select);   
            resp_heading{k}(i, j) = spike_temp( j );           
        end
        resp_mat{k}(i) = mean(resp{k,i});  % the mean firing rate for each heading 
        resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);
        resp_left_choose{k,i} = spike_rates(select & (choice == LEFT) );
        resp_right_choose{k,i} = spike_rates(select & (choice == RIGHT) );    
        if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)
      %  if (length(resp_left_choose{k,i}) / length(resp{k,i}) <0.25) |  (length(resp_left_choose{k,i}) / length(resp{k,i}) >0.75)   % make sure each condition has at least 3 data values   
            CP{k}(i) = 9999;
            Z_Spikes(select) = 9999;
        else
            CP{k}(i) = 0;
        end 
    end  
    
    % now across all data according to pre and null
    % decide preferred direction first
    [rr,pp] = corrcoef(unique_heading, resp_mat{k}(:)); % if rr>0, right is pre, vice versa
    line_re{k} = rr(1,2);
    line_p{k} = pp(1,2);
    %left_select = logical( (heading==unique_heading(1)) | (heading==unique_heading(2)) |(heading==unique_heading(3))|(heading==unique_heading(4)) );
    %right_select = logical( (heading==unique_heading(6)) | (heading==unique_heading(7)) |(heading==unique_heading(8))|(heading==unique_heading(9)) );
    left_select = logical(heading > 0);
    right_select = logical(heading < 0);
        
    if line_re{k} > 0
        resp_pre_left_all{k} = Z_Spikes( right_select & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
        resp_pre_right_all{k} = Z_Spikes( right_select & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
        resp_null_left_all{k} = Z_Spikes( left_select & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
        resp_null_right_all{k} = Z_Spikes( left_select & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
    else
        resp_pre_left_all{k} = Z_Spikes( left_select & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
        resp_pre_right_all{k} = Z_Spikes( left_select & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
        resp_null_left_all{k} = Z_Spikes( right_select & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
        resp_null_right_all{k} = Z_Spikes( right_select & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
    end
    
end

%now calculate propotion correct from area under ROC curves
for k = 1 : length(unique_stim_type)
    for i = 1 : length(unique_heading)
        if CP{k}(i)~=9999
           CP{k}(i) = rocN( resp_left_choose{k,i},resp_right_choose{k,i},100 );
        else
           CP{k}(i) = NaN;
        end
        if  line_re{k} > 0  
            CP{k}(i) = 1 - CP{k}(i);
        end
    end
    
     % prefer versus null
     if (length(resp_pre_left_all{k}) > 3) & (length(resp_pre_right_all{k}) > 3)
         CP_pre{k} = rocN( resp_pre_left_all{k},resp_pre_right_all{k},100 );         
     else
         CP_pre{k} = NaN;
     end
     if (length(resp_null_left_all{k}) > 3) & (length(resp_null_right_all{k}) > 3)
         CP_null{k} = rocN( resp_null_left_all{k},resp_null_right_all{k},100 );
     else
         CP_null{k} = NaN;
     end
     if  line_re{k} > 0  
         CP_pre{k} = 1 - CP_pre{k};
         CP_null{k} = 1 - CP_null{k};
     end
end

% now plot neuronal respose distribution according to monkey's choice
num_count_left = [];
num_count_right = [];
for k = 1:length(unique_stim_type)
    xbin_low{k} = min( min(resp_heading{k}(:, :)) );   % set a xbin range
    xbin_high{k} =  max( max(resp_heading{k}(:, :)) );
    xbin_step{k} = floor( ( max(max(resp_heading{k}(:, :)))-min(min(resp_heading{k}(:, :))) ) / 10 );   % temporarilly take 1/10 step of the total range;
    if xbin_step{k} ==0
        xbin_step{k} = 1;
    end
    xbin{k} = xbin_low{k} : xbin_step{k} : xbin_high{k};
    xbin_tick{k} = floor(xbin_low{k}) : (xbin_step{k}*2) : floor(xbin_high{k});  % for later x label use
    xbin{k}
    for i = 1:length(unique_heading)
        [n_left,count_left] = hist( resp_left_choose{k,i}, xbin{k} );
        [n_right,count_right] = hist( resp_right_choose{k,i}, xbin{k} );
        num_count_left{k,i} = n_left;
        num_count_right{k,i} = n_right;
    end
end
% now combined data set with left and right choice
for k = 1:length(unique_stim_type)
    for i = 1 : length(unique_heading)
        try
            %num_count{k,i} =  [num_count_left{k,i}', num_count_right{k,i}'];    
            num_count{k,i} =  [num_count_left{k,i}(:), num_count_right{k,i}(:)];    % HH20130824
        catch 
            keyboard
        end
    end  
end

% %--------------------------------------------------
% plot distribution here
h{1} = 'vestibular';
h{2} = 'visual';
h{3} = 'combined';
legend_txt{1}='c-left';
legend_txt{2}='c-right';
figure(4);
set(4,'Position', [5,25, 980,650], 'Name', 'psycho_neurometic function');
orient landscape;
xoffset=0;
yoffset=0;
for k = 1:length(unique_stim_type) 
    for i = 1 : (length(unique_heading)+1)/2
        axes('position',[0.05+xoffset 0.84-yoffset, 0.13 0.10]);    % plot the left heading data set
%         bar( xbin{k}, num_count{k,i},'group');
        bar( xbin{k}, num_count{k,i},1, 'group');colormap Gray;
        xlim([min( min(resp_heading{k}(:, :)) )-2, max( max(resp_heading{k}(:, :)))+2 ]);
        ylim([0,max(max(num_count{k,i}))+1]);
        if i==1
           title(h{k});
        end
        if i==5 & k==1
           xlabel('spikes/s');
        end
        if k==1
           ylabel( num2str(unique_heading(i)) );
        end 
        if i==5
            set(gca, 'xtick', xbin_tick{k});
        else
            set(gca, 'xtick', [] );
        end
          
        axes('position',[0.2+xoffset 0.84-yoffset, 0.12 0.10]);     % plot the right heading data set
%         bar( xbin{k}, num_count{k,length(unique_heading)+1-i},'group');
        bar( xbin{k}, num_count{k,length(unique_heading)+1-i},1, 'group'); colormap Gray;  
        xlim([min( min(resp_heading{k}(:, :)) )-2, max( max(resp_heading{k}(:, :)))+2 ]);
        ylim([0,max(max(num_count{k,length(unique_heading)+1-i}))+1]);
        if k==length(unique_stim_type) & i==1
            legend(legend_txt{:},0);
        end
        if i==5
            set(gca, 'xtick', xbin_tick{k});
        else
            set(gca, 'xtick', [] );
        end
        
        % show CP in the figure
        axes('position',[0.05+xoffset 0.8-yoffset, 0.13 0.14]); 
        plot(0,0);
        xlim([-10,10]);
        ylim([-10,10]);
        text(2, 6, num2str(CP{k}(i)) );
        axis off;
        
        axes('position',[0.2+xoffset 0.8-yoffset, 0.12 0.14]);  
        plot(0,0);
        xlim([-10,10]);
        ylim([-10,10]);
        text(2, 6, num2str(CP{k}(length(unique_heading)+1-i)) );
        axis off;
        
        yoffset=yoffset+0.14;
    end
    yoffset=0;
    xoffset=xoffset+0.32;
    
end
axes('position',[0.02 0.95, 0.1 0.1]);
xlim([1,10]);
ylim([0,1]);
text(0,0.2,FILE);
axis off;
%---------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 50 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
buff = sprintf(sprint_txt, FILE, CP{1}(:), CP_pre{1}, CP_null{1} );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum_CPdist.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         CP\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%--------------------------------------------------------------------------

return;