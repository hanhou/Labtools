% DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	YONG, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_3D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,StartEventBin, StopEventBin, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_spike_data = data.spike_data(SpikeChan,:);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_tri = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_tri ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_tri) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition_num = munique(condition_num');

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));

%-------------------------------------------------------------------------
%ANOVA modified by Aihua, it does not require whole trials, it does not matter if trial stopped during repetition
elevation_trials = find(unique_elevation~=-90 & unique_elevation~=90);
trials_per_rep = (length(unique_azimuth)*length(elevation_trials)+2) * length(unique_condition_num) +1;
repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);

% first parse raw data into repetitions, including null trials
for q = 1:repetitions
   azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   elevation_rep{q} = temp_elevation(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   condition_num_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
end

resp_mat_anova = [];
for k=1: length(unique_condition_num)
   clear select_rep;
   for q=1:1:repetitions
       n = 0;
       for i=1:length(unique_azimuth)
           for j=1:length(unique_elevation)
               select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==unique_elevation(j) & condition_num_rep{q}==unique_condition_num(k) );
               if (sum(select_rep{q}) > 0)
                   n = n+1;
                   resp_mat_anova{k}(q,n) = spike_rates_rep{q}(select_rep{q})';
               end
           end
       end
   end
   [p_anova, table, stats] = anova1(resp_mat_anova{k},[],'off');
   P_anova(k) = p_anova;
   anova_table{k} = table;
   F_val(k) = anova_table{k}(2,5);
end
F_val = cell2mat(F_val);

%% ADD CODE HERE FOR PLOTTING
resp_mat = [];
for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=1: length(unique_condition_num)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(spike_rates(select));
                resp_mat_vector(k, j, i) = mean(spike_rates(select)); % for vector sum only
                for t = 1 : length(spike_rates(select));              % this is to calculate response matrix based on each trial
                    spike_temp = spike_rates(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(spike_rates(select));     % calculate std between trials for later DSI usage
                resp_mat_ste(k, j, i) = resp_mat_std(k, j, i)/ sqrt(length(find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
            else
%                resp_mat_trial{k}(t, j, i) = 0;
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,1) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
                resp_mat_ste(k, j, i) = 0;
            end
        end        
    end
end

% calculate maximum and minimum firing rate
max_res = max(max(max(resp_mat)));
min_res = min(min(min(resp_mat)));

% Define figure
xoffset=0;
yoffset=0;
figure(2);
set(2,'Position', [5,15 980,650], 'Name', '3D Direction Tuning');
orient landscape;
%set(0, DefaultAxesXTickMode, 'manual', DefaultAxesYTickMode, 'manual', 'DefaultAxesZTickMode', 'manual');
axis off;
% for cosine plot
%---------YOng's Cosine Plot------------Now disable by Katsu, 05/22/06
% azi_cos = [1,2,3,4,5,6,7,8,9];
% ele_sin = [-1,-0.707,0,0.707,1];
%
for k=1: length(unique_condition_num) 
    
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.4;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
%---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
%     contourf( azi_cos, ele_sin, squeeze( resp_mat_tran(k,:,:)) );
%--------- Old fashioned sountourf same as Rotation 3D -------------------
    contourf( squeeze( resp_mat(k,:,:)) );% temporaly commented by
    colorbar;
    % make 0 correspond to rightward and 180 correspond to leftward
    set(gca, 'ydir' , 'reverse');
    set(gca, 'xtick', [] );
    set(gca, 'ytick', [] );    
    title( num2str(unique_stim_type(k)) );

    % plot 1-D for mean respond as a function of elevation
    % notice that elevation scale is transformed by consine
    axes('position',[0.06+xoffset 0.54+yoffset 0.04 0.24]);
    for j=1:length(unique_elevation)
        y_elevation_mean(1,j)=mean(resp_mat(k,j,:));
        y_elevation_std(1,j) =std( spike_rates([find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )]) );
        y_elevation_ste(1,j) =y_elevation_std(1,j)/ sqrt(length(find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
    end
 %---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
%    x_elevation=[-1,-0.707,0,0.707,1];
    x_elevation=unique_elevation;%05/22/06 Katsu changed not cosin axis
    errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'ko-');%-----------Temporaly disable
    xlabel('Elevation');
    view(90,90);
    set(gca, 'xtick',x_elevation);
    xlim([-90, 90]);%12/20 Katsu changed, not cosine axis
    % xlim([-1, 1]);%---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
    ylim([min(y_elevation_mean(1,:))-max(y_elevation_ste(1,:)), max(y_elevation_mean(1,:))+max(y_elevation_ste(1,:))]); %----------Now disable temporaly

    % plot 1-D for mean respond as a function of azimuth
    axes('position',[0.11+xoffset 0.46+yoffset-0.2 0.274 0.2]);
    for i=1:(length(unique_azimuth) )
    %    y_azimuth_mean(1,i)=mean(resp_mat_tran(k,:,i));
        y_azimuth_mean(1,i)=resp_mat(k,2,i);
        y_azimuth_std(1,i) =std( spike_rates([find( (azimuth==unique_azimuth(i))&(elevation==unique_elevation(2))&(condition_num==unique_condition_num(k)) )]) );
        y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (azimuth==unique_azimuth(i))&(elevation==unique_elevation(2))&(condition_num==unique_condition_num(k)) )) );    
    end
    for i=1:( length(unique_azimuth))  
        y_azimuth_ste(1,i) = y_azimuth_ste(1,i);
    end
    x_azimuth=unique_azimuth;
    errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste,'ko-');%----------------temporaly disable
    xlim([0,unique_azimuth(end)]);
    hold on;
    plot([0,360],[spon_resp,spon_resp],'--');
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',unique_azimuth);     
    xlabel('Azimuth');
    ylim([min(y_azimuth_mean(1,:))-max(y_azimuth_ste(1,:)), max(y_azimuth_mean(1,:))+max(y_azimuth_ste(1,:))]); %----------Now disable temporaly

    xoffset=xoffset+0.48;
    
    % calculate min and max firing rate, standard deviation, HTI, Vectorsum
    Min_resp(k) = min( min( resp_mat(k,:,:)) );
    Max_resp(k) = max( max( resp_mat(k,:,:)) );
    M=squeeze(resp_mat_vector(k,:,:));     % notice that here DSI should use resp_temp without 0 value set manually
    % this part is to calculate vestibular gain
%    resp_onedim{k} = [M(1,1),M(2,:),M(3,:),M(4,:),M(5,1)]';     % hard-code temperarilly    
    N=squeeze(resp_mat_vector(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually 
    [Azi, Ele, Amp] = vectorsum(N);
    Vec_sum{k}=[Azi, Ele, Amp];
    % Heading Tuning Index
%     r(k) = HTI(M,spon_resp);   % call HTI function    
end

% % %-------------------------------------------------------------------
% %check significance of HTI and calculate p value, do bootstrap at the same time to test value varience
% perm_num=1000;
% bin = 0.005;
% spike_rates_perm = [];
% for n=1: perm_num
%     % this is permuted based on trials
%     for k=1:length(unique_condition_num)   
%         spike_rates_pe{k} = spike_rates( find( condition_num==unique_condition_num(k) ) );
%         spike_rates_pe{k} = spike_rates_pe{k}( randperm(length(spike_rates_pe{k})) );
%     end
% 
%     % put permuted data back to spike_rates
%     spike_rates_perm(length(spike_rates))=0;
%     for k=1:length(unique_condition_num) 
%         ii = find(stim_type == unique_stim_type(k));
%         spike_rates_perm(ii) = spike_rates_pe{k};
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
%     
%     % re-calculate HTI now
%     for k=1: length(unique_condition_num)
%  %       resp_perm_std(k) = sum( sum(resp_mat_perm_std(k,:,:)) ) / vector_num; 
%         M_perm=squeeze(resp_mat_perm(k,:,:));
%         r_perm(k,n) = HTI(M_perm, spon_resp);                                  
%     end
%     % do bootstrap now
%     % first to bootstap raw data among trils 
%     repetition = 5;
%     for k=1:length(unique_stim_type)
%         for i=1:length(unique_azimuth)
%             for j=1:length(unique_elevation)
%                     select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (stim_type==unique_stim_type(k)) );
%                     if (sum(select) > 0)
%                         spike_select = spike_rates(select);
%                         for b=1:repetition    % use 5 repetitions temporarilly, should doesn't matter whether have 5 repetitions or not actually
%                             spike_select = spike_select( randperm(length(spike_select)) );
%                             spike_bootstrap(b) = spike_select(1);   % always take the first one element
%                         end 
%                         resp_mat_boot(k, j, i) = mean(spike_bootstrap);
%                     else
%                         resp_mat_boot(k, j, i) = 0;
%                     end       
%              end
%          end
%      end
%      % now recalculate values
%      for k=1: length(unique_stim_type)   
%          Mb=squeeze(resp_mat_boot(k,:,:));
%          r_boot(k,n) = HTI(Mb,spon_resp) ; 
%          % also calculate the significant different angles between preferred headings, for only vestibualr condition
%          [Azi_boot, Ele_boot, Amp_boot] = vectorsum(Mb);
%          Vec_sum_boot{k}=[Azi_boot, Ele_boot, Amp_boot];
%          Angle_boot(n)=(180/3.14159) * acos( sin(Vec_sum{1}(2)*3.14159/180) * sin(Vec_sum_boot{1}(2)*3.14159/180)  +  cos(Vec_sum_boot{1}(2)*3.14159/180) * sin(Vec_sum_boot{1}(1)*3.14159/180) * cos(Vec_sum{1}(2)*3.14159/180) * sin(Vec_sum{1}(1)*3.14159/180) + cos(Vec_sum_boot{1}(2)*3.14159/180) * cos(Vec_sum_boot{1}(1)*3.14159/180) * cos(Vec_sum{1}(2)*3.14159/180) * cos(Vec_sum{1}(1)*3.14159/180) );
%      end  
% end
% % now calculate p value or significant test
% x_bin = 0 : bin : 1;
% for k = 1 : length(unique_condition_num)
%     hist_perm(k,:) = hist( r_perm(k,:), x_bin );  % for permutation
%     hist_boot(k,:) = hist( r_boot(k,:), x_bin );  % for bootstrap
%     [hist_boot_angle(k,:),x_angle] = hist( Angle_boot(:), 200 );  % for bootstrap, set 200 bins temporarilly
%     bin_sum = 0;
%     n = 0;
%     while ( n < (r(k)/bin) )
%           n = n+1;
%           bin_sum = bin_sum + hist_perm(k, n);
%           p{k} = (perm_num - bin_sum)/ perm_num;    % calculate p value for HTI
%     end 
%     bin_sum = 0;
%     n = 0;
%     while ( bin_sum < 0.025*sum( hist_boot(k,:)) )   % define confidential value to be 0.05, now consider one side only which is 0.025 of confidence
%           n = n+1;
%           bin_sum = bin_sum + hist_boot(k, n);      
%           HTI_boot(k) = r(k) - n * bin ;    % calculate what HTI value is thought to be significant different
%     end 
% %     bin_sum = 0;
% %     n = 0;
% %     while ( bin_sum < 0.975*sum( hist_boot_angle(k,:)) )   % define confidential value to be 0.05, now consider one side only which is 0.025 of confidence
% %           n = n+1;
% %           bin_sum = bin_sum + hist_boot_angle(k, n);      
% %           Angle_boot(k) = x_angle(n) ;    
% %     end 
% end
%%------------------------------------------------------------------
% DDI Direction discrimination index   by Katsu 05/18/06
%--------------------------------------------------------------------
%spike_rates_sqrt = sqrt(spike_rates);That is not original, it makes value bigger by 
each_stim_trials=repetitions*(length(unique_azimuth)*length(elevation_trials)+2) % 26 trajectory  -90+90 -45 0 + 45 *0 45 90 .......
%
SSE_term = [];
for k=1:length(unique_stim_type)
    n=0;
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)
            clear select;
            select=logical((azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (stim_type==unique_stim_type(k)));
            if (sum(select) > 0)
                   n = n+1;
                   SSE_term(k,j,i)=sum((spike_rates(select)-mean(spike_rates(select))).^2);
               else
                   SSE_term(k,j,i)=0;
            end
         end
     end
%      SSE_azimth_sum(k,j)=sum(SSE_term(k,j,:))  %it is not correct
%      SSE_total(k)=sum(SSE_azimth_sum(k,:))
     SSE_total(k)=sum(sum(SSE_term(k,:,:)));
     max_min_term(k)=(Max_resp(k)-Min_resp(k))/2;
     var_term(k)=sqrt(SSE_total(k)/(each_stim_trials-n));
     DDI(k)=max_min_term(k)/(max_min_term(k)+var_term(k))
 end
 %----------------------------------------------------------------------------
% Now show vectorsum, DSI, p and spontaneous at the top of figure
axes('position',[0.05,0.85, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_condition_num)] );
h_spon = num2str(spon_resp);
text(0, length(unique_condition_num), FILE);
%text(15,length(unique_condition_num),'Spon            Minimum        Maximum       Azi             Ele                Amp           Std             HTI             HTIerr             p');
text(10,length(unique_condition_num),'Protocol          Spon    Minimum   Maximum    Azi      Ele       Amp           F-val        p-ANOVA          DDI');
for k=1:length(unique_condition_num) 
    h_text{k}=num2str( [spon_resp, Min_resp(k), Max_resp(k), Vec_sum{k}, F_val(k), P_anova(k), DDI(k)]);
    text(0,length(unique_condition_num)-k,num2str(unique_stim_type(k)));
    text(10,length(unique_condition_num)-k,'Translation');
    text(20,length(unique_condition_num)-k, h_text{k} );
end

axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSTH
% add parameters here
% timebin for plot PSTH
timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_rates)*5000/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% find spontaneous trials which azimuth,elevation,stim_type=-9999
spon_found = find(null_trials==1);    

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(null_trials==1 | temp_spike_rates > 3000 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;
for k=1: length(unique_condition_num)
    for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );            
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                resp{k}(j,i) = mean(spike_rates(select));
                act_found = find( select==1 );
                % count spikes per timebin on every same condition trials
                for repeat=1:length(act_found) 
                    for n=1:(x_length)
                        temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                        time_step=time_step+timebin;
                    end
                    time_step=1;                    
                end
                count_y_trial{k,i,j}(:,:) = temp_count;  % each trial's PSTH 
                % get the average of the total same conditions if repetion is > 1
           %     if (length(act_found) > 1);
                dim=size(temp_count);
                if dim(1) > 1;
                   count_y{i,j,k} = mean(temp_count);
                else
                   count_y{i,j,k}= temp_count;     % for only one repetition cases
                end
               
             else
                resp{k}(j,i) = 0; 
                count_y{i,j,k}=count_y{1,j,k};
             end   
             % normalize count_y
             if max(count_y{i,j,k})~=0;
                count_y_norm{i,j,k}=count_y{i,j,k} / max(count_y{i,j,k});
             else
                count_y_norm{i,j,k}=0;
             end
        end
    end  
    % now find the peak
    [row_max, col_max] = find( resp{k}(:,:)==max(max(resp{k}(:,:))) );
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_m{k}=row_max(1);
    col_m{k}=col_max(1);
    if max(count_y{col_max(1), row_max(1), k})~=0;
       count_y_max{k} = count_y{col_max(1), row_max(1), k} / max(count_y{col_max(1), row_max(1), k});
    else
       count_y_max{k} =0;
    end
    % find the largest y to set scale later
    if max(count_y{col_max(1), row_max(1), k}) > max_count
        max_count = max(count_y{col_max(1), row_max(1), k});
    end
end

% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
% define figure
figure(3);
set(3,'Position', [5,5 1000,700], 'Name', '3D Direction Tuning');
orient landscape;
title(FILE);
axis off;

xoffset=0;
yoffset=0;

% now plot
for k=1: length(unique_condition_num) 
    
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.42;
        xoffset = 0;
    end
    % output some text 
    axes('position',[0 0 1 0.9]); 
    xlim([-50,50]);
    ylim([-50,50]);
    text(-30+xoffset*100,52+yoffset*110, num2str(unique_stim_type(k)) );
    text(25,-40, 'Translation');
    axis off;
    hold on;
    for i=1:length(unique_azimuth)                    % aizmuth 270 are plotted two times in order to make circular data
        for j=1:length(unique_elevation)
            axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]);                              
            bar( x_time,count_y{i,j,k}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
            hold on;
            plot( x_start, y_marker, 'r-');
            plot( x_stop,  y_marker, 'r-');
            set( gca, 'xticklabel', ' ' );            
            % set the same scale for all plot
            xlim([0,x_length]);
            ylim([0,max_count]);
        end   
        xlabel(num2str(unique_azimuth(i)));
    end 

    xoffset=xoffset+0.46;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Also, write out some summary data to a cumulative summary file
% buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
%      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std,  F_val, P_anova, DDI );
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Direction_Que_Labyrinthectomy_DDI_Katsu.dat'];%051906 Katsu for DDI
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
% %     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t');
% fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);

return;