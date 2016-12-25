function P_anova_hor=DirectionTuningPlot_3D_Supple(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);   
%DirectionTuningPlot_3D_Supple: 3D Translation add the plus/minus 22.5 deg in the horizontal plane. 

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);

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
spon_resp = mean(temp_spike_rates(spon_found))
% added by Katsu 111606
spon_std = std(temp_spike_rates(spon_found))
%-------------------------------------------------------------------------
%% ADD CODE HERE FOR PLOTTING
resp_mat = [];
for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=1: length(unique_condition_num)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
%             select2 = logical( (azimuth==unique_azimuth(i)) & (elevation==0) & (condition_num==unique_condition_num(k)) );
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
                resp_mat_vector(k,j,i) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
                resp_mat_ste(k, j, i) = 0;
            end
        end        
    end
end

%% creat a real 3-D based plot where the center correspond to forward and both lateral edges correspond to backward
%%%% Usually, axis azimuth from left is 270-225-180-135-90--0---90 %%%% 
unique_azimuth_s=[0 45 90 135 180 225 270 315];
unique_azimuth_plot=[270 225 180 135 112.5 90 67.5 45 0 315 270];%unique_azimuth_plot=[270 225 180 135 90 45 0 315 270];
for i=1:length(unique_azimuth_plot)
    Index(i)=find(unique_azimuth==unique_azimuth_plot(i));
    resp_mat_tran(:,:,i) = resp_mat(:,:,Index(i));
end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate maximum and minimum firing rate
max_res = max(max(max(resp_mat)));
% max_vi_ve = max(max(max(resp_mat(2:3,:,:))));
min_res = min(min(min(resp_mat_tran)));

vector_num = length(unique_azimuth) * (length(unique_elevation)-2) + 2;
repeat = floor( length(spike_rates) / vector_num );

% Define figure
xoffset=0;
yoffset=0;
figure(2);
set(2,'Position', [5,15 980,650], 'Name', '3D Direction Tuning');
orient landscape;
axis off;

for k=1: length(unique_condition_num) 
    
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.4;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
%---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
%     contourf( azi_cos, ele_sin, squeeze( resp_mat_tran(k,:,:)) );
%--------- Old fashioned sountourf same as Rotation 3D -------------------
    contourf( squeeze( resp_mat_tran(k,:,:)) ); % commented by Katsu
    % set the same scale for visual and combined conditions but here assuming vestibular response is always smaller than that in visual and
    % combined conditions
    colorbar;
    % make 0 correspond to rightward and 180 correspond to leftward
    set(gca, 'ydir' , 'reverse');
    set(gca, 'xtick', [] );
    set(gca, 'ytick', [] );    
    title( h_title{k} );

    % plot 1-D for mean respond as a function of elevation
    % notice that elevation scale is transformed by consine
    axes('position',[0.06+xoffset 0.54+yoffset 0.04 0.24]);
%     set(gca, 'xtick', [] );
%     set(gca, 'ytick', [] );  
    for j=1:length(unique_elevation)
        y_elevation_mean(1,j)=mean(resp_mat_tran(k,j,:));
        y_elevation_std(1,j) =std( spike_rates([find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )]) );
        y_elevation_ste(1,j) =y_elevation_std(1,j)/ sqrt(length(find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
    end
 %---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
%     x_elevation=[-1,-0.707,0,0.707,1];
    x_elevation=unique_elevation;%05/22/06 Katsu changed not cosin axis
    errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'ko-');%-----------Temporaly disable
%     errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'k-');% Katsu for paper settle for m3c294
    xlabel('Elevation');
     view(90,90);
    set(gca, 'xtick',x_elevation);
    xlim([-90, 90]);%12/20 Katsu changed, not cosine axis
%     xlim([-1, 1]);%---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
%     xlim([-1.1, 1.1]);%---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
    ylim([min(y_elevation_mean(1,:))-max(y_elevation_ste(1,:)), max(y_elevation_mean(1,:))+max(y_elevation_ste(1,:))]);%axis off %----------Now add axis off
%


    % plot 1-D for mean respond as a function of azimuth
    axes('position',[0.11+xoffset 0.46+yoffset 0.274 0.06]);
    for i=1:(length(unique_azimuth_s))
        y_azimuth_mean(1,i)=mean(resp_mat_tran(k,:,i));
        y_azimuth_std(1,i) =std( spike_rates([find( (azimuth==unique_azimuth_s(i))&(condition_num==unique_condition_num(k)) )]) );
        y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (azimuth==unique_azimuth_s(i))&(condition_num==unique_condition_num(k)) )) );    
    end
    y_azimuth_mean(1,9) = mean(resp_mat_tran(k,:,1));
    
    for i=1:( length(unique_azimuth_s)+1 )
        if (i < 8)        
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8-i);
        elseif (i == 8)
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8);
        else
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,7);
        end
    end

    x_azimuth=1:(length(unique_azimuth_s)+1);
    errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'ko-');%----------------temporaly disable
%     errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'k-');% Katsu for paper settle for m3c294
%     xlim( [1, length(unique_azimuth)+1] );
    xlim( [0.9, length(unique_azimuth_s)+1.1] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); % Katsu
%     set(gca, 'xticklabel','0|45|90|135|180|225|270|315|360'); 
    xlabel('Azimuth');
    ylim([min(y_azimuth_mean(1,:))-max(y_azimuth_ste(1,:)), max(y_azimuth_mean(1,:))+max(y_azimuth_ste(1,:))]);%axis off %----------Now add axis off

    xoffset=xoffset+0.48;
    
    % calculate min and max firing rate, standard deviation, HTI, Vectorsum
    Min_resp(k) = min( min( resp_mat_tran(k,:,:)) );
    Max_resp(k) = max( max( resp_mat_tran(k,:,:)) );
    resp_std(k) = sum( sum(resp_mat_std(k,:,:)) ) / vector_num;  % notice that do not use mean here, its 26 vectors intead of 40
    M=squeeze(resp_mat(k,:,:));     % notice that here DSI should use resp_temp without 0 value set manually
    % this part is to calculate vestibular gain
    resp_onedim{k} = [M(1,1),M(2,:),M(3,:),M(4,:),M(5,1)]';     % hard-code temperarilly    
    N=squeeze(resp_mat_vector(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually 
    [Azi, Ele, Amp] = vectorsum(N);
    Vec_sum{k}=[Azi, Ele, Amp];
    % Heading Tuning Index
    r(k) = HTI(M,spon_resp);   % call HTI function    
end

% %-------------------------------------------------------------------
%check significance of HTI and calculate p value, do bootstrap at the same time to test value varience
perm_num=1000;
bin = 0.005;
spike_rates_perm = [];
for n=1: perm_num
    % this is permuted based on trials
    for k=1:length(unique_condition_num)   
        spike_rates_pe{k} = spike_rates( find( condition_num==unique_condition_num(k) ) );
        spike_rates_pe{k} = spike_rates_pe{k}( randperm(length(spike_rates_pe{k})) );
    end

    % put permuted data back to spike_rates
    spike_rates_perm(length(spike_rates))=0;
    for k=1:length(unique_condition_num) 
        ii = find(stim_type == unique_stim_type(k));
        spike_rates_perm(ii) = spike_rates_pe{k};
    end
    
    % re-creat a matrix similar as resp_mat              
    resp_vector_perm = [];
    for i=1:length(unique_azimuth_s)
        for j=1:length(unique_elevation)
            for k=1:length(unique_condition_num)
                select = logical((azimuth==unique_azimuth_s(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
                if (sum(select) > 0)
                    resp_mat_perm(k,j,i) = mean(spike_rates_perm(select));
                    resp_mat_perm_std(k,j,i) = std(spike_rates_perm(select));
                else
                    resp_mat_perm(k,j,i) = 0;
                    resp_mat_perm_std(k,j,i) = 0;
                end
            end        
        end
    end
    
    % re-calculate HTI now
    for k=1: length(unique_condition_num)
 %       resp_perm_std(k) = sum( sum(resp_mat_perm_std(k,:,:)) ) / vector_num; 
        M_perm=squeeze(resp_mat_perm(k,:,:));
        r_perm(k,n) = HTI(M_perm, spon_resp);                                  
    end
    % do bootstrap now
    % first to bootstap raw data among trils 
    repetition = 5;
    for k=1:length(unique_stim_type)
        for i=1:length(unique_azimuth)
            for j=1:length(unique_elevation)
                    select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (stim_type==unique_stim_type(k)) );
                    if (sum(select) > 0)
                        spike_select = spike_rates(select);
                        for b=1:repetition    % use 5 repetitions temporarilly, should doesn't matter whether have 5 repetitions or not actually
                            spike_select = spike_select( randperm(length(spike_select)) );
                            spike_bootstrap(b) = spike_select(1);   % always take the first one element
                        end 
                        resp_mat_boot(k, j, i) = mean(spike_bootstrap);
                    else
                        resp_mat_boot(k, j, i) = 0;
                    end       
             end
         end
     end
     % now recalculate values
     for k=1: length(unique_stim_type)   
         Mb=squeeze(resp_mat_boot(k,:,:));
         r_boot(k,n) = HTI(Mb,spon_resp) ; 
         % also calculate the significant different angles between preferred headings, for only vestibualr condition
         [Azi_boot, Ele_boot, Amp_boot] = vectorsum(Mb);
         Vec_sum_boot{k}=[Azi_boot, Ele_boot, Amp_boot];
         Angle_boot(n)=(180/3.14159) * acos( sin(Vec_sum{1}(2)*3.14159/180) * sin(Vec_sum_boot{1}(2)*3.14159/180)  +  cos(Vec_sum_boot{1}(2)*3.14159/180) * sin(Vec_sum_boot{1}(1)*3.14159/180) * cos(Vec_sum{1}(2)*3.14159/180) * sin(Vec_sum{1}(1)*3.14159/180) + cos(Vec_sum_boot{1}(2)*3.14159/180) * cos(Vec_sum_boot{1}(1)*3.14159/180) * cos(Vec_sum{1}(2)*3.14159/180) * cos(Vec_sum{1}(1)*3.14159/180) );
     end  
end
% now calculate p value or significant test
x_bin = 0 : bin : 1;
for k = 1 : length(unique_condition_num)
    hist_perm(k,:) = hist( r_perm(k,:), x_bin );  % for permutation
    hist_boot(k,:) = hist( r_boot(k,:), x_bin );  % for bootstrap
    [hist_boot_angle(k,:),x_angle] = hist( Angle_boot(:), 200 );  % for bootstrap, set 200 bins temporarilly
    bin_sum = 0;
    n = 0;
    while ( n < (r(k)/bin) )
          n = n+1;
          bin_sum = bin_sum + hist_perm(k, n);
          p{k} = (perm_num - bin_sum)/ perm_num;    % calculate p value for HTI
    end 
    bin_sum = 0;
    n = 0;
    while ( bin_sum < 0.025*sum( hist_boot(k,:)) )   % define confidential value to be 0.05, now consider one side only which is 0.025 of confidence
          n = n+1;
          bin_sum = bin_sum + hist_boot(k, n);      
          HTI_boot(k) = r(k) - n * bin ;    % calculate what HTI value is thought to be significant different
    end 
%     bin_sum = 0;
%     n = 0;
%     while ( bin_sum < 0.975*sum( hist_boot_angle(k,:)) )   % define confidential value to be 0.05, now consider one side only which is 0.025 of confidence
%           n = n+1;
%           bin_sum = bin_sum + hist_boot_angle(k, n);      
%           Angle_boot(k) = x_angle(n) ;    
%     end 
end
%%------------------------------------------------------------------

%ANOVA modified by Aihua, it does not require whole trials, it does not matter if trial stopped during repetition
trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-2*(length(unique_azimuth)-1)) * length(unique_condition_num) + 1;
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


resp_mat_anova_hor = [];
for k=1: length(unique_condition_num)
   clear select_rep;
   for q=1:1:repetitions
       n = 0;
       for i=1:length(unique_azimuth) 
           select_rep_hor{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==0 & condition_num_rep{q}==unique_condition_num(k) );
           if (sum(select_rep_hor{q}) > 0)
               n = n+1;
               resp_mat_anova_hor{k}(q,n) = spike_rates_rep{q}(select_rep_hor{q})';
           end       
       end
   end
   [p_anova_hor, table_hor, stats_hor] = anova1(resp_mat_anova_hor{k},[],'off');
   P_anova_hor(k) = p_anova_hor;
   anova_table_hor{k} = table_hor;
   F_val_hor(k) = anova_table_hor{k}(2,5);
end
F_val_hor = cell2mat(F_val_hor);

%%------------------------------------------------------------------
% DDI Direction discrimination index   by Katsu 05/18/06
%--------------------------------------------------------------------
%each_stim_trials=repetitions*(length(unique_azimuth_s)*length(unique_elevation)-14) % 26 trajectory  -90+90 -45 0 + 45 *0 45 90 .......
each_stim_trials=repetitions*(length(unique_azimuth)*length(unique_elevation)-2*(length(unique_azimuth)-1)); % 26 trajectory  -90+90 -45 0 + 45 *0 45 90 .......
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
text(10,length(unique_condition_num),'Protocol       Spon   Minimum  Maximum   Azi     Ele      Amp    Std       HTI        HTIerr        p         F-val       p-ANOVA         DDI         p-ANOVA-Hor');
for k=1:length(unique_condition_num) 
    h_text{k}=num2str( [spon_resp, Min_resp(k), Max_resp(k), Vec_sum{k}, resp_std(k), r(k), HTI_boot(k), p{k}, F_val(k), P_anova(k), DDI(k),P_anova_hor(k)], 4);
    text(0,length(unique_condition_num)-k,h_title{k});
    text(10,length(unique_condition_num)-k,'Translation');
    text(20,length(unique_condition_num)-k, h_text{k} );
end

axis off;

% print;
% close;

