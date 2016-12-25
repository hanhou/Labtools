function Rotation3D_DDIperm(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);

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
spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');


condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition_num = munique(condition_num');

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));
% added by Katsu 111606
spon_std = std(temp_spike_rates(spon_found))




% -------------------------------------------------------------------------
%ANOVA modified by Aihua, it does not require whole trials, it does not matter if trial stopped during repetition
trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-14) * length(unique_condition_num) + 1;
repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);

% first parse raw data into repetitions, including null trials
for q = 1:repetitions
   azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   elevation_rep{q} = temp_elevation(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   condition_num_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Comment 090506 %%%%%%%%%%%Uncomment 090606

%%%%%%%%%Comment 090706 %%%%%%%%%%%Uncomment 090806

%%%%%%%%%Comment 111406 %%%%%%%%%%%Uncomment 111406
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
                resp_mat_vector(k, j, i) = mean(spike_rates(select));
                for t = 1 : length(spike_rates(select));              % this is to calculate response matrix based on each trial
                    spike_temp = spike_rates(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(spike_rates(select));     % calculate std between trials for later DSI usage
                resp_mat_ste(k, j, i) = resp_mat_std(k, j, i)/ sqrt(length(find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
            else
%                resp_mat_trial{k}(t, j, i) = 0; % From the begining
%                (071306) this sentence was commented by Katsu
%Following was wrong, I think from 07/10/06 by Katsu
%                 resp_mat(k, j, i) = 0;
%                 resp_mat_vector(k, j, i) = resp_mat_vector(k, j, 1);
%                 resp_mat_std(k, j, i) = 0;
%                 resp_mat_ste(k, j, i) = 0;
%So corrected by Katsu 07/15/06
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,i) =0; % for vector sum only % once this value was 1, not i %%011707 by Katsu
                resp_mat_std(k, j, i) = 0;
                resp_mat_std(k, j, i) = 0;
                resp_mat_ste(k, j, i) = 0;
            end
        end        
    end
end

% % creat a real 3-D based plot where the center correspond to forward and
resp_mat_tran(:,:,1) = resp_mat(:,:,7);
resp_mat_tran(:,:,2) = resp_mat(:,:,6);
resp_mat_tran(:,:,3) = resp_mat(:,:,5);
resp_mat_tran(:,:,4) = resp_mat(:,:,4);
resp_mat_tran(:,:,5) = resp_mat(:,:,3);
resp_mat_tran(:,:,6) = resp_mat(:,:,2);
resp_mat_tran(:,:,7) = resp_mat(:,:,1);
resp_mat_tran(:,:,8) = resp_mat(:,:,8);
resp_mat_tran(:,:,9) = resp_mat_tran(:,:,1);


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
set(2,'Position', [5,15 980,650], 'Name', 'Rotation_3D');
orient landscape;
%set(0, DefaultAxesXTickMode, 'manual', DefaultAxesYTickMode, 'manual', 'DefaultAxesZTickMode', 'manual');
axis off;

for k=1: length(unique_condition_num) 
    
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.4;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
    contourf( squeeze( resp_mat_tran(k,:,:)) );
    % set the same scale for visual and combined conditions but here assuming vestibular response is always smaller than that in visual and
    % combined conditions
%     if ( k==2 | k==3 )
%    caxis([min_res, max_res]);
%    caxis([0, 60]);
%     end

    colorbar;
    % make 0 correspond to rightward and 180 correspond to leftward
    set(gca, 'ydir' , 'reverse');
    set(gca, 'xtick', [] );
    set(gca, 'ytick', [] );    
    title( h_title{k} );
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % plot 1-D for mean respond as a function of elevation
    axes('position',[0.06+xoffset 0.54+yoffset 0.04 0.24]);
    for j=1:length(unique_elevation)
        y_elevation_mean(1,j)=mean(resp_mat_tran(k,j,:));
        y_elevation_std(1,j) =std( spike_rates([find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )]) );
        y_elevation_ste(1,j) =y_elevation_std(1,j)/ sqrt(length(find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
    end
    x_elevation=unique_elevation;
    errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'ko-');%----------Now disable temporaly
%     errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'k-');% Katsu for paper settle for m3c294
    xlabel('Elevation');
    view(90,90);
    set(gca, 'xtick',unique_elevation);
    xlim([-90, 90]);
    ylim([min(y_elevation_mean(1,:))-max(y_elevation_ste(1,:)), max(y_elevation_mean(1,:))+max(y_elevation_ste(1,:))]);%----------Now disable temporaly

    % plot 1-D for mean respond as a function of azimuth
    axes('position',[0.11+xoffset 0.46+yoffset 0.274 0.06]);
    for i=1:(length(unique_azimuth) )
        y_azimuth_mean(1,i)=mean(resp_mat_tran(k,:,i));
        y_azimuth_std(1,i) =std( spike_rates([find( (azimuth==unique_azimuth(i))&(condition_num==unique_condition_num(k)) )]) );
        y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (azimuth==unique_azimuth(i))&(condition_num==unique_condition_num(k)) )) );    
    end
    y_azimuth_mean(1,9) = mean(resp_mat_tran(k,:,1));
    for i=1:( length(unique_azimuth)+1 )
        if (i < 8)        
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8-i);
        elseif (i == 8)
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8);
        else
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,7);
        end
    end
    x_azimuth=1:(length(unique_azimuth)+1);
    errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'ko-');%----------Now disable temporaly
%     errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'k-');% Katsu for paper settle for m3c294
    xlim( [1, length(unique_azimuth)+1] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); 
    xlabel('Azimuth');
    ylim([min(y_azimuth_mean(1,:))-max(y_azimuth_ste(1,:)), max(y_azimuth_mean(1,:))+max(y_azimuth_ste(1,:))]);%----------Now disable temporaly

    xoffset=xoffset+0.48;
    
    % calculate min and max firing rate, standard deviation, HTI, Vectorsum
    Min_resp(k) = min( min( resp_mat_tran(k,:,:)) );
    Max_resp(k) = max( max( resp_mat_tran(k,:,:)) );
    resp_std(k) = sum( sum(resp_mat_std(k,:,:)) ) / vector_num;  % notice that do not use mean here, its (length(unique_azimuth)*length(unique_elevation)-14) vectors intead of 40
    
    Ave_26{k}=sum(sum(squeeze(resp_mat_vector(k,:,:))))/26; %resp_mat_vec includes 0, so sum*sum devided by 26 trajectories
    
    M=squeeze(resp_mat_vector(k,:,:));     % notice that here DSI should use resp_temp without 0 value set manually
    % this part is to calculate vestibular gain
    resp_onedim{k} = [M(1,1),M(2,:),M(3,:),M(4,:),M(5,1)]';     % hard-code temperarilly    
    N=squeeze(resp_mat_vector(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually 
    [Azi, Ele, Amp] = vectorsum(N);
    Vec_sum{k}=[Azi, Ele, Amp];
    % Heading Tuning Index
    r(k) = HTI(M,spon_resp);   % call HTI function    
end
% calculate vestibular gain by a*ves=comb-vis
if (length(unique_stim_type) == 3)
   [bb,bint,rr,rint,stats] = regress( (resp_onedim{3}-resp_onedim{2}), [ones(vector_num,1),(resp_onedim{1}-spon_resp)] );    % with offset
   gain = bb(2);
else
   gain = NaN;
end
unique_stim_type
Max_resp
Min_resp
Ave_26;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % linear sum model Modified by Yong (for Katsu) 06/19/07
% % R(comb) = a1 x R(vest) + a2 x R(visu) + a3  ; This is the model
% % -------------------------------------------------------------------------
% % %%%%%%%%         Note; Run only 3 conditions neurons         %%%%%%%%%5
% % ----------------------------------------------------------------------------                  
% %
% A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
% OPTIONS = optimset('fmincon');
% OPTIONS = optimset('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
% % yy1 =
% % @(x)sum((resp_onedim{1}*x(1)+resp_onedim{2}*x(2)+x(3)-resp_onedim{3}).^2);  1st original
% yy1 = @(x)sum( ((resp_onedim{1}-spon_resp)*x(1)+(resp_onedim{2}-spon_resp)*x(2)+x(3)-(resp_onedim{3}-spon_resp)).^2 );  %w1*ve+w2*vi
% es1 = [0.5,0.5,0];
% LB1 = [0,0,-100];% once [0,0,0] once [-5, -5, 0]
% UB1 = [5,5,100];
% 
% v1 = fmincon(yy1,es1,A,B,Aeq,Beq,LB1,UB1, NONLCON, OPTIONS); % fminsearch   
% v1
% gain_ratio = v1(1) / v1(2)
% VVR = (Max_resp(2)-spon_resp) / (Max_resp(1)-spon_resp)
% VVR_new = (Max_resp(2)-Min_resp(2)) / (Max_resp(1)-Min_resp(1))
% % 
% % %%% To plot model fit is correcto or not
% % % pred = resp_onedim{1}*v1(1)+resp_onedim{2}*v1(2)+v1(3);
% % % figure;% to make sure fit (green) is good or not
% % % plot(resp_onedim{1},'k-');
% % % hold on;
% % % plot(resp_onedim{2},'r-');
% % % plot(resp_onedim{3},'g-');
% % % plot(pred,'b-');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%------------------------------------------------------------------
% DDI Direction discrimination index   by Katsu 05/18/06
%--------------------------------------------------------------------
%spike_rates_sqrt = sqrt(spike_rates);That is not original, it makes value bigger by 
each_stim_trials=repetitions*(length(unique_azimuth)*length(unique_elevation)-14) % (length(unique_azimuth)*length(unique_elevation)-14) trajectory  -90+90 -45 0 + 45 *0 45 90 .......
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
                   SSE_term(k,j,i)=0;%This is correct extra -90 and 90 should be 0, to be summed
            end
         end
     end
%      SSE_azimth_sum(k,j)=sum(SSE_term(k,j,:))  %it is not correct
%      SSE_total(k)=sum(SSE_azimth_sum(k,:))
     SSE_total(k)=sum(sum(SSE_term(k,:,:)));
     max_min_term(k)=(Max_resp(k)-Min_resp(k))/2;
     var_term(k)=sqrt(SSE_total(k)/(each_stim_trials-n));
     DDI(k)=max_min_term(k)/(max_min_term(k)+var_term(k));
end
r
DDI
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %check significance of DDI and calculate p value, Not !!! do bootstrap at the same time to test value varience
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
% 
%%%%%%   commended 1000 boot strap

perm_num=1000;
bin = 0.005;
spike_rates_perm = [];
SSE_term_perm=[];
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
    for k=1:length(unique_condition_num)
        m=0;
        
        for j=1:length(unique_elevation)
            for i=1:length(unique_azimuth)
                
                select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
                if (sum(select) > 0)
                    m=m+1;
                    
                    resp_mat_perm(k,j,i) = mean(spike_rates_perm(select));
%                     resp_mat_perm_std(k,j,i) = std(spike_rates_perm(select));
                    SSE_term_perm(k,j,i)=sum((spike_rates_perm(select)-mean(spike_rates_perm(select))).^2);
                else
                    resp_mat_perm(k,j,i) = resp_mat_perm(k,j,1);% Note! This is not =0. In order to calcurate Max or Min, line 349-350.
%                     resp_mat_perm_std(k,j,i) = 0;
                    SSE_term_perm(k,j,i) = 0;% This will be sum, so need to be 0
                end
            end        
        end
    end
    
    
    
    % re-calculate HTI and DDI now
    for k=1: length(unique_condition_num)
 %       resp_perm_std(k) = sum( sum(resp_mat_perm_std(k,:,:)) ) / vector_num; 
        M_perm=squeeze(resp_mat_perm(k,:,:));
        r_perm(k,n) = HTI(M_perm, spon_resp);
        
        Min_resp_perm(k,n) = min( min( resp_mat_perm(k,:,:) ));
        Max_resp_perm(k,n) = max( max( resp_mat_perm(k,:,:) ));
    
        SSE_total_perm(k,n)=sum(sum(SSE_term_perm(k,:,:)));
        max_min_term_perm(k,n)=(Max_resp_perm(k,n)-Min_resp_perm(k,n))/2;
        var_term_perm(k,n)=sqrt(SSE_total_perm(k,n)/(each_stim_trials-m));
        DDI_perm(k,n)=max_min_term_perm(k,n)/(max_min_term_perm(k,n)+var_term_perm(k,n));
        
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
    
    hist_DDI_perm(k,:) = hist( DDI_perm(k,:), x_bin );  % for DDI permutation
    
    
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
    while ( n < (DDI(k)/bin) )
          n = n+1;
          bin_sum = bin_sum + hist_DDI_perm(k, n);
          DDI_p{k} = (perm_num - bin_sum)/ perm_num;    % calculate p value for HTI
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



%%% There is  no HTI p and DDI p, so commented


 %---------------------------------------------------------------------------
% Now show vectorsum, DSI, p and spontaneous at the top of figure--------and DDI by Katsu 05/18/06
axes('position',[0.05,0.85, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_condition_num)] );
h_spon = num2str(spon_resp);
text(0, length(unique_condition_num), FILE);
text(10,length(unique_condition_num),'Protocol  Spon   Minimum   Maximum    Azi      Ele       Amp     Std        HTI          HTIerr      p-HTI       F-val        p-ANOVA       DDI   p-DDI');
for k=1:length(unique_condition_num) 
    h_text{k}=num2str( [spon_resp, Min_resp(k), Max_resp(k), Vec_sum{k}, resp_std(k), r(k), HTI_boot(k), p{k}, F_val(k), P_anova(k), DDI(k), DDI_p{k}], 4);
    text(0,length(unique_condition_num)-k,h_title{k});
    text(9,length(unique_condition_num)-k,'Rotation');
    text(17,length(unique_condition_num)-k, h_text{k} );
end

axis off;



% %--------------------------------------------------------------------------
% %  Katsu 
% %--------------------------------------------------------------------------

% %% VVR gain analysis  version 062007
% sprint_txt = ['%s'];
% for i = 1 : 1000
%      sprint_txt = [sprint_txt, ' %4.4f'];    
% end
% % buff = sprintf('%s\t %4.3f\t  %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %4.3f\t  %4.3f\t %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t %4.3f\t   %4.3f\t %4.3f\t   %4.3f\t %4.3f\t   %4.3f\t', ...
% %      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot,  resp_std, gain, F_val, DDI, var_term, p{:} , P_anova, DDI_p{:});
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_DDIperm.dat'];
% %
% %  VVR nerw gain analysis for Katsu 062007
% buff = sprintf(sprint_txt,FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, gain, F_val, DDI,  P_anova,  v1, gain_ratio, VVR, VVR_new);
% %%   Original with 26 ave.
% % buff = sprintf(sprint_txt,FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot,  resp_std, gain, F_val, DDI, var_term, p{:} , P_anova, DDI_p{:}, Ave_26{:});
% %
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Translation_p_DDI.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_MST_Lothar.dat']
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_MT_p_DDI.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Laby_26ave.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_3combined0604_35%.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_OVER100.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_gain_VVR_new_0625.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
% %     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t');
% fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_min\t c_max\t Vis_max\t c_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t c_azi\t c_ele\t c_amp\t gain\t Veb_F_anova\t Vis_F_anova\t c_F_anova\t Veb_DDI\t Vis_DDI\t c_DDI\t Veb_P_anova\t Vis_P_anova\t c_P_anova\t a1\t a2\t a3\t gain_ratio\t VVR\t VVR_new\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);


%%%% original including HTI p DDI p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
sprint_txt = ['%s'];
for i = 1 : 1000
     sprint_txt = [sprint_txt, ' %4.4f'];    
end
% buff = sprintf('%s\t %4.3f\t  %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %4.3f\t  %4.3f\t %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t  %4.3f\t   %4.3f\t %4.3f\t   %4.3f\t %4.3f\t   %4.3f\t %4.3f\t   %4.3f\t', ...
%      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot,  resp_std, gain, F_val, DDI, var_term, p{:} , P_anova, DDI_p{:});
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_DDIperm.dat'];
%
%  VVR nerw gain analysis for Katsu
% buff = sprintf(sprint_txt,FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot,  resp_std, gain, F_val, DDI, var_term, p{:} , P_anova, DDI_p{:}, Ave_26{:}, v1, gain_ratio, VVR);
%%%   Original with 26 ave.
buff = sprintf(sprint_txt,FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot,  resp_std, gain, F_val, DDI, var_term, p{:} , P_anova, DDI_p{:}, Ave_26{:});
%
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Translation_p_DDI.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_MST_Lothar.dat']
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_MT_Lase_all_2.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Laby_26ave.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_3combined0604_35%.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_OVER100.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_gain_VVR.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
%     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t');
fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_DDI\t Vis_DDI\t Veb_var_term\t Vis_var_term\t Veb_P\t Vis_P\t Veb_P_anova\t Vis_P_anova\t Veb_p_DDI\t Vis_p_DDI\t 26Ave\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
% %---------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------
return;