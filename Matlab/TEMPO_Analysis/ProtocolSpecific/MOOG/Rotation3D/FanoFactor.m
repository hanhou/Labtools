function FanoFactor(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
% h_title{1}='Vestibular';
% h_title{2}='Visual';
% h_title{3}='Combined';
unique_condition_num = munique(condition_num');

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));
% added by Katsu 111606
% spon_std = std(temp_spike_rates(spon_found));




% -------------------------------------------------------------------------
% %ANOVA modified by Aihua, it does not require whole trials, it does not matter if trial stopped during repetition
% trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-14) * length(unique_condition_num) + 1;
% repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);
% 
% % first parse raw data into repetitions, including null trials
% for q = 1:repetitions
%    azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%    elevation_rep{q} = temp_elevation(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%    condition_num_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%    spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% resp_mat_anova = [];
% for k=1: length(unique_condition_num)
%    clear select_rep;
%    for q=1:1:repetitions
%        n = 0;
%        for i=1:length(unique_azimuth)
%            for j=1:length(unique_elevation)
%                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==unique_elevation(j) & condition_num_rep{q}==unique_condition_num(k) );
%                if (sum(select_rep{q}) > 0)
%                    n = n+1;
%                    resp_mat_anova{k}(q,n) = spike_rates_rep{q}(select_rep{q})';
%                end
%            end
%        end
%    end
%    [p_anova, table, stats] = anova1(resp_mat_anova{k},[],'off');
%    P_anova(k) = p_anova;
%    anova_table{k} = table;
%    F_val(k) = anova_table{k}(2,5);
% end
% F_val = cell2mat(F_val);




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
                resp_mat_var(k, j, i) = var(spike_rates(select));   
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
                resp_mat_vector(k,j,1) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
                resp_mat_var(k, j, i) = 0;
                resp_mat_ste(k, j, i) = 0;
            end
        end        
    end
end

% % creat a real 3-D based plot where the center correspond to forward and
% % both lateral edges correspond to backward
%
% resp_mat_tran(:,:,1) = resp_mat(:,:,7);
% resp_mat_tran(:,:,2) = resp_mat(:,:,6);
% resp_mat_tran(:,:,3) = resp_mat(:,:,5);
% resp_mat_tran(:,:,4) = resp_mat(:,:,4);
% resp_mat_tran(:,:,5) = resp_mat(:,:,3);
% resp_mat_tran(:,:,6) = resp_mat(:,:,2);
% resp_mat_tran(:,:,7) = resp_mat(:,:,1);
% resp_mat_tran(:,:,8) = resp_mat(:,:,8);
% resp_mat_tran(:,:,9) = resp_mat_tran(:,:,1);


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate maximum and minimum firing rate
% max_res = max(max(max(resp_mat)))
% % max_vi_ve = max(max(max(resp_mat(2:3,:,:))));
% min_res = min(min(min(resp_mat_tran)))

% vector_num = length(unique_azimuth) * (length(unique_elevation)-2) + 2
% repeat = floor( length(spike_rates) / vector_num )

%%%%%%%%%%%%%%%%%%%%MATRIX%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:8
    mean_resp2(m)=resp_mat(1,2,m);% k(condition)=1, only for vestibuar
    var_resp2(m)=resp_mat_var(1,2,m);%:,2,m 2=-45 degrees for ele.
end
for m=1:8
    mean_resp3(m)=resp_mat(1,3,m);% ele = 0 deg.
    var_resp3(m)=resp_mat_var(1,3,m);
end
for m=1:8
    mean_resp4(m)=resp_mat(1,4,m);% ele. = 45  deg.
    var_resp4(m)=resp_mat_var(1,4,m);
end

    mean_resp1=resp_mat(1,1,1);% 1,1,1 is ele. -90 deg
    var_resp1=resp_mat_var(1,1,1);
    
    mean_resp5=resp_mat(1,5,1);% 1,5,1 = ele= 90 deg.
    var_resp5=resp_mat_var(1,5,1);
    
    mean_resp=[mean_resp1 mean_resp2 mean_resp3 mean_resp4 mean_resp5];
    var_resp=[var_resp1 var_resp2 var_resp3 var_resp4 var_resp5];
  
   
    
    for n=1:26 %%%%% To avoid log(0) turned to infinity...and To avoid varience is 0 x1*x2*..0=0
        if mean_resp(n)==0
            mean_resp(n)=0.001
        end
        if var_resp(n)==0
            var_resp(n)=0.001
        end
    end
    
    
    Min_resp=min(mean_resp);
    Max_resp=max(mean_resp);
    Max_var=max(var_resp);
    Min_var=min(var_resp);
    
    if Max_resp>Max_var
        Max_axis=Max_resp;
    else
        Max_axis=Max_var;
    end
    
% %%%%%%%%%%%%%%%%% calculate LOG Fano%%%%%%%%%%%%%%%%
% a=polyfit(log10(mean_resp),log10(var_resp),1)
% Fano=10^a(2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% calculate Geometric mean %%%%%%%%%%%%%%%%

% EachFano=var_resp./mean_resp;
% numbers=length(mean_resp);
% ff=[1];
% for i=1:26
% ff=ff*EachFano(i);
% end
% 
% GeoFano=ff^(1/26)

% There was MATLAB Command
EachFano=var_resp./mean_resp;
GeoFano=geomean(EachFano);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%% plot %%%%%%%%%%
figure(5)
TITLEgeo=[FILE  '   ;  Fano Factor = ' num2str(GeoFano)];
plot(mean_resp,var_resp,'o');xlim([0 Max_axis]);ylim([0 Max_axis]);xlabel('mean');ylabel('variance');title(TITLEgeo);axis square;


% figure(6)
% TITLE=[FILE  '   ;  Fano Factor = ' num2str(Fano)];
% plot(log10(mean_resp),log10(var_resp),'o');xlim([0 3]);ylim([0 3]);xlabel('log10 scale mean');ylabel('log10 scale variance');title(TITLE);axis square;


%--------------------------------------------------------------------------
%  Output file for Geometric Mean Fano
%--------------------------------------------------------------------------
buff = sprintf('%s\t %4.2f\t   %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
     FILE, spon_resp, Min_resp, Max_resp, Min_var, Max_var, GeoFano);
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_GeoFano.dat'];
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Darkness_GeoFano.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
%     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t');
fprintf(fid, 'FILE\t SPon\t minFR\t maxFR\t minV\t maxV\t Fanofactor\t ');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
% %------------------------

%---------------------------------------------------------------------------------------
% 
% %--------------------------------------------------------------------------
% %  Output file
% %--------------------------------------------------------------------------
% buff = sprintf('%s\t %4.2f\t   %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
%      FILE, spon_resp, Min_resp, Max_resp, a(2), Fano);
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Rotation_Fano.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\Darkness_Fano.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
% %     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t');
% fprintf(fid, 'FILE\t SPon\t min\t max\t polyfit\t Fanofactor\t ');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% %---------------------------------------------------------------------------------------
% %------------------------
return;