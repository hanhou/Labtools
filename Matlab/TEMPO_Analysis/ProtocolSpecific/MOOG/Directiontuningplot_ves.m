%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	GCD, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------
function Directiontuningplot_ves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_trials ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
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


% decide whether protocol is sigma fixed or not
if (length(unique_stim_type) > 1 )
    condition_num = stim_type;
    h_title{1}='Vestibular';
    h_title{2}='Visual';
    h_title{3}='Combined';
else
    condition_num = amplitude;
    h_title{1}='amp=0.03, sig=6';
    h_title{2}='amp=0.08, sig=6';
    h_title{3}='amp=0.13, sig=6';
end
unique_condition_num = munique(condition_num');

%% ADD CODE HERE FOR PLOTTING
resp_mat = [];
for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=1: length(unique_condition_num)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(spike_rates(select));
                resp_temp(k, j, i) = mean(spike_rates(select));
                resp_mat_std(k, j, i) = std(spike_rates(select));    % calculate std between trials for later DSI usage                
            else
                resp_mat(k, j, i) = 0;
                resp_temp(k, j, i) = resp_temp(k,j,1);
                resp_mat_std(k, j, i) = 0;
            end
        end        
    end
end

%% creat a real 3-D based plot where the center correspond to forward and
%% both lateral edges correspond to backward
resp_mat_tran = [];
for i=1:( length(unique_azimuth)+1 )        % add a azimuth 360 to make it circularlly continuously
    for j=1:length(unique_elevation)
        for k=1:length(unique_condition_num)
            if ( j == 1 | j==5 )                          % fill NaN point with data
                resp_mat_tran(k, j, i) = resp_mat(k, j, 1);   
            else
                if (i < 8 )
                    resp_mat_tran(k, j, i) = resp_mat(k,j, 8-i);
                elseif(i==8)
                    resp_mat_tran(k, j, i) = resp_mat(k,j, 8);
                else
                    resp_mat_tran(k, j, i) = resp_mat(k,j, 7);
                end
            end
        end        
    end
end
% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));

% calculate maximum and minimum firing rate
max_res = max(max(max(resp_mat)));
min_res = min(min(min(resp_mat_tran)));

% deal with only one vector having response, numvector should be 6
% judge method:
% for k=1: length(unique_condition_num)
%     [row_max, col_max] = find( squeeze(resp_mat(k,:,:)) == max(max(squeeze(resp_mat(k,:,:)))) );
%     test_val = [ resp_mat(k,row_max-1,col_max), resp_mat(k,row_max+1,col_max),resp_mat(k,row_max,col_max-1),resp_mat(k,row_max,col_max+1)];
%     if test_val / max(max(squeeze(resp_mat(k,:,:)))) < 
% end
% row_max
% col_max
% number of vectors in each condition 
vector_num = length(unique_azimuth) * (length(unique_elevation)-2) + 2;
%vector_num = 6 ;

%------------------------------------------------------------------
% Define figure

xoffset=0;
yoffset=0;
figure(2);
set(2,'Position', [5,15 1200,900], 'Name', '3D Direction Tuning');
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
    caxis([min_res, max_res]);    
    colorbar;
    % make 0 correspond to rightward and 180 correspond to leftward
    set(gca, 'ydir' , 'reverse');
    set(gca, 'xtick', [] );
    set(gca, 'ytick', [] );
    title( h_title{k} );

    % plot 1-D for mean respond as a function of elevation
    axes('position',[0.06+xoffset 0.54+yoffset 0.04 0.24]);
    for j=1:length(unique_elevation)
        y_elevation_mean(1,j)=mean(resp_mat_tran(k,j,:));
        y_elevation_std(1,j) =std( spike_rates([find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )]) );
        y_elevation_ste(1,j) =y_elevation_std(1,j)/ sqrt(length(find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
    end
    x_elevation=-90:45:90;
    errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'o-');
    xlabel('Elevation');
    view(90,90);
    set(gca, 'xtick',[-90,-45,0,45,90]);
    xlim([-90, 90]);
    ylim([min(y_elevation_mean(1,:))-max(y_elevation_ste(1,:)), max(y_elevation_mean(1,:))+max(y_elevation_ste(1,:))]);

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
    errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'o-');
    xlim( [1, length(unique_azimuth)+1] );
    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
    xlabel('Azimuth');
    ylim([min(y_azimuth_mean(1,:))-max(y_azimuth_ste(1,:)), max(y_azimuth_mean(1,:))+max(y_azimuth_ste(1,:))]);
    
    xoffset=xoffset+0.48;
    
    % calculate min and max firing rate, standard deviation, DSI, Vectorsum
    Min_resp(k) = min( min( resp_mat_tran(k,:,:)) );
    Max_resp(k) = max( max( resp_mat_tran(k,:,:)) );
    resp_std(k) = sum( sum(resp_mat_std(k,:,:)) ) / vector_num;  % notice that do not use mean here, its 26 vectors intead of 40
    M=squeeze(resp_mat(k,:,:));     % notice that here DSI should use resp_temp without 0 value set manually
    DSI_temp(k) = DSI(M,spon_resp,resp_std(k));
    N=squeeze(resp_mat(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually 
    [Azi, Ele, Amp] =vectorsum(N);
    Vec_sum{k}=[Azi, Ele, Amp];
    
end

%-------------------------------------------------------------------
%check significance of DSI and calculate p value
perm_num=1000;
bin = 0.005;
for n=1: perm_num 
    for k=1:length(unique_condition_num)   
        spike_rates_pe{k} = spike_rates( find( condition_num==unique_condition_num(k) ) );
        spike_rates_pe{k} = spike_rates_pe{k}( randperm(length(spike_rates_pe{k})) );
    end
    % get the permuted spikerate to re-calculate DSI for each condition
    spike_rates_perm=[spike_rates_pe{1},spike_rates_pe{2},spike_rates_pe{3}];
    % re-creat a matrix similar as resp_mat              
    resp_vector_perm = [];
    for i=1:length(unique_azimuth)
        for j=1:length(unique_elevation)
            for k=1:length(unique_condition_num)
                select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
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
    % re-calculate DSI now
    for k=1: length(unique_condition_num)
        resp_perm_std(k) = sum( sum(resp_mat_perm_std(k,:,:)) ) / vector_num; 
        M_perm=squeeze(resp_mat_perm(k,:,:));
        DSI_perm(k,n) = DSI(M_perm, spon_resp, resp_perm_std(k) );
    end
    
end
x_bin = 0 : bin : 1;
for k = 1 : length(unique_condition_num)
    histo(k,:) = hist( DSI_perm(k,:), x_bin );
    bin_sum = 0;
    n = 0;
    while ( n < (DSI_temp(k)/bin) )
          n = n+1;
          bin_sum = bin_sum + histo(k, n);
          p{k} = (perm_num - bin_sum)/ perm_num;    % calculate p value
    end 
end

%------------------------------------------------------------------

% Now show vectorsum, DSI, p and spontaneous at the top of figure
axes('position',[0.05,0.85, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_condition_num)] );
h_spon = num2str(spon_resp);
text(0, length(unique_condition_num), FILE);
text(15,length(unique_condition_num),'Spon            Minimum        Maximum       Azi             Ele                Amp           Std             DSI                   p');
for k=1:length(unique_condition_num) 
    h_text{k}=num2str( [spon_resp, Min_resp(k), Max_resp(k), Vec_sum{k}, resp_std(k), DSI_temp(k), p{k} ] );
    text(0,length(unique_condition_num)-k,h_title{k});
    text(15,length(unique_condition_num)-k, h_text{k} );
end

axis off;

%---------------------------------------------------------------------------------------
%ALso, write out some summary data to a cumulative summary file

buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
     FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, DSI_temp, p{:} , resp_std );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\DirectionTuning_accel.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_DSI\t Vis_DSI\t Comb_DSI\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);

%---------------------------------------------------------------------------------------

return;