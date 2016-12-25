% DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	YONG, 12/10/08  
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_3D_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
plfp_chan = 2;

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_duration = data.moog_params(DURATION,:,MOOG);
temp_spike_data = data.spike_data(SpikeChan, :);

if abs(sum(sum(sum(data.plfp_data)))) >0
   temp_plfp_data(1,1:5000,1:length(temp_azimuth)) = data.plfp_data(plfp_chan,1:2:10000,:);
else
   temp_plfp_data(1,1:5000,1:length(temp_azimuth)) = 1;
end
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
duration = temp_duration(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);
plfp_data(:,:) = squeeze(temp_plfp_data(1,:,~null_trials & select_trials));
plfp_data_null(:,:) = squeeze(temp_plfp_data(1,:,null_trials));
plfp_null(1,:) = median(plfp_data_null(:,:),2);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_duration = munique(duration');

condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition_num = munique(condition_num');
StartEventBin(1)=996;
% calculate LFP value for each trial, use RMS
for i=1:length(azimuth)
    offset_temp = mean( plfp_data(StartEventBin(1)+115:StartEventBin(1)+115+unique_duration, i) );
    plfp_temp = plfp_data(StartEventBin(1)+115:StartEventBin(1)+115+unique_duration, i)-offset_temp;
    plfp_rates(i) = sqrt(sum(plfp_temp.^2));
    % use fft    
    yy = fft(plfp_temp,512); % 1k Hz
    Pyy = yy.* conj(yy) / 512; % only use the first 155 point, which ends at 300.78 Hz
    Pyy = Pyy./max(Pyy); % normalized
    ff = 1000*(0:256)/512;
    for j = 1:15
        lowend = 0+(j-1)*20;
        highend = 20+(j-1)*20;
        fre_in = find(ff>lowend & ff < highend);
        plfp_power(j,i) = sum(Pyy(fre_in));
    end    
    % the following is for the time course and frequecy components with
    % running mean
    startbin = 200;
    stopbin = unique_duration;
    for timecourse = 1:(stopbin-startbin)/10+1 % 200ms interval sliding with 10ms
        offset_temp = mean( plfp_data(StartEventBin(1)+115+(timecourse-1)*10:StartEventBin(1)+115+(timecourse-1)*10+startbin, i) );
        plfp_temp = plfp_data(StartEventBin(1)+115+(timecourse-1)*10:StartEventBin(1)+115+(timecourse-1)*10+startbin, i)-offset_temp;
        yy = fft(plfp_temp,512); % 1k Hz
        Pyy = yy.* conj(yy) / 512; % only use the first 155 point, which ends at 300.78 Hz
        Pyy = Pyy./max(Pyy); % normalize
        for freqency = 1: 29 % 20 interval sliding with 10
            lowend = 0 + (frequency-1)*10;
            highend = 20 + (frequency-1)*10;
            fre_in = find(ff>lowend & ff < highend);
            plfp_power{i}(frequency, timecourse) = sum(Pyy(fre_in));
        end
    end        
end

% % for spontaneous
% offset_temp = mean( plfp_null(1,StartEventBin(1)+115:StartEventBin(1)+115+unique_duration) );
% plfp_null = sum((plfp_null(1,StartEventBin(1)+115:StartEventBin(1)+115+unique_duration)-offset_temp).^2);
% plfp_null = sqrt(plfp_null);

% plot option, if regular plot, set to 0, if lambert plot, set to 1
% lamber_plot = 0;  % regular plot with elevation in a linear step
lamber_plot = 1; % lambert plot with elevation in a sin-transformed way
repetition = floor( length(azimuth) / (26*length(unique_stim_type)) );
%%%% Usually, axis azimuth from left is 270-225-180-135-90--0---90 %%%% 

unique_azimuth_s=[0 45 90 135 180 225 270 315];
temp=fliplr(unique_azimuth');
unique_azimuth_plot=[temp(2:end) temp(1:2)];clear temp
if lamber_plot == 1
    x_elevation=[-1,-0.707,0,0.707,1]; %uncomment for cosine plot
elseif lamber_plot ==0
    x_elevation=unique_elevation;%05/22/06 Katsu changed not cosin axis
end
% for cosine plot
%---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
azi_cos = [1,2,3,4,5,6,7,8,9];
ele_sin = [-1,-0.707,0,0.707,1];  

resp_mat = []; 
% for k=1: length(unique_condition_num) 
%     figure(k+1);
%     set(k+1,'Position', [5,15 980,650], 'Name', '3D Direction Tuning');
%     orient landscape;
%     %set(0, DefaultAxesXTickMode, 'manual', DefaultAxesYTickMode, 'manual', 'DefaultAxesZTickMode', 'manual');
%     axis off;
%     xoffset = 0;
%     yoffset = 0;
%     for f=1:15
%         spike_rates = plfp_power(f,:); % use each band           
% 
%         pc=0;
%         for j=1:length(unique_elevation)
%             for i=1:length(unique_azimuth)        
%                 select = find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
%                 if (sum(select) > 0)  
%                     pc = pc+1;
%                     resp_mat(k, j, i) = mean(spike_rates(select));       
% 
%                     spike_temp = spike_rates(select);
%                     resp_trial_anova1{k}(1:repetition,pc) =  spike_temp(1:repetition);    % for later anova use
%                     tempsum = 0;
%                     for tr = 1:length(select)
%                         tempsum = tempsum + plfp_power{select(tr)}(:,:));                        
%                     end
%                     resp_timecoursefrequency{i,j,k}(:,:) = tempsum./length(select);
%                 else
%     %                resp_mat_trial{k}(t, j, i) = 0;
%                     resp_mat(k, j, i) = resp_mat(k,j,1);
%                 end
%             end        
%         end    
%         P_anova(f,k) =  anova1( resp_trial_anova1{k}(:,:),'','off' ); 
% 
%         for i=1:length(unique_azimuth_plot)
%             Index(i)=find(unique_azimuth==unique_azimuth_plot(i));
%             resp_mat_tran(k,:,i) = resp_mat(k,:,Index(i));
%         end     
% 
%         axes('position',[0.05+xoffset 0.75-yoffset 0.20 0.20]);
% 
%         if lamber_plot ==1
%              contourf( azi_cos, ele_sin, squeeze( resp_mat_tran(k,:,:)) );
%         elseif lamber_plot ==0
%              contourf( squeeze( resp_mat_tran(k,:,:)) ); 
%         end
%         colorbar;
%      %   caxis([min(min(spike_rates)), max(max(spike_rates))]);
%         % make 0 correspond to rightward and 180 correspond to leftward
%         set(gca, 'ydir' , 'reverse');  
%         title( num2str(P_anova(f,k)) );
%         set(gca, 'xtick',x_elevation);
%         if lamber_plot == 1
%            ylim([-1, 1]);
%            set(gca, 'YTickMode','manual');
%            set(gca, 'ytick',[-1,-0.707,0,0.707,1]);
%            set(gca, 'yticklabel','-90|-45|0|45|90'); 
%         elseif lamber_plot ==0
%            ylim([-90, 90]);
%         end
%         set(gca, 'XTickMode','manual');
%         set(gca, 'xtick',[]);
%    %     set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
%    %     set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); % Katsu
%  
%         xoffset = xoffset+0.25;
%         if xoffset>0.9
%             xoffset = 0;
%             yoffset = yoffset+0.25;
%         end
% 
%     end
% end

for k=1: length(unique_condition_num) 
    figure(k+1);
    set(k+1,'Position', [5,15 980,650], 'Name', '3D Direction Tuning');
    orient landscape;
    %set(0, DefaultAxesXTickMode, 'manual', DefaultAxesYTickMode, 'manual', 'DefaultAxesZTickMode', 'manual');
    axis off;       

    pc=0;
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)        
            select = find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
            if (sum(select) > 0)  
                for tr = 1:length(select)
                    tempsum = tempsum + plfp_power{select(tr)}(:,:));                        
                end
                resp_timecoursefrequency{i,j,k}(:,:) = tempsum./length(select);
                % plot figure
                axes('position',[0.05+(i-1)*0.1 0.95-(j-1)*0.2 0.1 0.2]);
                contourf(resp_timecoursefrequency{i,j,k}(:,:));
            end
        end        
    end
end


% %---------------------------------------------------------------------------------------
% %Also, write out some summary data to a cumulative summary file
% buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
%      FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, p(:), DDI(:), P_anova );
% 
% outfile = ['Z:\Users\Yong\FEF\3D_Direction_Tuning_batch.dat']; 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
% %     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t Veb_DDI\t Vis_DDI\t Com_DDI\t Veb_var_term\t Vis_var_term\t Com_var_term\t');
%     fprintf(fid, 'FILE\t SPon\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% 
% return;
