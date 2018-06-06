% DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	YONG, 12/10/08  

% edited by LBY 2016/06
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_3D_LBY(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);

%contains protocol specific keywords - 1/4/01 BJP
TEMPO_Defs;
Path_Defs;
ProtocolDefs; 


% 以下为可以进行改动的一些声明
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';

StartEventBin(1)=996;% an estimate value

%SpikeChan = 1;
plfp_chan = 1;
analyze_lfp = 0; % firing rate based on spikes，此参数有两种模式，当lfp=0即为与lfp无关，发放率与spike的个数有关
%analyze_lfp = 1; % firing rate based on LFP

% plot option, if regular plot, set to 0, if lambert plot, set to 1
lamber_plot = 1; % lambert plot with elevation in a sin-transformed way

%% load the data and choose the trials to analyse
% to select trials fall between BegTrials & EndTrials and not spontenous 

%get the column of values for azimuth, elevation, stim_type, amplitude,
%duration, spike data and spike rates for temporary
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);%有3种
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);%前进距离
temp_duration = data.moog_params(DURATION,:,MOOG);%给刺激的时间
temp_spike_data = data.spike_data(SpikeChan, :);%发放的数据,有的1000，有的2000，分别表示1s和2s
temp_spike_rates = data.spike_rates(SpikeChan, :); %发放率 ，排除坏点用的                                                                                                                   

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical((temp_azimuth == data.one_time_params(NULL_VALUE)));

%now, remove trials that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_tri = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_tri ~= NaN)%有坏点
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_tri) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

% select trials fall between BegTrials & EndTrials and not spontenous
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
duration = temp_duration(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

% to acquire the only parameters
unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
% unique_elevation = [0]; % for Fisher information CZX 20130619
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_duration = munique(duration');
% the following is just for convenient
% condition_num is equal to stim_type
condition_num = stim_type;
unique_condition_num = munique(condition_num');

% to select trials for spike data (select trials fall between BegTrials &
% EndTrials and not spontenous)
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 99;%给丢弃的数据打标签
end
spike_data(1,:) = temp_spike_data( 1, find(temp_spike_data(1,:)~=99) );%丢弃打了标签的数据
spike_data(1, find(spike_data>10) ) = 1; %丢弃spike_data大于10这种不符合的数据

%% to acquire spike data in trials we need and do simple analysis for plot

% to calculate the spike rates during the middle 1s 
for ss =  1 : length(spike_rates) % ss marks the index of trial
    if unique_duration == 1500 % Added by CZX to be compatible with the stimulus duration 1500ms
        % use the whole 1 second 
        % suppose the delay is 115ms, to get the middle 1s, we have to
        % start analyzing data 250ms after stimuli onset and stop analyzing
        % 250ms before stimuli offset. CZX 2013-06-05
        spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+68+250 + 5000*(ss-1) : StartEventBin(1)+1068 +250+ 5000*(ss-1)) ) ; 
    end
end

% calculate spontaneous firing rate计算自发发放率
spon_found = find(null_trials(select_trials)==1); %null_trials=1，表示自发发放
spon_resp = mean(temp_spike_rates(spon_found))
spon_resp_all = temp_spike_rates(spon_found)%所有的自发组合在一起
% added by Katsu 111606
spon_std = std(temp_spike_rates(spon_found))

% to get the repetition number
% 这段程序到底有什么用呢 后面明明有一样的程序，太冗余了
pc=0;
for k=1: length(unique_condition_num)
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)        
            select = find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
            if (sum(select) > 0) % means there are trials in this condition
                pc=pc+1;
                trialrepeat(pc) = length(select);
            end
        end
    end
end
repetition = min(trialrepeat) % to get the min repeat number

% to know which repetitions they are
resp_mat = [];
for k=1: length(unique_condition_num)
    pc=0;
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)        
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
            % 判断是否有这种情况下的trial
            select_hori = find( (azimuth==unique_azimuth(i)) & (elevation==0) & (condition_num==unique_condition_num(k)) );
            % to find trials under horizontal conditions
            if (sum(select) > 0)  
                pc = pc+1;
                if analyze_lfp==1
                    resp_mat(k, j, i) = median(spike_rates(select)); 
                    % avoid large values
                else
                    resp_mat(k, j, i) = mean(spike_rates(select)); 
                end
                % to calculate mean spike rates for every conditon 
                % 看来在LFP情况下用的是中位数，SU情况下是平均值
                resp_mat_vector(k, j, i) = mean(spike_rates(select)); 
                % for vectorsum to calculate the preferred direction
                resp_sse(k,pc) = sum( (spike_rates(select)-mean(spike_rates(select))).^2 );
                % 分别计算每个condition下（n个重复）的误差平方和
                resp_trialnum(k,pc)= length(spike_rates(select));
                % 列出每个condition下的重复数
                spike_temp = spike_rates(select); % 暂时储存用
                resp_trial_anova1{k}(1:repetition,pc) =  spike_temp(1:repetition); 
                % for later anova use & plot
                spike_temp_hori = spike_rates(select_hori);
                resp_trial_hori{k}(1:repetition,i) = spike_temp_hori(1:repetition);
                % for later anova use & plot ( horizontal)
                resp_mat_std(k, j, i) = std(spike_rates(select));     
                resp_mat_ste(k, j, i) = resp_mat_std(k, j, i)/ sqrt(length(find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
                % calculate std & ste
                % ste = std/sqrt(n)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                resp_mean(k,pc) = mean(spike_rates(select));
                resp_SD(k,pc) = std(spike_rates(select));
                resp_ste(k,pc) = resp_SD(k,pc)/ sqrt(repetition); 
                % CZX for fine tuning%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%20130620%%%%%%%%%%%%%%% 
                % 求firing rate的最大最小值
%                 resp_max(k,pc) = max(mean(spike_rates(select)));
%                 resp_min(k,pc) = min(mean(spike_rates(select)));     
                %%%%%%%%%%%%20130620%%%%%%%%%%%%%%%
               
                % z-score data  这里都是什么…
                z_dist = spike_rates(select);
                if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
                   z_dist = (z_dist - mean(z_dist))/std(z_dist);
                else
                   z_dist = 0; 
                end
                Z_Spikes(select) = z_dist; 
                resp_z{j,i} = z_dist;
            else
%               resp_mat_trial{k}(t, j, i) = 0;
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,i) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
                resp_mat_ste(k, j, i) = 0;
                resp_z{j,i} = resp_z{j,1};
            end
        end        
    end
    
    P_anova(k) =  anova1(resp_trial_anova1{k}(:,:),'','off' );  %检测均值的可信度
    P_anova_hori(k) =  anova1(resp_trial_hori{k}(:,:),'','off' ); 
    
    resp_sse_DDI(k)=sum(resp_sse(k,:))/(sum(resp_trialnum(k,:))-26); % CZX
    DDI_new(k)=(max(resp_mean(k,:))-min(resp_mean(k,:)))/(max(resp_mean(k,:))-min(resp_mean(k,:))+2 *sqrt(resp_sse_DDI(k)));
    % 计算DDI
    % 以上为在每种stim type下计算的
end

P_anova % show ANOVA for different stim types

% % for fine tuning 
% if length(unique_azimuth)  ==2
%     disp('Fine Tuning!!')
%     Fine_Tuning
%     return
% elseif pc < 26
%     Horizontal_Fine_Tuning
%     disp('It should be a horizontal fine tuning!')
%     return
% %      keyboard
% else
%     disp('Normal condition and Go on !')
% end

%% Plot

% to plot every trials including repetitions
m{1}='b.-'; m{2}='rx-'; m{3}='ro-'; m{4}='r+-'; m{5}='rd-'; m{6}='rs-'; 
m{7}='rv-'; m{8}='r*-'; m{9}='r.-';m{10}='gx-'; m{11}='go-'; m{12}='g+-'; 
m{13}='gd-'; m{14}='gs-'; m{15}='gv-'; m{16}='g*-'; m{17}='g.-';m{18}='kx-'; 
m{19}='ko-'; m{20}='k+-'; m{21}='kd-'; m{22}='ks-'; m{23}='kv-'; m{24}='k*-'; 
m{25}='k.-';m{26}='m.-'; 

% figure(2);
% clf;
% n=0;
% for j=1:26 % 26 can be replaced by pc
%    n = n+1;
%    plot(resp_trial_anova1{1}(:, j), m{n});
%    hold on;
% end
% title(FILE);

% to plot the contour figure (3D tuning)

% to transform the data for plot
% Usually, axis azimuth from left is 270-225-180-135-90--0---90 %
unique_azimuth_s=[0 45 90 135 180 225 270 315];
temp=fliplr(unique_azimuth');
unique_azimuth_plot=[temp(2:end) temp(1:2)];
% Usually, axis azimuth from left is 270-225-180-135-90--0---90 %
clear temp
for i=1:length(unique_azimuth_plot)
    Index(i)=find(unique_azimuth==unique_azimuth_plot(i));
    resp_mat_tran(:,:,i) = resp_mat(:,:,Index(i)); % resp_mat -> mean firing rate
end

% calculate maximum and minimum firing rate
max_res = max(max(max(resp_mat)));
min_res = min(min(min(resp_mat_tran)));


% xoffset=0;
% yoffset=0;
% set(2,'Position', [15,15 1150,910], 'Name', '3D Direction Tuning');
% orient landscape;
% axis off;
% for cosine plot
azi_cos = [1,2,3,4,5,6,7,8,9];
ele_sin = [-1,-0.707,0,0.707,1];

for k=1: length(unique_condition_num) 
%     if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
%         yoffset = yoffset-0.41;
%         xoffset = 0;
%     end
%     axes('position',[0.11+xoffset 0.58+yoffset 0.32 0.24]);
    axes;
    set(gca,'units','points','pos',[520*k-380,340,400,300]);
    if lamber_plot ==1
         contourf( azi_cos, ele_sin, squeeze( resp_mat_tran(k,:,:)),'linecolor','w','linestyle','none' );
    elseif lamber_plot ==0
         contourf( squeeze( resp_mat_tran(k,:,:)),'linecolor','w','linestyle','none' ); 
    end
    
%     box off;
    % make 0 correspond to rightward and 180 correspond to leftward
    set(gca, 'ydir' , 'reverse');
    set(gca, 'xtick', [] );
    set(gca, 'ytick', [] );    
    box off;
%     axis off;
    title( h_title{unique_stim_type(k)} ,'fontsize',30 );
    % plot 1-D for mean respond as a function of elevation
    % notice that elevation scale is transformed by consine

%     set(gca, 'xtick', [] );
%     set(gca, 'ytick', [] );  
    for j=1:length(unique_elevation)
        y_elevation_mean(1,j)=mean(resp_mat_tran(k,j,:));
        y_elevation_std(1,j) =std( spike_rates([find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )]) );
        y_elevation_ste(1,j) =y_elevation_std(1,j)/ sqrt(length(find( (elevation==unique_elevation(j))&(condition_num==unique_condition_num(k)) )) );
    end
    if lamber_plot == 1
        x_elevation=[-1,-0.707,0,0.707,1]; %uncomment for cosine plot
    elseif lamber_plot ==0
        x_elevation=unique_elevation;%05/22/06 Katsu changed not cosin axis
    end
    axes;
    set(gca,'units','points','pos',[520*k-480,340,100,300]);
    errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'k.-','markerSize',30,'linewidth',3);
%     box off;
    axis off;
    xlabel('Elevation');
    view(90,90);
    set(gca, 'xtick',x_elevation);
    if lamber_plot == 1
       xlim([-1, 1]);
       set(gca, 'XTickMode','manual');
       set(gca, 'xtick',[-1,-0.707,0,0.707,1]);
       set(gca, 'xticklabel','-90|-45|0|45|90'); 
    elseif lamber_plot ==0
       xlim([-90, 90]);
    end
    ylim([min(y_elevation_mean(1,:))-max(y_elevation_ste(1,:)), max(y_elevation_mean(1,:))+max(y_elevation_ste(1,:))]);%axis off %----------Now add axis off

    % plot 1-D for mean respond as a function of azimuth
%       axes('position',[0.11+xoffset 0.40+yoffset 0.32 0.1]);
    for i=1:(length(unique_azimuth_s) )
        y_azimuth_mean(1,i)=mean(resp_mat_tran(k,:,i));
        y_azimuth_std(1,i) =std( spike_rates([find( (azimuth==unique_azimuth_s(i))&(condition_num==unique_condition_num(k)) )]) );
        y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (azimuth==unique_azimuth_s(i))&(condition_num==unique_condition_num(k)) )) );    
    end
    y_azimuth_mean(1,9) = mean(resp_mat_tran(k,:,1)); % an extra value for circle
    for i=1:( length(unique_azimuth_s)+1 )
        if (i < 8)        
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8-i);
        elseif (i == 8)
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8);
        else
            y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,7);
        end
    end
    % above is because transformation of 270|225|180|135|90|45|0|-45|-90
    x_azimuth=1:(length(unique_azimuth_s)+1);
    axes;
    set(gca,'units','points','pos',[520*k-380,220,400,100]);
    errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'k.-','markerSize',30,'linewidth',3);
%     set(gca,'color','w');
%     errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'k-');% Katsu for paper settle for m3c294
    xlim( [1, length(unique_azimuth)+1] );
%     xlim( [0.9, length(unique_azimuth_s)+1.1] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); % Katsu
%     set(gca, 'xticklabel','0|45|90|135|180|225|270|315|360'); 
    xlabel('Azimuth');
    ylim([min(y_azimuth_mean(1,:))-max(y_azimuth_ste(1,:)), max(y_azimuth_mean(1,:))+max(y_azimuth_ste(1,:))]);%axis off %----------Now add axis off
    axis off;  
    % calculate min and max firing rate, standard deviation, HTI, Vectorsum
    Min_resp(k) = min( min( resp_mat_tran(k,:,:)) );
    Max_resp(k) = max( max( resp_mat_tran(k,:,:)) );
    resp_std(k) = sum(resp_sse(k,:))/(sum(resp_trialnum(k,:))-26);
    M=squeeze(resp_mat(k,:,:));     % notice that here DSI should use resp_temp without 0 value set manually
    % this part is to calculate vestibular gain
    resp_onedim{k} = [M(1,1),M(2,:),M(3,:),M(4,:),M(5,1)]';     % hard-code temperarilly    
    N=squeeze(resp_mat_vector(k,:,:));      % notice that here vectorsum should use resp_mat with 0 value set manually 
    [Azi, Ele, Amp] = vectorsum(N); % preferred direction and r
    Vec_sum{k}=[Azi, Ele, Amp]; % preferred direction in each condition 
    % Heading Tuning Index
    r(k) = HTI(M,spon_resp);   % call HTI function  
    DDI(k) = (Max_resp(k)-Min_resp(k))/(Max_resp(k)-Min_resp(k)+2*sqrt(resp_std(k)));
end



%% analyse

%check significance of HTI and calculate p value, do bootstrap at the same time to test value varience
perm_num=1000;
bin = 0.005;
spike_rates_perm = [];
for n=1: perm_num
    % this is permuted based on trials
    % according to conditions
    for k=1:length(unique_condition_num)   
        spike_rates_pe{k} = spike_rates( find( condition_num==unique_condition_num(k)));
        spike_rates_pe{k} = spike_rates_pe{k}( randperm(length(spike_rates_pe{k})) );
    end

    % put permuted data back to spike_rates
    spike_rates_perm(length(spike_rates))=0; % set all values = 0
    for k=1:length(unique_condition_num) 
        ii = find(stim_type == unique_stim_type(k));
        spike_rates_perm(ii) = spike_rates_pe{k};
    end
    % 将数据打乱之后放入 spike_rates_perm 
    
    % re-creat a matrix similar as resp_mat but use data in spike_rates_perm              
    resp_vector_perm = [];
    for i=1:length(unique_azimuth)
        for j=1:length(unique_elevation)
            for k=1:length(unique_condition_num)
                select = logical((azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
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
    
    % re-calculate HTI & DDI now
    for k=1: length(unique_condition_num)
 %       resp_perm_std(k) = sum( sum(resp_mat_perm_std(k,:,:))) / vector_num; 
        M_perm=squeeze(resp_mat_perm(k,:,:));
        r_perm(k,n) = HTI(M_perm, spon_resp); 
        Min_resp_perm = min( min( resp_mat_perm(k,:,:)) );
        Max_resp_perm = max( max( resp_mat_perm(k,:,:)) );
        resp_std_perm = sum( sum(resp_mat_perm_std(k,:,:)) ) / 26;  
        % notice that do not use mean here, its 26 vectors intead of 40
        DDI_perm(k,n) = (Max_resp_perm-Min_resp_perm)/(Max_resp_perm-Min_resp_perm+resp_std_perm);
    end

end
% now calculate p value or significant test
for k = 1 : length(unique_condition_num)
    p_HTI(k) = length(find(r_perm(k,:)>=r(k)) )/perm_num; 
    p_DDI(k) = length(find( DDI_perm(k,:)>=DDI(k)) )/perm_num;
end


%% 

% VVR=( Max_resp(2)-Min_resp(2))/(Max_resp(1)-Min_resp(1))% for illustration, CZX 2012-11-11
% Re1=[resp_mean(1,:)' resp_mean(2,:)' ones(26,1)];
% Re2=resp_mean(3,:)';
% [b,bint,r,rint,stats] = regress(Re2,Re1);
% b

% %% text on the figure and save them
% 
% % Now show vectorsum, DSI, p and spontaneous at the top of figure
% axes;
% set(gca,'units','points','pos',[20,20,500*iik+100,100]);
% xlim( [0,100] );
% ylim( [0,length(unique_condition_num)+1] );
% h_spon = num2str(spon_resp);
% text(0, length(unique_condition_num)+1, FILE);
% text(20, length(unique_condition_num)+1, 'SpikeChan=');
% text(40, length(unique_condition_num)+1, num2str(SpikeChan));
% text(10,length(unique_condition_num),'Protocol      Spon      Minimum      Maximum       Azi        Ele             DDI    HTI      p-ANOVA');
% for k=1:length(unique_condition_num) 
%     h_text{k}=num2str( [spon_resp, Min_resp(k), Max_resp(k), Vec_sum{k}, DDI(k), r(k),P_anova(k)] );
%     text(0,length(unique_condition_num)-k,h_title{unique_stim_type(k)});
%     text(10,length(unique_condition_num)-k,'Translation');
%     text(25,length(unique_condition_num)-k, h_text{k} );
% end
% 
% axis off;
% set(gcf,'color','white')
% set(findall(gcf,'FontSize',10),'FontSize',13);
% 
% % save figure 2 & 3
% FileName_Temp = num2str(FILE);
% FileName_Temp =  FileName_Temp(1:end-4);%%%To get rid of the '.htb', or the formate of pictures can not be manipulated. 
% str1 = [FileName_Temp,'_Ch' num2str(SpikeChan)];
% str2 = [str1, '_raw'];
% saveas(3,['Z:\LBY\Recording Data\Qiaoqiao\3D_T_Tuning\' str1], 'pdf')
% saveas(2,['Z:\LBY\Recording Data\Qiaoqiao\3D_T_Tuning\' str2], 'pdf')
%%  

%Also, write out some summary data to a cumulative summary file
%Vec_sum{k}=[Azi, Ele]; % preferred direction in each condition
sprint_txt = ['%s'];
for i = 1 : 100
     sprint_txt = [sprint_txt, ' %1.3f'];    
end
%buff = sprintf(sprint_txt,FILE, unique_stim_type,P_anova_hori );
if length(unique_stim_type)==1
    buff = sprintf(sprint_txt,FILE,SpikeChan,unique_stim_type,P_anova,Vec_sum{1},DDI );
    outfile = ['Z:\LBY\Recording data\Qiaoqiao\3D_T_Tuning\3D_T_Tuning_V.txt']; 
elseif length(unique_stim_type)==2
    buff = sprintf(sprint_txt,FILE,SpikeChan,unique_stim_type,P_anova,Vec_sum{1},Vec_sum{2}, DDI );
    outfile = ['Z:\LBY\Recording data\Qiaoqiao\3D_T_Tuning3D_T_Tuning_VV.txt']; 
else
    buff = sprintf(sprint_txt,FILE,SpikeChan,unique_stim_type,P_anova,Vec_sum{1},Vec_sum{2},Vec_sum{3}, DDI );
    outfile = ['Z:\LBY\Recording data\Qiaoqiao\3D_T_Tuning\3D_T_Tuning_VVC.txt']; 
end
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
%    fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t Veb_DDI\t Vis_DDI\t Com_DDI\t Veb_var_term\t Vis_var_term\t Com_var_term\t');
%    fprintf(fid, 'FILE\t SPon\t');
    fileTag = ['FILE\t SpikeChan\t stim_type\t P_anova\t preferred_direction\t DDI\t '];
    fprintf(fid, fileTag);
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
return;
