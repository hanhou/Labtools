%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
% shift analyze window across the stimulus duration
%--	02/22/06 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cum_shiftwindow_HH(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);

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
temp_spike_data = data.spike_data(SpikeChan,:);   % spike rasters
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
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
disc_heading = unique_heading( floor(length(unique_heading)/2)+1 : end );

% remove trials selected
Discard_trials = find(trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );
spike_data( find(spike_data>100) ) = 1; % something is absolutely wrong

spike_rates_copy = spike_rates; % save a copy for spike_rates

% shift window, currently using 100ms shift with 1000ms width
window_interval = 1000;
step=100;
endloop = floor((3000-window_interval/2)/100);

if isempty(batch_flag) ;progressbar('ShiftWindows'); end

for s = 1 : endloop 
    % replace spike_rates with spike_data based on analize window set
    for ss =  1 : length(spike_rates) % ss marks the index of trial
         spike_rates(ss) = sum( spike_data(1,StartEventBin(1)-0.5*window_interval+step*(s-1)+5000*(ss-1) : StartEventBin(1)-0.5*window_interval+window_interval+step*(s-1)+5000*(ss-1)) ) / 1.01; % 996~3006 every 1000ms
%        spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+500-50*(s-1)+5000*(ss-1) : StopEventBin(1)-500+50*(s-1)+5000*(ss-1)))/(0.001*(StopEventBin(1)-StartEventBin(1)-1000+100*(s-1)));%Enlarge the window
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %neurometric dataset and calculate ROC, Choice Probability(CP)
    %determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
    LEFT = 1;
    RIGHT = 2;
    for i= 1 : length(spike_rates) 
        temp = data.event_data(1,:,i + BegTrial-1);
        events = temp(temp>0);  % all non-zero entries
        if (sum(events == IN_T1_WIN_CD) > 0)
            choice(i) = RIGHT;
        elseif (sum(events == IN_T2_WIN_CD) > 0)
            choice(i) = LEFT;
        else
            disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
        end
    end
%     if FILE=='m2c384r2.htb'
%        choice(889) =2; % for cell m2c384r2 % for some reason the choice is 0 for
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % neural database
    one_repetition = length(unique_heading)*length(unique_stim_type);
    repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition
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
    Z_Spikes_Ori = Z_Spikes; % keep a Z_Spikes unchanged for later use
    % now group neuronal data into two groups according to monkey's choice
    for k = 1:length(unique_stim_type)    % notice that the condition is double than disc_heading    
        for i = 1:length(unique_heading)
            select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
            resp{k,i} = spike_rates(select);   
            resp_mat{k}(i) = mean(resp{k,i});  % the mean firing rate for each heading 
            resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);
            % calculate CP, group data based on monkey's choice 
            resp_left_choose{k,i} = spike_rates(select & (choice == LEFT) );
            resp_right_choose{k,i} = spike_rates(select & (choice == RIGHT) );
            if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)   % make sure each condition has at least 3 data values
          %  if (length(resp_left_choose{k,i}) / length(resp{k,i}) <0.25) |  (length(resp_left_choose{k,i}) / length(resp{k,i}) >0.75)  
                Z_Spikes(select) = 9999;   % similar to NaN, just make a mark            
           %     Z_Spikes( (heading == unique_heading(length(unique_heading)+1-i)) & (stim_type == unique_stim_type(k)) ) = 9999; % the corresponding heading is also excluded
                CP{k}(i) = 9999;
          %      CP{k}(length(unique_heading)+1-i) = 9999;
            else
                CP{k}(i) = 0;
            end 
        end  
        % now across all data 
        resp_left_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
        resp_right_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
        resp_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) ); 
        
    end

    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decide sign of slope, decide preferred direction 
    kk{1}='b.-'; kk{2}='r.-'; kk{3}='k.-'; kk{4}='g.-'; kk{5}='y.'; kk{6}='y.-'; 
    nu{1}='b.--'; nu{2}='r.--'; nu{3}='k.--'; nu{4}='g.--'; nu{5}='y.--'; nu{6}='y.--';  
    for k = 1 : length(unique_stim_type)
        % decide whether ves and vis is congruent tuning. Fit line by linear
        % regression first and compare the sign of each condition to decide whether
        % congruent or opposite, this is used to check whether congruent cells lead
        % to better neuronal performance in combined condition, and vice versa
        [rr,pp] = corrcoef(unique_heading, resp_mat{k}(:));
        line_re{k} = rr(1,2);
        line_p{k} = pp(1,2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
    % neurothreshold 
    fit_data_neuro = [];
    fit_data_neuro_cut = [];       
    for k = 1 : length(unique_stim_type)
        fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
        for i = 1 : length(unique_heading)-1   % subtract the 0 heading
            trials_n =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;            
            if i < (1+length(unique_heading))/2
                 % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction            
                 Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i+1},resp{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
             else
                 Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i}, resp{k,(i+1)},100 ); % compare to the 0 heading condition, which is straght ahead
             end
     %        if  resp_mat{k}(1) < resp_mat{k}(end)
             if line_re{k} > 0 % if left heading is not the larger responses, then linear regression will be positive 
                 % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
                 % here we asume if left response larger than right response then asign the left to be preferred direction
                 Neuro_correct{k}(i) = 1 - Neuro_correct{k}(i);            
             end  
             fit_data_neuro_cum{k}(i,2) = Neuro_correct{k}(i); 
             fit_data_neuro_cum{k}(i,3) = sum(trials_n);  % for later function fit use
         end
    end
    % choice probability
    for k = 1 : length(unique_stim_type)
        for i = 1 : length(unique_heading)  
            if CP{k}(i)~=9999
               CP{k}(i) = rocN( resp_left_choose{k,i},resp_right_choose{k,i},100 );
            else
               CP{k}(i) = NaN;
            end
            if (length(resp_left_all{k}) > 3) | (length(resp_right_all{k}) > 3)
            %if  (length(resp_left_all{k}) / length(resp_all{k}) >0.25) |  (length(resp_left_all{k}) / length(resp_all{k}) <0.75)
                CP_all{k} = rocN( resp_left_all{k},resp_right_all{k},100 );
            else
                CP_all{k} = NaN; 
            end
            if  line_re{k} > 0
                CP{k}(i) = 1 - CP{k}(i);
                CP_all{k} = 1 - CP_all{k};
            end
        end
    end
    % %--------------------------------------------------------------------------
    %do permutation to test the significance of CP_all{k}, re-calculate CP 1000 times
%     perm_num = 1000;
%     Z_Spikes_perm = Z_Spikes;
%     bin = 0.005;
%     x_bin = 0 : bin : 1;
%     Z_Spike_con = cell(1,perm_num);
%     
%     for k = 1:1    
%  %   for k = 1:length(unique_stim_type)
%         for n = 1 : perm_num
%             % temperarilly only use near-threshold heading angles where monkey make a guess mainly
%             select = logical( (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
%             Z_Spikes_con{k} = Z_Spikes_perm( select );
%             Z_Spikes_con{k} = Z_Spikes_con{k}(randperm(length(Z_Spikes_con{k})));   % permute spike_rates
%             Z_Spikes_perm(select) = Z_Spikes_con{k};    % now in spike_rates, the corresponding data were permuted already
% 
%             resp_left_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
%             resp_right_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
% 
%             if  (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
%                 CP_all_perm{k}(n) = rocN( resp_left_all_perm{k}, resp_right_all_perm{k},100 );
%             else
%                 CP_all_perm{k}(n) = NaN; 
%             end
%             if  line_re{k} > 0  
%                 CP_all_perm{k}(n) = 1 - CP_all_perm{k}(n);             
%             end  
%         end
%         % now calculate p value or significant test
%         if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
% 
%             hist_perm(k,:) = hist( CP_all_perm{k}(:), x_bin );  % for permutation
%             bin_sum = 0;
%             n = 0;
%             while ( n < (CP_all{k}/bin) )
%                  n = n+1;
%                  bin_sum = bin_sum + hist_perm(k, n);
%                  if CP_all{k} >= 0.5                  % note it's two tail test
%                     p{k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
%                  else
%                     p{k} = 2* bin_sum / perm_num;
%                  end
%             end
%         else
%             p{k} = NaN;
%         end 
%     end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% use Wichman's MLE method to estimate threshold and bias
   
    for k = 1:length(unique_stim_type)        
        
   %     fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
    %    fit_data_neuro_cum{k}(:,2) = Neuro_correct{k}(:);
%         if FILE=='m2c247r2.htb' & s==8   % for unkown reason, for cell m2c247r2, when s==8,k=2, the data cannot be fit 
%    %-9.0000         0   18.0000
%   % -3.4700    0.0139   18.0000
%    %-1.3300    0.1852   18.0000
%    %-0.5100    0.5556   19.0000
%    % 0.5100    0.4444   18.0000
%    % 1.3300    0.8148   18.0000
%    % 3.4700    0.9861   18.0000
%    % 9.0000    1.0000   18.0000
%             [bb,tt] = cum_gaussfit_max1(fit_data_neuro_cum{k});
%             Thresh_neu{k} = tt;
%         else
%             wichman_neu = pfit(fit_data_neuro_cum{k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
%             Thresh_neu{k} = wichman_neu.params.est(2);
            [bb,tt] = cum_gaussfit_max1(fit_data_neuro_cum{k});
            Thresh_neu{k} = tt;
            % negative and positive infinite value means flat tuning
            if Thresh_neu{k}<0 | Thresh_neu{k}> 300
                Thresh_neu{k} = 300;
%                wichman_neu.params.est(2) = 300;
            end
%             Bias_neu{k} = wichman_neu.params.est(1);
%             neu_perf{k} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
 %       end
    end
    for k = 1:length(unique_stim_type)
        % calculate cp and threshold based on different analyze window
        CP_all_shift(k,s) = CP_all{k};
%         p_shift(1,s) = p{1};
        [~,p_shift(k,s)] = ttest2(resp_left_all{k},resp_right_all{k}); % Use ttest2 instead of permutation for speed (the outcomes are quite similar). Hh(k,2)0140601
        Thresh_neu_shift(k,s) = Thresh_neu{k};
    end
    
    if isempty(batch_flag)  ;progressbar(s/endloop); end
end


% Hh(k,2)0140511 plot CP time course

stim_type_names = {'Vest','Vis','Comb'}; % stim_type = 0, 1, 2, 3
colors = {'k','r','g'};
xx = 0:step:step*(endloop-1);

set(figure(777),'color','w','position',[5 50 1400 700]);
if ~isempty(batch_flag); set(777,'visible','off'); end

clf;

for k = 1:length(unique_stim_type)
    if k==1
        [ax,h(k,1),h(k,2)] = plotyy(xx,CP_all_shift(k,:),xx,Thresh_neu_shift(k,:));
    else
        axes(ax(1)); hold on; h(k,1) = plot(xx,CP_all_shift(k,:));
        axes(ax(2)); hold on; h(k,2) = plot(xx,Thresh_neu_shift(k,:));
    end
    
    set(h(k,1),'linewidth',2,'marker','o','markersize',10,'color',colors{unique_stim_type(k)})
    set(h(k,2),'linewidth',2,'linestyle','--','marker','^','markersize',10,'color',colors{unique_stim_type(k)})
    
    axes(ax(1)); hold on;
    plot(xx(p_shift(k,:)<0.05),CP_all_shift(k,p_shift(k,:)<0.05),['o' colors{unique_stim_type(k)}],'markerfacecolor',colors{unique_stim_type(k)},'markersize',10);   % Sig. CP
end

set(ax,'xlim',[0 xx(end)]);

axes(ax(1)); hold on;
lims = axis; hold on;
plot([lims(1) lims(2)],[0.5 0.5],'k--');
ylim([0.3 1]);
set(ax(1),'xtick',[],'ytick',0.3:0.1:1,'ycolor','k');
ylabel(ax(1),'Grand CP')

axes(ax(2));
ylim([1 200]);
xlabel(ax(2),['Center of ' num2str(window_interval) ' ms time window (ms)']);
ylabel(ax(2),'Neuronal Threshold')
set(ax(2),'yscale','log','yticklabel',[1 10 100],'ytick',[1 10 100],'ycolor','k');

visDur = round(mean(StopEventBin-StartEventBin)/100)*100;
plot(ax(1),[visDur visDur],[0.3 1],'k');
title([FILE '_' num2str(SpikeChan)]);
SetFigure(15);

%%%%%%%%%%%%%%%%%%%%%%%  Batch Output.   Hh(k,2)0140510  %%%%%%%%%%%%%%%%%

if ~isempty(batch_flag)
    
    outpath = 'Z:\Data\Tempo\Batch Output\CP_shift_window\';
    
    % Save figures
    orient landscape;
    %set(gcf,'PaperUnits','inches','PaperSize',[3,6],'PaperPosition',[0 0 3 6]);
    saveas(gcf,[outpath FILE '_' num2str(SpikeChan) '.png'],'png');
    % print('-dpng','-r100',[outpath FILE '.png'])
    
    % Print results
    sprint_txt = '%s\t';
    for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
        sprint_txt = [sprint_txt, ' %4.3f\t'];
    end
    
    outfile = [outpath 'CP_shift.dat'];
    printHead = 0;
    if (exist(outfile, 'file') == 0)   % file does not yet exist
        printHead = 1;
    end
    
    % This line controls the output format
    buff = sprintf(sprint_txt, [FILE '_' num2str(SpikeChan)], xx, Thresh_neu_shift(1,:), CP_all_shift(1,:), p_shift(1,:));
    
    fid = fopen(outfile, 'a');
    if (printHead)
        fprintf(fid, 'FILE\t   NumWindows   Neurothres   CPs    ps  ');
        fprintf(fid, '\r\n');
    end
    
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Also, write out some summary data to a cumulative summary file
% sprint_txt = ['%s']; 
% for i = 1 : 1000 % this should be large enough to cover all the data that need to be exported
%      sprint_txt = [sprint_txt, ' %4.3f'];    
% end
% %buff = sprintf(sprint_txt, FILE, Thresh_neu_shift(1,:),Thresh_neu_shift(2,:),Thresh_neu_shift(3,:), CP_all_shift(1,:), CP_all_shift(2,:),CP_all_shift(3,:) );
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum_shiftwindowlong.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum_shiftwindow_cuecombined250ms.dat'];
% buff = sprintf(sprint_txt, FILE, Thresh_neu_shift(1,:),CP_all_shift(1,:) );
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum_shiftwindow_1000ms_cah.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)   % file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t         Ve_N_thr\t ves_CP\t ves_p\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% %---------------------------------------------------------------------------------------
return;