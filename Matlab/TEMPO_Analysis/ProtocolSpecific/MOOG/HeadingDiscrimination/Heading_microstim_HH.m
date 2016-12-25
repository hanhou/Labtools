%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%   HH2013
%-----------------------------------------------------------------------------------------------------------------------

function Heading_microstim_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, batch_flag);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME, :);
temp_mask_status = data.moog_params(MASK_STATUS,:,MOOG);
temp_mask_radius = data.moog_params(MASK_RADIUS,:,MOOG);
temp_microstim = data.misc_params(MICROSTIM,:);
%temp_spike_rates = data.spike_rates(SpikeChan, :); 

trials = 1:length(temp_heading);		% a vector of trial indices
select_trials =  (trials >= BegTrial) & (trials <= EndTrial) ;
%select_trials =  ((trials >= BegTrial) & (trials<=540)) | ((trials >= 1081) & (trials <= EndTrial));

stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
total_trials = temp_total_trials( select_trials);
mask_status= temp_mask_status( select_trials );
mask_radius= temp_mask_radius( select_trials );
microstim = temp_microstim(select_trials );
%spike_rates = temp_spike_rates( select_trials);

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
unique_mask_status = munique(mask_status');
unique_mask_radius = munique(mask_radius');
unique_microstim = munique(microstim');

% decide whether its vary stim_type or amplitude or mask
if find(unique_mask_status == 1) > 1   % seldom used
    condition = mask_status;
    con_txt = 'mask_status';
elseif find(unique_microstim == 1) > 1
    condition = microstim;
    con_txt = 'mircostim';    
end
unique_condition = munique(condition');
if length(unique_motion_coherence) >1
    stim_type = motion_coherence;
    unique_stim_type = unique_motion_coherence;   
end 

one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_condition);
repetition = floor( length(heading)/one_repetition ); % take minimum repetition

%determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
LEFT = 1;
RIGHT = 2;
for i= 1 : length(total_trials) 
    temp = data.event_data(1,:,i + BegTrial-1);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = RIGHT;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = LEFT;
    else
     %   choice(i) = RIGHT;
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end
end

correct_rate = [];
yy=[];
for k = 1:length(unique_stim_type)
    count = 1;
    for n = 1:length(unique_condition)
        for i = 1:length(unique_heading)
             trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & (condition == unique_condition(n)) ) ;
             rightward_trials = (trials_select & (choice == RIGHT) );  % For 0 degree, this is not necessarily the "correct_trials" like the on-line version.
             
%              % Tolerance. HH20141222
%              if (i == 1) && (sum(rightward_trials) <= 1)
%                  rightward_trials(trials_select) = false;
%              end
%              
%              if (i == length(unique_heading)) && (sum(rightward_trials) >= sum(trials_select)-1)
%                  rightward_trials(trials_select) = true;
%              end
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Set the smallest and largest value to 0 and 1. HH20140609
             % CAUTION !! This is just for testing the reasonability of
             % different fitting methods and SHOULD NOT be used for routine
             % data analysis.
             
%              if i == 1 ;rightward_trials = []; end
%              if i == length(unique_heading) ;rightward_trials = (trials_select); end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
             correct_rate{n}(k,i) = 1*sum(rightward_trials) / sum(trials_select);
             fit_data_psycho_cum{n,k}(i, 1) = unique_heading( i );
             fit_data_psycho_cum{n,k}(i, 2) = correct_rate{n}(k,i);
             fit_data_psycho_cum{n,k}(i, 3) = sum(trials_select); 
             
             % for logistic regression later
             yy{k}(count,1) = unique_heading(i);
             yy{k}(count,2) = unique_condition(n);
             yy{k}(count,3) = sum(rightward_trials);   
             yy{k}(count,4) = sum(trials_select);  
             count = count + 1;    
         end
         
    end
end 
%%%%%% use Wichman's MLE method to estimate threshold and bias
for n = 1:length(unique_condition)
    for k = 1:length(unique_stim_type)    
%         wichman_psy = pfit(fit_data_psycho_cum{n,k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
%         Thresh_psy{n,k} = wichman_psy.params.est(2);
%         Bias_psy{n,k} = wichman_psy.params.est(1);
%         psy_perf{n,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
        
        % I have some problem in using pfit toolbox under X64 system, so I use
        % maximal likelihood fitting in stead. HH20130826
        [bbb,ttt] = cum_gaussfit_max1(fit_data_psycho_cum{n,k});
        Thresh_psy{n,k} = ttt;
        Bias_psy{n,k} = bbb;
        psy_perf{n,k} = [Bias_psy{n,k},Thresh_psy{n,k}];
    end
end

% %%%%%%%%% Logistic test  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%logit in matlab  
if (length(unique_condition) == 1)		% there are no microstim trials, so don't do logit
    LOGIST_REG = 0;
else 
    LOGIST_REG = 1;
end


if (LOGIST_REG == 1)
    thresh = []; slope = []; 
    
    for k = 1:length(unique_stim_type) 
        % use logit fit, probit fit might be better, need to compare sometime
        [b, dev, stats] = glmfit([yy{k}(:,1) yy{k}(:,2) yy{k}(:,1).*yy{k}(:,2)],[yy{k}(:,3) yy{k}(:,4)],'binomial','link','probit');
        % first, the no-stim case
        bias(k,1) = (norminv(0.5)-b(1))/b(2);		% 50 pct PD threshold
     %   threshold(k,1) = abs( (norminv(0.84)-b(1))/b(2) - bias(k,1) );   % this is actually equal to abs( norminv(0.84)/b(2) )
        threshold(k,1) = abs( norminv(0.84)/b(2) ); % 84% gaussian threshold
        
        % now, the stim case
        bias(k,2) = (norminv(0.5)-(b(1) + b(3)))/(b(2) + b(4));  
    %    threshold(k,2)=abs( (norminv(0.84)-(b(1) + b(3)))/(b(2) + b(4))-bias(k,2) ); % for probit fit, which is cumulative gaussian
        threshold(k,2) = abs( norminv(0.84)/(b(2)+b(4)) );
        
        % p value
        P_bias(k) = stats.p(3);  % P value for shift (bias)
        P_slope(k) = stats.p(4);	% P value for slope (threshold)
    end
else
    disp('no microstim');
    return;
end
bias(:,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot psychometric function here
symbo{1,1} = 'bo';    symbo{1,2} = 'ro';    symbo{1,3} = 'go'; % non-microstimulation trials
symbo{2,1} = 'b+';    symbo{2,2} = 'r+';    symbo{2,3} = 'g+'; % stimulation trials
fitline{1,1} = 'b-';    fitline{1,2} = 'r-';    fitline{1,3} = 'g-'; % non-microstimulation trials
fitline{2,1} = 'b--';    fitline{2,2} = 'r--';    fitline{2,3} = 'g--'; % stimulation trials

% I redefine colors here temporarily because HeTao can only do visual trials.
% symbo{1,1} = 'ko';   % symbo{1,2} = 'ro';    symbo{1,3} = 'go'; % non-microstimulation trials
% symbo{2,1} = 'ro';   %  symbo{2,2} = 'r+';    symbo{2,3} = 'g+'; % stimulation trials
% fitline{1,1} = 'k-';  %  fitline{1,2} = 'r-';    fitline{1,3} = 'g-'; % non-microstimulation trials
% fitline{2,1} = 'r--';  %  fitline{2,2} = 'r--';    fitline{2,3} = 'g--'; % stimulation trials



figure(2); clf;  
set(gcf,'color','white');
set(2,'Position', [50,50, 1000,740], 'Name', 'Heading Discrimination-Vestibular');
axes('position',[0.07,0.35, 0.5,0.6] ); hold on; box on;
line([min(unique_heading),max(unique_heading)],[0.5 0.5],'LineStyle','--','color','k');
line([0 0], [0 1],'LineStyle','--','color','k');

orient landscape;


% fit data with cumulative gaussian and plot both raw data and fitted curve
legend_txt = [];
xi = min(unique_heading) : 0.1 : max(unique_heading);
for n = 1:length(unique_condition)
    for k = 1:length(unique_stim_type)      
        set(plot(unique_heading, correct_rate{n}(k,:), symbo{n,unique_stim_type(k)},  xi, cum_gaussfit(psy_perf{n,k}, xi),  fitline{n,unique_stim_type(k)} ),...
            'LineWidth',3,'markersize',13,'markerfacecolor',symbo{n,unique_stim_type(k)}(1));
        xlabel('Heading Angle');   
        ylim([0,1]);
        ylabel('Rightward Choices');
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
        legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
        legend_txt{k*2} = [''];
    end
end

% HH20130829
cs=['b','r','m'];
text(min(unique_heading),0.9,FILE);
text(min(unique_heading),0.84,['rep = ' num2str(repetition)]);

for k = 1:length(unique_stim_type)
    text(max(unique_heading)/3,0.33 - 0.12*(k-1),sprintf('  \\it\\Delta_\\mu '),'color',symbo{2,unique_stim_type(k)}(1));
    if (P_bias(k)<0.05) ;col = symbo{2,unique_stim_type(k)}(1); else col = 'k'; end
    text(max(unique_heading)/3 ,0.336 - 0.12*(k-1),sprintf('       %.2f (%.2f)',abs(Bias_psy{2,k}-Bias_psy{1,k}),P_bias(k)),'color',col);
    %     text(max(unique_heading)/3,0.21,sprintf('\\itp_\\mu = \\rm%.2g',P_bias(k)),'color',cs((P_bias(k)<0.05)+1));
    text(max(unique_heading)/3,0.33 - 0.12*(k-1)-0.06,sprintf('  \\itr_\\sigma '),'color',symbo{2,unique_stim_type(k)}(1));
    if (P_slope(k)<0.05) ;col = symbo{2,unique_stim_type(k)}(1); else col = 'k'; end
    text(max(unique_heading)/3,0.336 - 0.12*(k-1)-0.06,sprintf('       %.2f (%.2f)',Thresh_psy{2,k}/Thresh_psy{1,k},P_slope(k)),'color',col);
    %     text(max(unique_heading)/3,0.09,sprintf('\\itp_\\sigma =\\rm %.2g',P_slope(k)),'color',cs((P_slope(k)<0.05)*2+1));
end

xlim([min(unique_heading)*1.1 max(unique_heading)*1.1]);

% output some text of basic parameters in the figure
axes('position',[0.1,0.01, 0.6,0.25] );
xlim( [0,50] );
ylim( [2,13] );
text(0,9,sprintf('%s,  Coh%% = %g,  rep = %g\n    ', FILE, unique_motion_coherence, repetition));
text(0,8, 'Stimtype       u0            sigma0         u1            sigma1             pu          psigma');

for n = 1:length(unique_condition)        
    for k = 1:length(unique_stim_type)  
        text(0,7-k, num2str(unique_stim_type(k)));  % non-microstim
        text(8+(n-1)*17,7-k,num2str(Bias_psy{n,k}) );
        text(18+(n-1)*17, 7-k,num2str(Thresh_psy{n,k}) );
        text(45,7-k,num2str(P_bias(k)) ); % p value for bias shift
        text(55,7-k,num2str(P_slope(k)) ); % p value for slope change (threshold)
    end
end
axis off;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sprint_txt = ['%s']; 
% for i = 1 : 1000 % this should be large enough to cover all the data that need to be exported
%      sprint_txt = [sprint_txt, ' %4.3f'];    
% end
% buff = sprintf(sprint_txt, FILE, Bias_psy{1,:}, Thresh_psy{1,:},Bias_psy{2,:}, Thresh_psy{2,:}, P_bias(:), P_slope );
% % buff = sprintf(sprint_txt, FILE, psy_bias_shift{1,2}(1)-psy_bias_shift{1,1}(1),psy_bias_shift{2,2}(1)-psy_bias_shift{2,1}(1),psy_bias_shift{3,2}(1)-psy_bias_shift{3,1}(1) );
% %outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\Psychome_combined2.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\microstim_timecourse_5.dat'];
% % buff = sprintf(sprint_txt, FILE, psy_bias_shift{2,1}, psy_bias_shift{2,2} );
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\microstim_timecourse.dat'];
% 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t       10_bias\t 20_bias\t 30_bias\t 10_thresh\t 20_thresh\t 30_thresh\t 11_bias\t 21_bias\t 31_bias\t 11_thresh\t 21_thresh\t 31_thresh\t 1_shfit_p\t 2_shfit_p\t 3_shfit_p\t 1_slope_p\t 2_slope_p\t 3_slope_p\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);


%%% calculate and plot shift over time
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
span = 3;  % calculate threshod every ? repeats;
slide = 1;  % slide threshod with increment of ? repeats;
BegTrial_shift = BegTrial;
EndTrial_shift = BegTrial_shift + span*one_repetition-1;
nn=0;

while EndTrial_shift <= EndTrial    
    nn = nn + 1;
    
    select_trials_shift = ( (trials >= BegTrial_shift) & (trials <= EndTrial_shift) );
    stim_type_shift = temp_stim_type( select_trials_shift );
    mask_status_shift = temp_mask_status( select_trials_shift );
    heading_shift = temp_heading( select_trials_shift );
    microstim_shift = temp_microstim( select_trials_shift );    
    motion_coherence_shift = temp_motion_coherence( select_trials_shift );    
       
    unique_stim_type_shift = munique(stim_type_shift');
    unique_mask_status_shift = temp_mask_status( select_trials_shift );
    unique_heading_shift = munique(heading_shift');
    unique_microstim_shift = munique(microstim_shift'); 
    unique_motion_coherence_shift = munique(motion_coherence_shift');
    
    if length(unique_motion_coherence) >1        
        stim_type_shift = motion_coherence_shift; 
        unique_stim_type_shift = unique_motion_coherence_shift;
    end  
    
    total_trials_shift = temp_total_trials( select_trials_shift);
    
    % choice_shift = choice( select_trials_shift );
    choice_shift = choice( select_trials_shift (BegTrial : end));  % This fixes the failure of shifting window when BegTrial > 1. HH20140506
    
	for k = 1:length(unique_stim_type)  
        for n = 1:length(unique_microstim)
            for i = 1:length(unique_heading)
                 trials_shift =logical( (heading_shift == unique_heading(i)) & (microstim_shift == unique_microstim_shift(n)) & (stim_type_shift == unique_stim_type_shift(k)) ) ;
                 correct_trials_shift = (trials_shift & (choice_shift == RIGHT) );     
                 correct_rate_shift = 1*sum(correct_trials_shift) / sum(trials_shift);   
    
                 fit_data_psycho_cum_shift{k,n}(i, 1) = unique_heading( i );  
                 fit_data_psycho_cum_shift{k,n}(i, 2) = correct_rate_shift;
                 fit_data_psycho_cum_shift{k,n}(i, 3) = sum(trials_shift);
            end
%             wichman_psy = pfit(fit_data_psycho_cum_shift{k,1},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
%             psy_thresh_shift(k,nn) = wichman_psy.params.est(2);

            [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum_shift{k,n}); % to save time, use a different fit method
            psy_bias_shift{k,n}(nn) = bb;
            psy_thresh_shift{k,n}(nn) = tt;
        end    
    end
    BegTrial_shift = BegTrial_shift + slide*one_repetition;
    EndTrial_shift = EndTrial_shift + slide*one_repetition; 
end

% now plot
% bias
axes('position',[0.65 0.7 0.3 0.25] );
for n = 1:length(unique_microstim)
	for k = 1:length(unique_stim_type)
        set(plot(psy_bias_shift{k,n}(:), fitline{n,unique_stim_type(k)}),'LineWidth',3);
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
       % semilogy(psy_thresh_shift(k,:), f{k});
        hold on;
        xlabel('Repetition');  
        ylabel('u');
        xlim([1, nn]);
	end
end

% threshold
axes('position',[0.65 0.35 0.3 0.25] );
for n = 1:length(unique_microstim)
	for k = 1:length(unique_stim_type)
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
        set(plot(psy_thresh_shift{k,n}(:), fitline{n,unique_stim_type(k)}),'LineWidth',3);
       % semilogy(psy_thresh_shift(k,:), f{k});
        hold on;
        xlabel('Repetition');  
        ylabel('threshold');
       xlim([1, nn]);
	end
end

set(findall(gcf,'FontSize',10),'FontSize',15);



% for k = 1:length(unique_stim_type)
%     fprintf('%g\t %g\t %g\t %g\t %g\t %g\t %g\n',repetition, Bias_psy{1,k},Thresh_psy{1,k},Bias_psy{2,k},Thresh_psy{2,k},P_bias(k),P_slope(k));
% end


%% Data Saving

% Reorganized. HH20141124
config.batch_flag = batch_flag;

%%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = PackResult(FILE, SpikeChan, repetition, unique_stim_type, ... % Obligatory!!
    Bias_psy,Thresh_psy,P_bias,P_slope,psy_bias_shift,psy_thresh_shift);
config.save_figures = gcf;

config.suffix = 'ustim';


% Only once
config.sprint_once_marker = 'g';
config.sprint_once_contents = 'result.repetition';
% Loop across each stim_type
config.sprint_loop_marker = {'ggggggss'};
config.sprint_loop_contents = {['result.Bias_psy{1,k},result.Thresh_psy{1,k},result.Bias_psy{2,k},result.Thresh_psy{2,k},result.P_bias(k),result.P_slope(k),',...
    'num2str([result.psy_bias_shift{k,1}(:)'' result.psy_bias_shift{k,2}(:)'']),num2str([result.psy_thresh_shift{k,1}(:)'' result.psy_thresh_shift{k,2}(:)'' ])']};         

config.append = 1; % Overwrite or append
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveResult(config, result);



% %%%%%%%%%%%%%%%%%%%%%%%  Output   HH20140510 / HH20140621 %%%%%%%%%%%%%%%%%
% 
% if ~isempty(batch_flag)  % Figures
%     
%     %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     suffix = ['ustim'];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     outpath = ['Z:\Data\Tempo\Batch\' batch_flag(1:end-2) '\'];
%    
%     if ~exist(outpath,'dir')
%         mkdir(outpath);    
%     end
%     
%     % Save figures
%     orient landscape;
%     savefilename = [outpath [FILE '_' num2str(SpikeChan)] '_' suffix '.png'];
% 
%     if exist(savefilename,'file')
%         delete(savefilename);
%     end
%     saveas(gcf,savefilename,'png');
%         
% end
% 
% 
% % Print results
% 
% %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Only once
% sprint_once_marker_temp = 'g';
% sprint_once_contents = 'repetition';
% % Loop across each stim_type
% sprint_loop_marker_temp = 'ggggggss';  
% sprint_loop_contents = ['Bias_psy{1,k},Thresh_psy{1,k},Bias_psy{2,k},Thresh_psy{2,k},P_bias(k),P_slope(k),',...
%     'num2str([psy_bias_shift{k,1}(:)'' psy_bias_shift{k,2}(:)'']),num2str([psy_thresh_shift{k,1}(:)'' psy_thresh_shift{k,2}(:)'' ])']; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% sprint_once_marker = [];
% for i = 1:length(sprint_once_marker_temp)
%     sprint_once_marker = [sprint_once_marker '%' sprint_once_marker_temp(i) '\t '];
% end
% 
% sprint_loop_marker = [];
% for i = 1:length(sprint_loop_marker_temp)
%     sprint_loop_marker = [sprint_loop_marker '%' sprint_loop_marker_temp(i) '\t '];
% end
% 
% if ~isempty(batch_flag)  % Print to file
%     
%     outfile = [outpath suffix '.dat'];
%     printHead = 0;
%     if (exist(outfile, 'file') == 0)   % file does not yet exist
%         printHead = 1;
%     end
%     
%     fid = fopen(outfile, 'a');
%     % This line controls the output format
% 
%     if (printHead)
%         fprintf(fid, ['FILE\t ' sprint_once_contents '\t' sprint_loop_contents]);
%         fprintf(fid, '\r\n');
%     end
%     
%     fprintf(fid,'%s\t',[FILE '_' num2str(SpikeChan)]);
%     
% else  % Print to screen
%     fid = 1;  
% end
% 
% % Print once
% if ~isempty(sprint_once_marker_temp)
%     eval(['buff = sprintf(sprint_once_marker,' sprint_once_contents ');']);
%     fprintf(fid, '%s', buff);
% end
% 
% % Print loops
% for conditions = 1:3 % Always output 3 conditions (if not exist, fill with NaNs)
%     if sum(unique_stim_type == conditions)==0
%         buff = sprintf(sprint_loop_marker,ones(1,sum(sprint_loop_marker=='%'))*NaN);
%     else
%         k = find(unique_stim_type == conditions);
%         eval(['buff = sprintf(sprint_loop_marker,' sprint_loop_contents ');']);
%     end
%     fprintf(fid, '%s', buff);
% end
% 
% fprintf(fid, '\r\n');
% 
% if ~isempty(batch_flag)  % Close file
%     fclose(fid);
% end


% %%%%% Eye traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for k = 1:length(unique_stim_type)
    for n = 1:length(unique_condition)
        for i = 1:length(unique_heading)
            select = find( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & (condition == unique_condition(n)) ) ;
            for e = 1 : length(select) 
%                 DC_L = median( data.eye_data(1,600:620,select(e)) ); % only use left eye temporally        
%                 L_trace_trial{k,n,i}(e,:) = data.eye_data(1,610:700,select(e)) - DC_L; % show saccade period               
                  DC_L = median( data.eye_data(3,300:320,select(e)) ); % only use right eye temporally          % HH20130826
                  L_trace_trial{k,n,i}(e,:) = data.eye_data(3,400:500,select(e)) - DC_L; % show saccade period               

            end
        end
    end
end

% plot now
for k = 1 : length(unique_stim_type)
    for i = 1 : length(unique_heading)    
        axes('position',[0.6+(k-1)*0.12,0.85-(i-1)*0.105, 0.12,0.095] );
        plot( L_trace_trial{k,1,i}(:,:)','b' ); % nonstimulated trials
        hold on;
        plot( L_trace_trial{k,2,i}(:,:)','r' ); % stimulated trials
        %xlim([20, 65]); % should be a reasonable rangle
        %ylim([-2,2]);
        xlim([0 100]);  % HH20130826
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);   
        if k==1
            ylabel(num2str(unique_heading(i)));
        end
        if k==1 & i ==1
            title('ves');
        elseif k==2 & i ==1
            title('vis');
        elseif k==3 & i==1
            title('com');
        end
    end
end
%} 
%---------------------------------------------------------------------------------------
return;