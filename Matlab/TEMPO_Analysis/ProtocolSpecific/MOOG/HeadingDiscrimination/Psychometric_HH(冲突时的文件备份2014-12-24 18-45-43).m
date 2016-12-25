%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------
%% HH20141026

function Psychometric_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);

TEMPO_Defs;
Path_Defs;

%get the column of values for azimuth and elevation and stim_type
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME, :);
temp_mask_status = data.moog_params(MASK_STATUS,:,MOOG);

trials = 1:length(temp_heading);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
total_trials = temp_total_trials( select_trials);
mask_status= temp_mask_status( select_trials );

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_motion_coherence = munique(motion_coherence');
unique_mask_status = munique(mask_status');

if length(unique_motion_coherence)==1
    one_repetition = length(unique_heading)*length(unique_stim_type);
else
    one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_motion_coherence)-length(unique_heading);
end
repetitionN = floor( length(heading)/one_repetition ); % take minimum repetition

% whether to plot performance over time
overtimeplot = 1;  % compute and plot
%overtimeplot = 0;  % not compute and plot

%determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)
LEFT = 1;
RIGHT = 2;

% for i= 1 : length(total_trials)
%     temp = data.event_data(1,:,i + BegTrial-1);
%     events = temp(temp>0);  % all non-zero entries
%     if (sum(events == IN_T1_WIN_CD) > 0)
%         choice(i) = RIGHT;
%     elseif (sum(events == IN_T2_WIN_CD) > 0)
%         choice(i) = LEFT;
%     else
%      %   choice(i) = RIGHT;
%         disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
%     end
% end

% The previous one was awful. HH20140522
event_in_bin = squeeze(data.event_data(:,:,select_trials))';  % TrialNum * 5000

choice_per_trial = LEFT * squeeze(sum(event_in_bin == IN_T2_WIN_CD,2))' + RIGHT * squeeze(sum(event_in_bin == IN_T1_WIN_CD,2))';
if length(unique(choice_per_trial)) > 2  % This is safer
    disp('Neither T1 or T2 chosen / More than one target chosen.  This should not happen! File must be bogus.');
    fprintf('%g cases\n..', length(unique(choice_per_trial)));
    keyboard;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(FILE,'m5c174r1')
        choice_per_trial(choice_per_trial==3) = 2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
        for i = 1:length(unique_heading)
            if unique_stim_type(k) == 1 % for vestibular condition, take all the data regardless of visual coherence
                trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) ) ;
            else
                trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & (motion_coherence==unique_motion_coherence(c)) ) ;
            end
            
            rightward_trials = (trials_select & (choice_per_trial == RIGHT) );
            rightward_rate = 1*sum(rightward_trials) / sum(trials_select);
            
            fit_data_psycho_cum{c,k}(i, 1) = unique_heading(i);
            fit_data_psycho_cum{c,k}(i, 2) = rightward_rate;
            fit_data_psycho_cum{c,k}(i, 3) = sum(trials_select);
        end
        %          halfheading = length(unique_heading_nonzero/2);
        %          for j = 1: halfheading
        %              trials_left = find( (heading==unique_heading_nonzero(halfheading+1-j)) & (choice==LEFT) & (stim_type==unique_stim_type(k))  ) ;
        %              trials_right  = find( (heading==unique_heading_nonzero(halfheading+j)) & (choice==RIGHT) & (stim_type==unique_stim_type(k))  ) ;
        %              trials_all = find( ((heading==unique_heading_nonzero(halfheading+1-j)|(heading==unique_heading_nonzero(halfheading+j)) & (stim_type==unique_stim_type(k)) );
        %              correct_rate(k,j) = (length(trials_right)+length(trials_left))/length(trials_all);
        %              % for later weibull fit
        %              fit_valid_weibull{c,k}(j,1) = unique_heading_nonzero(halfheading+j);
        %              fit_valid_weibull{c,k}(j,2) = correct_rate(k,j);
        %              fit_valid_weibull{c,k}(j,3) = fit_data_psycho_cum{c,k}(aa,3);
        %          end
        
        % the correct rate does not take coherence into account,temporarily 05/29/09
        trials_rightward = find( (heading > 0) & (choice_per_trial==RIGHT) & (stim_type==unique_stim_type(k))  ) ;
        trials_leftward  = find( (heading < 0) & (choice_per_trial==LEFT) & (stim_type==unique_stim_type(k))  ) ;
        trials_all = find( ((heading < 0)|(heading > 0)) & (stim_type==unique_stim_type(k)) ); %exclude 0 headings
        correct_proportion(k) = (length(trials_rightward)+length(trials_leftward))/length(trials_all);
        
        aa = find(fit_data_psycho_cum{c,k}(:,2)>-99); % sometime it could be NaN due to the absence of that heading conditions
        fit_valid{c,k}(:,1) = fit_data_psycho_cum{c,k}(aa,1);
        fit_valid{c,k}(:,2) = fit_data_psycho_cum{c,k}(aa,2);
        fit_valid{c,k}(:,3) = fit_data_psycho_cum{c,k}(aa,3);
        
        %          % for later weibull fit use
        %          fit_valid_weibull{c,k}(:,1) = unique_heading( unique_heading>0) );
        %          fit_valid_weibull{c,k}(:,2) = correct_rate(k,:);
        %          fit_valid_weibull{c,k}(:,3) = fit_data_psycho_cum{c,k}(aa,3);
    end
end

for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
        %  Maximum likelihood fitting
        [bb,tt] = cum_gaussfit_max1(fit_valid{c,k});
        Thresh_psy{c,k} = tt;
        Bias_psy{c,k} = bb;
        psy_perf{c,k} =[bb,tt];
    end
end

% added by GY 12-04-07
% now this is the prediction when there are three stimuli conditions
if length(unique_motion_coherence)==1 && length(unique_stim_type) ==3
    Thresh_pred = sqrt( Thresh_psy{1}^2*Thresh_psy{2}^2/(Thresh_psy{1}^2+Thresh_psy{2}^2) );
    Bias_pred = (Bias_psy{1} * Thresh_psy{2}^2 + Bias_psy{2} * Thresh_psy{1}^2) / (Thresh_psy{1}^2 + Thresh_psy{2}^2);
end
% % this is the output, you can use it for plot of example cells
% xi = min(unique_heading) : 0.1 : max(unique_heading);
% for k = 1:length(unique_stim_type)
%     yi{k} = cum_gaussfit(psy_perf{k}, xi);
% end
% if length(unique_stim_type) ==3
%     yi_pred = cum_gaussfit([Bias_psy{3},Thresh_pred], xi); % smoothed line for prediction with PSE at actual combined condition
% end


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot psychometric, neurometric, CP over time
% % run the slide threshold over time, see whether performance fluctuate across time
% not work for coherence, temporarily 05/29/09

n=0;

if overtimeplot == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(unique_heading)>= 8 % Normal conditions
        if repetitionN >= 15
            span = 10;  % calculate threshod every ? repeats;
        else
            span = 5;
        end
        slide = 1;  % slide threshod with increment of ? repeats;
    else  % Too few angles in training sessions 
        span = round(80 / length(unique_heading));
        slide = round(span/5);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    BegTrial_shift = BegTrial;
    EndTrial_shift = BegTrial_shift + span*one_repetition-1;
    
    while EndTrial_shift <= EndTrial
        n = n + 1;
        
        select_trials_shift = ( (trials(select_trials) >= BegTrial_shift) & (trials(select_trials) <= EndTrial_shift) );
        
        stim_type_shift = stim_type( select_trials_shift );
        motion_coherence_shift = motion_coherence (select_trials_shift);
        mask_status_shift = mask_status( select_trials_shift );
        heading_shift = heading( select_trials_shift );
        unique_mask_status_shift = mask_status( select_trials_shift );
        choice_per_trial_shift = choice_per_trial (select_trials_shift);
        
        unique_stim_type_shift = munique(stim_type_shift');
        
        if find(unique_mask_status == 1) > 1
            condition_shift = mask_status_shift;
            unique_condition_shift = unique_mask_status_shift;
        else
            condition_shift = stim_type_shift;
            unique_condition_shift = unique_stim_type_shift;
        end
        
        for c = 1:length(unique_motion_coherence) % different coherence level
            for k = 1:length(unique_condition_shift)
                for i = 1:length(unique_heading)
                    
                    if unique_stim_type(k) == 1 % for vestibular condition, take all the data regardless of visual coherence
                        trials_shift =logical( (heading_shift == unique_heading(i)) & (stim_type_shift==unique_stim_type(k)) ) ;
                    else
                        trials_shift =logical( (heading_shift == unique_heading(i)) & (stim_type_shift==unique_stim_type(k)) & (motion_coherence_shift==unique_motion_coherence(c)) ) ;
                    end
                    
                    rightward_trials = (trials_shift & (choice_per_trial_shift == RIGHT) );
                    rightward_rate = 1*sum(rightward_trials) / sum(trials_shift);
                    
                    fit_data_psycho_cum_shift{c,k}(i, 1) = unique_heading(i);
                    fit_data_psycho_cum_shift{c,k}(i, 2) = rightward_rate;
                    fit_data_psycho_cum_shift{c,k}(i, 3) = sum(trials_shift);
                    
                    %                     trials_shift =logical( (heading_shift == unique_heading(i)) & (condition_shift == unique_condition_shift(k)) ) ;
                    %
                    %                     correct_trials_shift = (trials_shift & (total_trials_shift == CORRECT) );
                    %                     % make 'S' curve by using the rightward choice for y-axis
                    %                     if sum(trials_shift)>0
                    %                         if ( unique_heading(i) < 0 )
                    %                             correct_rate_shift(i) = 1 - 1*sum(correct_trials_shift) / sum(trials_shift);
                    %                         else
                    %                             correct_rate_shift(i) = 1*sum(correct_trials_shift) / sum(trials_shift);
                    %                         end
                    %                     end
                    %                     Trials_num(i) = sum(trials_shift);
                end
                
                %                 aa = find(correct_rate_shift >-1 );
                %                 for j = 1:length(aa)
                %                     fit_data_psycho_cum_shift{c,k}(j, 1) = fit_data_psycho_cum{k}(aa(j), 1);
                %                     fit_data_psycho_cum_shift{c,k}(j, 2) = correct_rate_shift(aa(j));
                %                     fit_data_psycho_cum_shift{c,k}(j, 3) = Trials_num(aa(j));
                %                 end
                %                 % this fixes a strange error: cum_gaussfit/pfit sometimes fail when pct choices are all 0's or 1's -CRF 8-13-08
                %                 if fit_data_psycho_cum_shift{c,k}(:,2)==0 | fit_data_psycho_cum_shift{c,k}(:,2)==1
                %                     fit_data_psycho_cum_shift{c,k}(fit_data_psycho_cum_shift{c,k}==0) = 0.001;
                %                     fit_data_psycho_cum_shift{c,k}(fit_data_psycho_cum_shift{c,k}==1) = 0.999;
                %                 end
                
                [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum_shift{c,k}); % to save time, use a different fit method
                psy_thresh_shift(c,k,n) = tt;
                psy_bias_shift(c,k,n) = bb;  % added Bias, CRF 11-5-09
                
            end
        end
        BegTrial_shift = BegTrial_shift + slide*one_repetition;
        EndTrial_shift = BegTrial_shift + span*one_repetition-1;
    end
end



%%  Plot psychometric function

symbo{1,1} = 'bo';    symbo{1,2} = 'ro';    symbo{1,3} = 'go';
symbo{2,1} = 'b^';    symbo{2,2} = 'r^';    symbo{2,3} = 'g^';
symbo{3,1} = 'bv';    symbo{3,2} = 'rv';    symbo{3,3} = 'gv';

fitline{1,1} = 'b-';    fitline{1,2} = 'r-';    fitline{1,3} = 'g-';
fitline{2,1} = 'b--';    fitline{2,2} = 'r--';    fitline{2,3} = 'g--';
fitline{3,1} = 'b:';    fitline{3,2} = 'r:';    fitline{3,3} = 'g:';

figure(2); clf;
set(2,'Position', [63 51 1080 500], 'Name', 'Heading Discrimination-Vestibular');
axes('position',[0.10 0.113 0.377 0.754] );
% fit data with cumulative gaussian and plot both raw data and fitted curve
legend_txt = [];

xi = min(unique_heading) : 0.05 : max(unique_heading);  noteN = 0;
for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
        plot(unique_heading, fit_valid{c,k}(:,2), symbo{c,unique_stim_type(k)},  xi, cum_gaussfit(psy_perf{c,k}, xi),  fitline{c,unique_stim_type(k)} ,...
            'Linewidth',3,'MarkerSize',9,'MarkerFaceColor',fitline{c,unique_stim_type(k)}(1));
        xlabel('Heading Angles');
        ylim([0,1]);
        ylabel('Rightward Choices');
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
        hold on;
        legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
        legend_txt{k*2} = [''];
        
        % Annotation
        noteN = noteN+1;
        
        if unique_stim_type(k) >= 2 % show coherencce
            text(max(unique_heading)*0.3,0.35-noteN*0.07,sprintf('%5.2f  %5.2f  (%g%%)',Bias_psy{c,k} ,Thresh_psy{c,k},unique_motion_coherence(c)),'color',symbo{c,unique_stim_type(k)}(1));
        else
            text(max(unique_heading)*0.3,0.35-noteN*0.07,sprintf('%5.2f  %5.2f',Bias_psy{c,k} ,Thresh_psy{c,k}),'color',symbo{c,unique_stim_type(k)}(1));
        end
        
    end
    
end

text(max(unique_heading)*0.3,0.35,sprintf('    \\mu       \\sigma'));

if length(unique_motion_coherence)==1 && length(unique_stim_type) ==3
    plot(xi, cum_gaussfit([Bias_pred, Thresh_pred], xi),  'g--' ,...
        'Linewidth',3);
    text(max(unique_heading)*0.3,0.35-(noteN+1)*0.07,sprintf('%5.2f  %5.2f  (pred)',Bias_pred ,Thresh_pred),'color','g');
end

%------ Output some text of basic parameters in the figure
% Motion parameters
params = data.moog_params(:,1,2);


axes('position',[0.062 0.9 0.895 0.148] );
xlim( [0,50] );
ylim( [2,10] );
text(0, 5, sprintf('%s,    reps = %g,   amp = %g m,   Nsigma = %g,   duration = %g ms', FILE,repetitionN,params(AMPLITUDE), params(NUM_SIGMAS),params(DURATION)));


% text(15,10,'coherence =');
% text(30,10,'repeats =');
%text(45,10,'maskradius =');
% text(25,10,num2str(unique_motion_coherence) );
% text(40,10,num2str(repetitionN) );
%text(55,10,num2str(unique_mask_radius) );
% text(10,8, 'u                   sigma             correct rate');

% for c = 1:length(unique_motion_coherence) % different coherence level
%     for k = 1:length(unique_stim_type)
%         text(0,8-k-(c-1)*3, num2str(unique_stim_type(k)));  % non-microstim
%         text(10,8-k-(c-1)*3,num2str(Bias_psy{c,k}) );
%         text(20,8-k-(c-1)*3,num2str(Thresh_psy{c,k}) );
%         text(30,8-k-(c-1)*3,num2str(correct_proportion(k)) );
%     end
% end

axis off;

%% Plot psycho over time

% 1. Threshold
a1 = axes('position',[  0.5620    0.5460    0.0580    0.3280] );
a2 = axes('position',[ 0.6280 0.546 0.317 0.3280]);

noteN = 0;

for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
        %         plot(3:n+2,psy_thresh_shift(k,:), fitline{1,k});
        %        semilogy(psy_thresh_shift(k,:), f{k});
        
        axes(a1);
        bar(fix(span/2)-0.3+noteN*0.2,Thresh_psy{c,k},0.2,'facecolor',fitline{c,unique_stim_type(k)}(1),'edgecolor','none');
        hold on;
        
        %%%%%%%%%%%%%%%%%%%%%%% add by 6/18/2010
        if n > 0
            axes(a2)
            m_threshold = squeeze(psy_thresh_shift(c,k,:));
            index = find(m_threshold > 150);
            m_threshold(index) = 150;
            plot((1:n) + fix(span/2),m_threshold, ['s' fitline{c,unique_stim_type(k)}],'LineWidth',2);
            hold on;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        noteN = noteN + 1;
    end
end

% Plot predicted threshold
if length(unique_motion_coherence)==1 && length(unique_stim_type) ==3
    axes(a1);
    bar(fix(span/2)-0.2+noteN*0.2,Thresh_pred,0.2,'facecolor','none','edgecolor',fitline{c,unique_stim_type(k)}(1),...
        'LineWidth',2,'LineStyle','--');
    
    if n > 0
        axes(a2);
        thres_pred_shift = squeeze(sqrt(psy_thresh_shift(c,1,:).^2 .* psy_thresh_shift(c,2,:).^2./(psy_thresh_shift(c,1,:).^2 + psy_thresh_shift(c,2,:).^2)));
        plot((1:n) + fix(span/2), thres_pred_shift, 'g--','LineWidth',2);
    end
end

axes(a2); grid on;

set(gca,'XTick',fix(span/2): max(1,fix(n/6)) : n+fix(span/2));
ll = get(gca,'XTickLabel');

if n>1 xlabel(sprintf('Repetition (span = %g)',span)); end
xlim([fix(span/2)+0.4 max(n,0)+fix(span/2)+0.5]);
ylims = ylim;
ylim( [0 ylims(2)*1.2] );
set(gca,'yTicklabel',[]);

axes(a1); grid on
set(gca,'XTick',[]); xlabel('All');
axis tight;
ylim( [0 ylims(2)*1.2] );
ylabel('Threshold');



% 2. Bias, CRF 11-5-09
noteN = 0;

a1 = axes('position',[  0.5620    0.11    0.0580    0.32] );
a2 = axes('position',[ 0.6280 0.11 0.317 0.32]);

for c = 1:length(unique_motion_coherence) % different coherence level
    for k = 1:length(unique_stim_type)
        %         plot(3:n+2,psy_bias_shift(k,:), fitline{1,k});
        % semilogy(psy_thresh_shift(k,:), f{k});
        axes(a1);
        
        bar(fix(span/2)-0.3+noteN*0.2,Bias_psy{c,k},0.2,'facecolor',fitline{c,unique_stim_type(k)}(1),'edgecolor','none');
        hold on;
        
        %%%%%%%%%%%%%%%%%%%%%%%add by 6/18/2010
        if n > 0
            axes(a2);

            m_bias = squeeze(psy_bias_shift(c,k,:));
            m_bias(m_bias > 100) = 100;
            m_bias(m_bias < -100) = -100;
            plot((1:n) + fix(span/2), m_bias, ['s' fitline{c,unique_stim_type(k)}],'LineWidth',2); hold on;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        noteN = noteN + 1;
        
    end
end

% Plot predicted bias
if length(unique_motion_coherence)==1 && length(unique_stim_type) ==3
    axes(a1);
    bar(fix(span/2)-0.2+noteN*0.2,Bias_pred,0.2,'facecolor','none','edgecolor',fitline{c,unique_stim_type(k)}(1),...
        'LineWidth',2,'LineStyle','--');
    
    if n > 0
        axes(a2);
        Bias_pred_shift = squeeze((psy_bias_shift(c,1,:).* psy_thresh_shift(c,2,:).^2 + psy_bias_shift(c,2,:).* psy_thresh_shift(c,1,:).^2)./(psy_thresh_shift(c,2,:).^2 + psy_thresh_shift(c,1,:).^2));
        plot((1:n) + fix(span/2), Bias_pred_shift, 'g--','LineWidth',2);
    end
end


axes(a2); grid on;
set(gca,'XTick',fix(span/2): max(1,fix(n/6)) :n+fix(span/2));
ll = get(gca,'XTickLabel');

if n>1 xlabel(sprintf('Repetition (span = %g)',span)); end
set(gca,'yTicklabel',[]);
plot(xlim,[0 0],'k:');

xlim([fix(span/2)+0.4 max(n,0)+fix(span/2)+0.5]);
ylims = ylim;
ylim([-max(abs(ylims))*1.1 max(abs(ylims))*1.1]);

axes(a1); grid on;
set(gca,'XTick',[]); xlabel('All');
axis tight;
ylim([-max(abs(ylims))*1.1 max(abs(ylims))*1.1]);

ylabel('Bias');
plot(xlim,[0 0],'k:');


orient tall;

SetFigure(15);


%% Data Saving

%%%%%%%%%%%%%%%%%%%%%  Output   HH20140510 / HH20140621 / HH20141003 %%%%%%%%%%%%%%%%%

if ~isempty(batch_flag)  % Figures and raw data (always in "result" structure)
    
    %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    suffix = ['Psycho'];
    
    % Figures to save
    save_figures = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    outpath = ['Z:\Data\Tempo\Batch\' batch_flag(1:end-2) '\'];
    
    % Check directory
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
    savefilename = [outpath [FILE '_' num2str(SpikeChan)] '_' suffix];
    
    % Delete existing data files
    if exist([savefilename '.mat'],'file')
        delete([savefilename '*.*']);
    end
    
    % Save raw data
    %     save(savefilename,'result');
    
    % Save figures
    for ff = 1:length(save_figures)
        %         orient landscape;
        set(save_figures(ff),'Visible','on');
        print(save_figures(ff),'-dbitmap',[savefilename '_fig_' num2str(ff) '.bmp']);
        close(save_figures(ff));
        %         saveas(save_figures(ff),[savefilename '_fig_' num2str(ff)],'bmp');
    end
    
end


% Print part of data to texts (clipboard or .dat file)

%%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only once
sprint_once_marker_temp = 'gg';
sprint_once_contents = 'repetitionN, unique_motion_coherence(1)';  % We only output the first coherence
% Loop across each stim_type
% sprint_loop_marker_temp = {'sss';
%                        };
% sprint_loop_contents = {'num2str(result.PSTH{1,2,1}.ys((k-1)*2+1,:)), num2str(result.PSTH{1,2,1}.ys((k-1)*2+2,:)), num2str(result.PSTH{1,2,1}.ps(k,:))';
%                        };

if overtimeplot == 1 && exist('psy_thresh_shift','var')
    sprint_loop_marker_temp = {'gg';
        'ss'};
    sprint_loop_contents = {'Thresh_psy{1,k}, Bias_psy{1,k} ';
        'num2str(squeeze(psy_thresh_shift(1,k,:)))'', num2str(squeeze(psy_bias_shift(1,k,:)))'''};
else
    sprint_loop_marker_temp = {'gg'; 'ss'};
    sprint_loop_contents = {' Thresh_psy{1,k}, Bias_psy{1,k} '; 'NaN, NaN'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HD_rep	HD_vest_psy_thres	HD_vis_psy_thres	HD_comb_psy_thres	HD_vest_CP_grand_center	HD_vest_CP_grand_sac	HD_vis_CP_grand_center	HD_vis_CP_grand_sac	HD_comb_CP_grand_center	HD_comb_CP_grand_sac


sprint_once_marker = [];
for i = 1:length(sprint_once_marker_temp)
    sprint_once_marker = [sprint_once_marker '%' sprint_once_marker_temp(i) '\t '];
end

if ~isempty(batch_flag)  % Print to file
    
    outfile = [outpath suffix '.dat'];
    printHead = 0;
    if (exist(outfile, 'file') == 0)   % file does not yet exist
        printHead = 1;
    end
    
    fid = fopen(outfile, 'a');
    % This line controls the output format
    
    if (printHead)
        fprintf(fid, ['FILE\t ' sprint_once_contents '|\t']);
        
        for ll = 1:length(sprint_loop_contents)
            fprintf(fid,[sprint_loop_contents{ll} '|\t']);
        end
        fprintf(fid, '\r\n');
    end
    
    fprintf(fid,'%s\t',[FILE '_' num2str(SpikeChan)]);
    
else  % Print to screen
    fid = 1;
end

toClip = [];

% Print once
if ~isempty(sprint_once_marker_temp)
    eval(['buff = sprintf(sprint_once_marker,' sprint_once_contents ');']);
    fprintf(fid, '%s', buff);
    toClip = [toClip sprintf('%s', buff)];
end

% Print loops
for ll = 1:length(sprint_loop_contents)
    
    sprint_loop_marker = [];
    for i = 1:length(sprint_loop_marker_temp{ll})
        sprint_loop_marker = [sprint_loop_marker '%' sprint_loop_marker_temp{ll}(i) '\t '];
    end
    
    for conditions = 1:3 % Always output 3 conditions (if not exist, fill with NaNs)
        if sum(unique_stim_type == conditions)==0
            buff = sprintf(sprint_loop_marker,ones(1,sum(sprint_loop_marker=='%'))*NaN);
        else
            k = find(unique_stim_type == conditions);
            eval(['buff = sprintf(sprint_loop_marker,' sprint_loop_contents{ll} ');']);
        end
        fprintf(fid, '%s', buff);
        toClip = [toClip sprintf('%s', buff)];
    end
    
end

fprintf(fid, '\r\n');
toClip = [toClip sprintf('\r\n')];
clipboard('copy',toClip);

if ~isempty(batch_flag)  % Close file
    fclose(fid);
end
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sprint_txt = ['%s'];
% for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
%      sprint_txt = [sprint_txt, ' %4.3f\t'];
% end
% %buff = sprintf(sprint_txt, FILE, fit_data_psycho_cum{1}(:,1),fit_data_psycho_cum{1}(:,2),fit_data_psycho_cum{1}(:,3),fit_data_psycho_cum{2}(:,2),fit_data_psycho_cum{2}(:,3),fit_data_psycho_cum{3}(:,2),fit_data_psycho_cum{3}(:,3) );
% %buff = sprintf(sprint_txt, FILE, unique_motion_coherence, Thresh_psy{1}, Thresh_psy{2},Thresh_psy{3} );
% if length(unique_stim_type)==3
%     buff = sprintf(sprint_txt, FILE, unique_motion_coherence, unique_stim_type, Thresh_psy{1,1},Thresh_psy{1,2},Thresh_psy{1,3});
% elseif length(unique_stim_type)==2
%     buff = sprintf(sprint_txt, FILE, unique_motion_coherence, unique_stim_type,Thresh_psy{1,1},Thresh_psy{1,2});
% else
%     buff = sprintf(sprint_txt, FILE, unique_motion_coherence, unique_stim_type,Thresh_psy{1,1});
% end
%
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\Psychome_combined.dat'];
% %outfile = ['Z:\Users\Yong\psy.dat'];
% %buff = sprintf(sprint_txt, FILE, unique_heading', fit_valid{1}(:,2), fit_valid{2}(:,2),fit_valid{3}(:,2) );
% % buff = sprintf(sprint_txt, FILE, Bias_psy{1,1},Bias_psy{1,2},Bias_psy{1,3});
% % outfile = ['Z:\Users\Yong\inactivation_MSTdBias.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t          coherence\t  1_bias\t 2_bias\t 3_bias\t 1_thresh\t 2_thresh\t 3_thresh\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% %---------------------------------------------------------------------------------------
return;