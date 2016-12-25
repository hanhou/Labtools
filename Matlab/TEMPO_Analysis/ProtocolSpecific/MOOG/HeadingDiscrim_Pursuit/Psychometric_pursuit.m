%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function Psychometric_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

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
temp_fix_x = data.moog_params(FIX_X,:,MOOG); 
pursuit_speed = data.moog_params(PURSUIT_VELOCITY,1,MOOG)

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
fix_x = temp_fix_x( select_trials);

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
unique_mask_status = munique(mask_status');
unique_mask_radius = munique(mask_radius');
unique_microstim = munique(microstim');
unique_fix_x = munique(fix_x');

one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_fix_x);
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
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end
end

correct_rate = [];
yy=[];
for n = 1:length(unique_fix_x)
    for k = 1:length(unique_stim_type)   
        for i = 1:length(unique_heading)
             trials_select =logical( (heading == unique_heading(i)) & (stim_type==unique_stim_type(k)) & (fix_x == unique_fix_x(n)) ) ;
             correct_trials = (trials_select & (choice == RIGHT) );
             correct_rate{n}(k,i) = 1*sum(correct_trials) / sum(trials_select);  
             fit_data_psycho_cum{n,k}(i, 1) = unique_heading( i );  
             fit_data_psycho_cum{n,k}(i, 2) = correct_rate{n}(k,i);
             fit_data_psycho_cum{n,k}(i, 3) = sum(trials_select);   
         end
    end
end 
%%%%%% use Wichman's MLE method to estimate threshold and bias
for n = 1:length(unique_fix_x)
    for k = 1:length(unique_stim_type)    
        wichman_psy = pfit(fit_data_psycho_cum{n,k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
        Thresh_psy{n,k} = wichman_psy.params.est(2);
        Bias_psy{n,k} = wichman_psy.params.est(1);
        psy_perf{n,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
        
%         [bbb,ttt] = cum_gaussfit_max1(fit_data_psycho_cum{n,k});
%         Thresh_psy{n,k} = ttt;
%         Bias_psy{n,k} = bbb;
%         psy_perf{n,k} = [Bias_psy{n,k},Thresh_psy{n,k}];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot psychometric function here
symbo{1} = 'bo';    symbo{2} = 'ro';    symbo{3} = 'go'; 
fitline{1} = 'b-';    fitline{2} = 'r-';    fitline{3} = 'g-'; 
figure(2);
set(2,'Position', [5,1 1000,740], 'Name', 'Heading Discrimination-Vestibular');
orient landscape;
% fit data with cumulative gaussian and plot both raw data and fitted curve
legend_txt = [];
xi = min(unique_heading) : 0.1 : max(unique_heading);
for n = 1:length(unique_fix_x)
    axes('position',[0.05+(n-1)*0.3,0.3, 0.28,0.4] );    
    for k = 1:length(unique_stim_type)      
        plot(unique_heading, correct_rate{n}(k,:), symbo{k},  xi, cum_gaussfit(psy_perf{n,k}, xi),  fitline{k} );
        xlabel('Heading Angles');   
        ylim([0,1]);
        ylabel('Rightward Choices');
        set(gca, 'YTickMode','auto');
        set(gca, 'xTickMode','auto');
        hold on;
        legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
        legend_txt{k*2} = [''];
    end
    axes('position',[0.1+(n-1)*0.3,0.72, 0.28,0.1] );
    xlim( [0,50] );
    ylim( [1,2] );
    text(0,1,'pursuit =');
    text(20,1,num2str(unique_fix_x(n)));
    axis off;
end

% output some text of basic parameters in the figure
axes('position',[0.1,0.8, 0.5,0.15] );
xlim( [0,50] );
ylim( [2,10] );
text(0, 10, FILE);
text(10,10,'coherence =');
text(25,10,'pursuit speed =');
text(40,10,'repeats =');
%text(45,10,'maskradius =');
text(20,10,num2str(unique_motion_coherence) );
text(35,10,num2str(pursuit_speed) );
text(48,10,num2str(repetition) );
%text(55,10,num2str(unique_mask_radius) );
text(0,8.5, 'Pursuit       Stimtype       bias       threshold');

for n = 1:length(unique_fix_x)     
    for k = 1:length(unique_stim_type) 
        text(0,8-(n-1)*length(unique_stim_type)-k, num2str(unique_fix_x(n))); 
        text(10,8-(n-1)*length(unique_stim_type)-k, num2str(unique_stim_type(k))); 
        text(20,8-(n-1)*length(unique_stim_type)-k,num2str(Bias_psy{n,k}) );
        text(30,8-(n-1)*length(unique_stim_type)-k,num2str(Thresh_psy{n,k}) );
    end
end
axis off;

%%% calculate and plot shift over time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
span = 5;  % calculate threshod every ? repeats;
slide = 5;  % slide threshod with increment of ? repeats;
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
    choice_shift = choice( select_trials_shift );
    
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
axes('position',[0.2,0.05, 0.5,0.18] );
for n = 1:length(unique_microstim)
	for k = 1:length(unique_stim_type)
        plot(psy_bias_shift{k,n}(:), fitline{n,k});
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
axes('position',[0.55,0.05,0.4,0.18] );
for n = 1:length(unique_microstim)
	for k = 1:length(unique_stim_type)
        plot(psy_thresh_shift{k,n}(:), fitline{n,k});
       % semilogy(psy_thresh_shift(k,:), f{k});
        hold on;
        xlabel('Repetition');  
        ylabel('threshold');
        xlim([1, nn]);
	end
end
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprint_txt = ['%s']; 
for i = 1 : 1000 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
%buff = sprintf(sprint_txt, FILE, Bias_psy{1,:}, Thresh_psy{1,:},Bias_psy{2,:}, Thresh_psy{2,:}, P_bias(:), P_slope );
if length(unique_fix_x)==3
   buff = sprintf(sprint_txt, FILE, unique_motion_coherence,pursuit_speed,Bias_psy{1,1},Bias_psy{2,1},Bias_psy{3,1},Thresh_psy{1,1},Thresh_psy{2,1},Thresh_psy{3,1});
%outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\Psychome_combined2.dat'];
   outfile = ['Z:\Users\Yong\headingpursuit.dat']; 

    % buff = sprintf(sprint_txt, FILE, psy_bias_shift{2,1}, psy_bias_shift{2,2} );
    % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\microstim_timecourse.dat'];

    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t       10_bias\t 20_bias\t 30_bias\t 10_thresh\t 20_thresh\t 30_thresh\t 11_bias\t 21_bias\t 31_bias\t 11_thresh\t 21_thresh\t 31_thresh\t 1_shfit_p\t 2_shfit_p\t 3_shfit_p\t 1_slope_p\t 2_slope_p\t 3_slope_p\t');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
end
%---------------------------------------------------------------------------------------
return;