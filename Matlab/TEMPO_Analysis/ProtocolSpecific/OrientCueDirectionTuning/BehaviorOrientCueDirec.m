%-----------------------------------------------------------------------------------------------------------------------
%-- BehaviorOrientCueDirec.m -- Plots polar performance curves
%--	VR, 6/28/05
%-----------------------------------------------------------------------------------------------------------------------
function BehaviorOrientCueDirec(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffset, StartEvent, StopEvent, PATH, FILE);

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);
unique_direction = munique(direction');

%get the motion coherences
coherence = data.dots_params(DOTS_COHER, :, PATCH1);
unique_coherence = munique(coherence');

Pref_direction = data.one_time_params(PREFERRED_DIRECTION);
Cue_direction = data.cue_params(CUE_DIREC);

%get the cue validity
cue_val = data.cue_params(CUE_VALIDITY,:,PATCH2);
cue_direc = data.cue_params(CUE_DIREC, :, PATCH1);
unique_cue_val = munique(cue_val');

%compute cue status - 0=off, 1=on, given 0=off, -1=invalid, 1=valid
cue_status = cue_val.^2; %maps neutral 0->0, valid/invalid [1 -1] -> 1 
unique_cue_status = munique(cue_status');

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (coherence == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(coherence);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%[direction' coherence' spike_rates' null_trials' select_trials']



% figure;
% set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Psychometric Function');
% subplot(2, 1, 2);

symbols = {'bo', 'rx', 'g*'};
line_colors = {'b-', 'r-', 'g:'};
names = {'Neutral Cue', 'Orientation Cue'};
%% *********** POLAR PLOTS ******************************
% Because mmpolar.m seems to hang up when replotting with hold on, the
% performance data variables are appended to a string and the command is
% evaluated at once.

unique_theta = unique_direction .* pi ./ 180; %convert directions to radians

for i = 1:length(unique_coherence)
   figure;  set(gcf, 'Name', sprintf('%s: %d%% coherence',FILE, unique_coherence(i)));
   polar_str = 'mmpolar(';
   legend_str = 'legend(';
   for j = 1:length(unique_cue_status)
      for k = 1:length(unique_direction)
          t(k,j,i) = length(find( (coherence == unique_coherence(i)) & (cue_status == unique_cue_status(j)) ...
                                & (direction == unique_direction(k)) ));
          c(k,j,i) = length(find( (coherence == unique_coherence(i)) & (cue_status == unique_cue_status(j)) ...
                                & (direction == unique_direction(k)) & (data.misc_params(OUTCOME, :) == CORRECT) ));
          if (t(k,j,i) > 0)
              p(k,j,i) = c(k,j,i)./t(k,j,i);
          else
              p(k,j,i) = 0.5;
          end
      end
      %plot on polar axes.  append first element to end to complete the loop
      polar_str = strcat(polar_str,'[unique_theta; unique_theta(1)], [p(:,',sprintf('%d',j),...
                         ',i); p(1,', sprintf('%d',j),...
                         ',i)], line_colors{', sprintf('%d',j),'},');
      legend_str = strcat(legend_str, sprintf(' ''%s'',', names{unique_cue_status(j)+1}));
      %keyboard
      hold on;
   end
   polar_str = strcat(polar_str,' ''RLimit'', [0 1], ''TTickValue'', squeeze_angle(unique_direction));');
   eval(polar_str); 
   legend_str = strcat(legend_str,' ''Location'', ''BestOutside'');');
   eval(legend_str);
   hold on
   title(sprintf('Coherence = %3.1f%%',unique_coherence(i)));
   plot ([cos(Cue_direction*pi/180) cos((Cue_direction+180)/180*pi)], [sin(Cue_direction*pi/180) sin((Cue_direction+180)/180*pi)],'c-')
   plot ([cos(Pref_direction*pi/180) cos((Pref_direction+180)/180*pi)], [sin(Pref_direction*pi/180) sin((Pref_direction+180)/180*pi)],'sk','MarkerFaceColor','k')
end


% %% *********** PSYCHOMETRIC ANALYSIS ****************************
% 
% pct_correct = []; N_obs = []; fit_data = [];
% monkey_alpha = []; monkey_beta = [];
% 
% for j = 1:length(unique_cue_val)
%     for i=1:length(unique_coherence)
%         trials = ((coherence == unique_coherence(i)) & (cue_val == unique_cue_val(j))& select_trials);
%         correct_trials = (trials & (data.misc_params(OUTCOME, :) == CORRECT) );
%         pct_correct(i) = sum(correct_trials)/sum(trials);
%         N_obs(i) = sum(trials);
%         % data for Weibull fit
%         fit_data(i, 1) = unique_coherence(i);
%         fit_data(i, 2) = pct_correct(i);
%         fit_data(i, 3) = N_obs(i);
%     end
%     hold on;
%     unique_coherence_hack = unique_coherence;
%     unique_coherence_hack(1) = 1.01;
%     Handl(1) = plot(unique_coherence_hack, pct_correct, symbols{j});
%     hold off;
%     
%     fit_x = unique_coherence(1):0.1: unique_coherence(length(unique_coherence));
%     [monkey_alpha(j) monkey_beta(j)]= weibull_fit(fit_data);
%     monkey_fit_y = weibull_curve(fit_x, [monkey_alpha(j) monkey_beta(j)]);
%  
%     hold on;
%     plot(fit_x, monkey_fit_y, lines{j});
%     hold off;
% end
% 
% xlabel('Coherence (% dots)');
% ylabel('Fraction Correct');
% 
% YLim([0.4 1]);
% %comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
% set(gca, 'XScale', 'log');
% XLim([1 100]);
% 
% %now, print out some useful information in the upper subplot
% subplot(2, 1, 1);
% PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
% 
% start_time = find(data.event_data(1, :, 1) == VSTIM_ON_CD);
% stop_time = find(data.event_data(1, :, 1) == VSTIM_OFF_CD);
% stim_duration = stop_time - start_time
% 
% %now, print out some specific useful info.
% xpos = 0; ypos = 10;
% font_size = 11;
% bump_size = 8;
% cuevals = {'Invalid','Neutral','Valid'};
% for j = 1:length(unique_cue_val)
%     line = sprintf('Monkey: CueStatus = %s, threshold = %6.3f %%, slope = %6.3f', cuevals{unique_cue_val(j)+2}, monkey_alpha(j), monkey_beta(j) );
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% end
% line = sprintf('Directions tested: %6.3f, %6.3f deg', unique_direction(1), unique_direction(2) );
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% line = sprintf('Stimulus Duration: %5d', stim_duration );
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% %line = sprintf('(3,%0.5g) (6,%0.5g) (12,%0.5g) (24,%0.5g) (48,%0.5g) (96,%0.5g)', ...
% %    pct_correct(1), pct_correct(2), pct_correct(3), pct_correct(4), pct_correct(5), pct_correct(6) );
% %text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% 
% output = 0;
% if (output)
%     
%     %------------------------------------------------------------------------
%     %write out all relevant parameters to a cumulative text file, GCD 8/08/01
%     outfile = [BASE_PATH 'ProtocolSpecific\DirecDiscrim\Psycho_Curve_summary.dat'];
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t MThr\t MSlp\t DspLo\t DspHi\t Ntrials\t HCorr\t Durat\t');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
%     buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.4f\t %6.3f\t %6.3f\t %4d\t %6.3f\t %5d\t', ...
%         FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%         monkey_alpha,monkey_beta,unique_direction(1), unique_direction(2), (1+ EndTrial - BegTrial), unique_coherence(length(unique_coherence)), stim_duration );
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
%     %------------------------------------------------------------------------
% end

return;