%-----------------------------------------------------------------------------------------------------------------------
%-- Microstim.m -- Computes the effect of microstimulation using logistic regression.
%--	TU, 03/05/00
%-----------------------------------------------------------------------------------------------------------------------
function DD_Microstim(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

    Path_Defs;
	TEMPO_Defs;		%defns like IN_T1_WIN_CD
	ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
  
   	Pref_HDisp = data.one_time_params(PREFERRED_HDISP);
   
   	%get the column of values of horiz. disparities in the dots_params matrix
   	h_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
   	unique_hdisp = munique(h_disp');
      
    %get the binocular correlations
    bin_corr = data.dots_params(DOTS_BIN_CORR, :, PATCH1);
    unique_bin_corr = munique(bin_corr');
    signed_bin_corr = bin_corr;
   	signed_bin_corr(h_disp ~= Pref_HDisp) = -bin_corr(h_disp ~= Pref_HDisp);
	unique_signed_bin_corr = munique(signed_bin_corr');

   %get the column of values of microstim in the dots_params matrix
   	mstim = data.misc_params(2,:);
	%mstim = data.misc_params(MICROSTIM,:);
   	unique_mstim = munique(mstim');
   
   	%get indices of any NULL conditions (for measuring spontaneous activity)
   	null_trials = logical( (h_disp == data.one_time_params(NULL_VALUE)) );
   
   	%now, select trials that fall between BegTrial and EndTrial
   	trials = 1:length(h_disp);		% a vector of trial indices
   	select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   
   	%now, determine the choice that was made for each trial, PREFERRED or NULL
   	%by definition, a preferred choice will be made to Target1 and a null choice to Target 2
   	%thus, look for the events IN_T1_WIN_CD and IN_T2_WIN_CD.  GCD, 5/30/2000
   	num_trials = length(h_disp);
   	PREFERRED = 1;
   	NULL = 2;
   	for i=1:num_trials
      	temp = data.event_data(1,:,i);
      	events = temp(temp>0);  % all non-zero entries
      	if (sum(events == IN_T1_WIN_CD) > 0)
         	choice(i) = PREFERRED;
      	elseif (sum(events == IN_T2_WIN_CD) > 0)
         	choice(i) = NULL;
      	else
         	disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    	end        
   	end
         
   	figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 150 500 473], 'Name', 'Microstim');
   	subplot(2, 1, 2);
      
   	% generate data in a form that Stata can process to do logistic regression
    % use this for logit in matlab too 05/30/02 TU
    
	yy=[];
	count = 1;
	for i = 1:length(unique_signed_bin_corr)
		for k = 1:length(unique_mstim)
			yy(count,1) = unique_signed_bin_corr(i);
			yy(count,2) = unique_mstim(k);
			yy(count,3) = sum((choice == PREFERRED) & (signed_bin_corr == unique_signed_bin_corr(i)) & (mstim == unique_mstim(k)) & select_trials);  % # preferred decisions
			yy(count,4) = sum((signed_bin_corr == unique_signed_bin_corr(i)) & (mstim == unique_mstim(k)) & select_trials);		% # trials
			count = count + 1;
		end
	end
   
    %logit in matlab  
   	if (length(unique_mstim) == 1)		% there are no microstim trials, so don't do logit
     		LOGIST_REG = 0;
   	else LOGIST_REG = 1;
    end
     
	if (LOGIST_REG == 1)
        [b, dev, stats] = glmfit([yy(:,1) yy(:,2) yy(:,1).*yy(:,2)],[yy(:,3) yy(:,4)],'binomial');
        
        % compute the best-fitting curves
		bin_corr_ramp = linspace(min(unique_signed_bin_corr), max(unique_signed_bin_corr), 500);
        thresh = []; slope = []; crv = [];
		% first, the no-stim case
		crv(1,:) = 1./(1+exp(-1*(b(1)+b(2).*bin_corr_ramp)));
		thresh(1) = -1*b(1)/b(2);		% 50 pct PD threshold
		slope(1) = b(2);
        se_thresh = (1/b(2))*sqrt(stats.se(1)^2 + (b(1)^2/b(2)^2)*stats.se(2)^2  );
		% now, the stim case
		crv(2,:) = 1./(1+exp(-1*(b(1)+b(3)+b(2).*bin_corr_ramp+b(4).*bin_corr_ramp)));
		thresh(2) = -1*(b(1) + b(3))/(b(2) + b(4));
    	slope(2) = b(2) + b(4);
    	P_thresh = stats.p(3);  % P value for shift
    	P_slope = stats.p(4);	% P value for slope
    else
        disp('no microstim');
        return;
    end

    %---------------------------------------------------------------------------------------
    %Write out threshold se data to a summary file
    buff = sprintf('%s\t %8.6f\t %8.6f\t %8.6f\t %8.6f\t', ...
        FILE, thresh(1), se_thresh, slope(1), 2*se_thresh*slope(1));
    outfile = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\DepthDiscrimMstimThresholdErrs.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t NSThr\t NSTHrSE\t\t NSslope\t\t CIslope\t\t ');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %---------------------------------------------------------------------------------------
    
    
%Old stuff for STATA--------------------------------------------------------------------    
%    % now create the Stata input file
%	cols = 4;
%	file_input = [PATH '..\analysis\stata\' FILE '.stata.inp'];
%   	fid = fopen(file_input,'w');
%   	fprintf(fid, 'disp\tmstim\tch\ttot\n');
%	for i = 1:(count-1)
%		for j = 1:cols
%			fprintf(fid,'%4d\t', yy(i,j));
%		end
%		fprintf(fid,'\n');		
%   	end
%   	fclose(fid);
%      
%    bin_corr_ramp = linspace(min(unique_signed_bin_corr), max(unique_signed_bin_corr), 500);
%   
%   	if (length(unique_mstim) == 1)		% there are no microstim trials, so don't do logit
%     		LOGIST_REG = 0;
%   	else LOGIST_REG = 1;
%    end
%     
%	if (LOGIST_REG == 1)
%    	unix(['copy ' file_input ' Z:\Stata_analysis\temp.dat']);
%		unix('del Z:\Stata_analysis\logistic_reg.log');  %just to make sure
%		unix('D:\Stata\wstata do Z:\Stata_analysis\logistic_reg');
%		[coeffs] = Stata_getLogistRegCoefs('Z:\Stata_analysis\logistic_reg.log');
%		file_output = [PATH '..\analysis\stata\' FILE '.stata.log']
%		unix(['copy Z:\Stata_analysis\logistic_reg.log ' file_output]);
%		unix('del Z:\Stata\logistic_reg.log');
%     
%		% compute the best-fitting curves
%		%thresh = []; slope = []; crv = [];
%		% first, the no-stim case for the first pedestal
%		crv(1,:) = 1./(1+exp(-1*((coeffs(4,1))+coeffs(1,1).*bin_corr_ramp)));
%		thresh(1) = -1*(coeffs(4,1))/coeffs(1,1);		% 50 pct PD threshold
%		slope(1) = coeffs(1,1);
%		% now, the stim case for the first pedestal
%		crv(2,:) = 1./(1+exp(-1*((coeffs(4,1))+coeffs(2,1)+coeffs(1,1).*bin_corr_ramp+coeffs(3,1).*bin_corr_ramp)));
%		thresh(2) = -1*(coeffs(4,1) + coeffs(2,1))/(coeffs(1,1) + coeffs(3,1));
%    	slope(2) = coeffs(1,1) + coeffs(3,1);
%    	P_thresh = coeffs(2,4);  % P value for shift
%    	P_slope = coeffs(3,4);	% P value for slope
%     end
%----------------------------------------------------------------------------------------------------

	% now compute and plot the psychometric functions
   	pdecs = []; y = [];
	for k = 1:length(unique_mstim)
		for i = 1:length(unique_signed_bin_corr)
			y(i,1) = unique_signed_bin_corr(i);
			y(i,2) = sum((choice == PREFERRED) & (signed_bin_corr == unique_signed_bin_corr(i)) & (mstim == unique_mstim(k)) & select_trials) / sum((signed_bin_corr == unique_signed_bin_corr(i)) & (mstim == unique_mstim(k)) & select_trials);  % fraction preferred decisions
			%y(i,3) = sum((ctr_disp == unique_ctr_disp(i)) & (mstim == unique_mstim(k)) & select_trials);		% # trials
         
         pdecs(i,k) = y(i,2);
		end
		%pdecs';
      
       	hold on;
       	if (k == 1) 
   			Handl(1) = plot(unique_signed_bin_corr, pdecs(:,k), 'ko', 'MarkerFaceColor', 'k');
           	%plot([unique_unsigned_hdisp(1) unique_unsigned_hdisp(length(unique_unsigned_hdisp))], [0.5 0.5], 'k-.');
            plot(bin_corr_ramp, crv(k,:), 'k-'); 
       	elseif (k == 2)
    		Handl(2) = plot(unique_signed_bin_corr, pdecs(:,k), 'ko');
           	%plot([unique_unsigned_hdisp(1) unique_unsigned_hdisp(length(unique_unsigned_hdisp))], [0.5 0.5], 'k-.');
            plot(bin_corr_ramp, crv(k,:), 'k--'); 
        end
        hold off;
      
    end	%for k

	xlabel('Binocular Correlation');
	ylabel('Proportion Preferred Decisions');
   	text(-50,1.15,[PATH '\\' FILE]);
   	YLim([0.0 1]);
   	legend(Handl, 'Non-Stim', 'Stim', 2);

  	%now, print out some useful information in the upper subplot
   	subplot(2, 1, 1);
   	PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
   
   	%now, print out some specific useful info.
    xpos = 0; 
    font_size = 9;
    bump_size = 8;

    ypos = 10;
    line = sprintf('Non-stim: threshold = %6.3f deg, slope = %6.3f', thresh(1), slope(1) );
   	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
  	line = sprintf('Stim: threshold = %6.3f deg, slope = %6.3f', thresh(2), slope(2) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('P values: threshold = %6.3f, slope = %6.3f', P_thresh, P_slope );
   	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    xpos = xpos + 50; 
      
%----------------------------------------------------------------------------------------
%now, print out data and fits to a file for external plotting purposes (e.g., in Origin)
char = size(PATH,2) - 1;
while PATH(char) ~='\'	%Analysis directory is one branch below Raw Data Dir
    char = char - 1;
end   
PATHOUT = [PATH(1:char) 'Analysis\Microstim\'];
char = size(FILE,2) - 1;
while FILE(char) ~='.'
    char = char - 1;
end

FILEOUT = [FILE(1:char) 'stim_curves'];
fileid = [PATHOUT FILEOUT];
fwriteid = eval(['fopen(fileid, ''w'')']);
fprintf(fwriteid,'BCorrI\tNstim_fit\tStim_fit\tBCorrR\tNstim_raw\tStim_raw\n');
for i=1:length(bin_corr_ramp)
    fprintf(fwriteid, '%6.3f\t%6.3f\t%6.3f\t', bin_corr_ramp(i), crv(1,i), crv(2,i));
    if (i <= length(unique_signed_bin_corr))
        fprintf(fwriteid, '%6.3f\t%6.3f\t%6.3f\n', unique_signed_bin_corr(i), pdecs(i,1), pdecs(i,2) );
    else
        fprintf(fwriteid, '\n');
    end
end
fclose(fwriteid);
%----------------------------------------------------------------------------------

return;