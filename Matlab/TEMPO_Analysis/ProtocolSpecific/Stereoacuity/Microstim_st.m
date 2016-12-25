%-----------------------------------------------------------------------------------------------------------------------
%-- Microstim_st.m -- Computes the effect of microstimulation using logistic regression.
%--	TU, 11/20/00
%-----------------------------------------------------------------------------------------------------------------------
function Microstim_st(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
TEMPO_Defs;		%defns like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

Pref_HDisp = data.one_time_params(PREFERRED_HDISP);

%get the column of values of horiz. disparities in the dots_params matrix
h_disp_p1 = [];
h_disp_p1 = data.dots_params(DOTS_HDISP,:,PATCH1);
unique_hdisp_p1 = munique(h_disp_p1');

%get the column of values of horiz. disparities in the dots_params matrix
h_disp_p4 = [];
h_disp_p4 = data.dots_params(DOTS_HDISP,:,PATCH4);
unique_hdisp_p4 = munique(h_disp_p4');

ctr_disp = [];
ctr_disp = h_disp_p1 - h_disp_p4;
unique_ctr_disp = munique(ctr_disp'); 

%compute the unsigned horizontal disparity
unsigned_hdisp = abs(ctr_disp);
unsigned_hdisp = (round(10000 * unsigned_hdisp)) / 10000;
unique_unsigned_hdisp = munique(unsigned_hdisp');


%get the column of values of microstim in the dots_params matrix
%something quirky with MATLAB. MICROSTIM CURRENTLY = 2 BJP 3/6/01
mstim = data.misc_params(2,:);
%      mstim = data.misc_params(MICROSTIM,:);
unique_mstim = munique(mstim');

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (h_disp_p1 == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(h_disp_p1);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%now, determine the choice that was made for each trial, PREFERRED or NULL
%by definition, a preferred choice will be made to Target1 and a null choice to Target 2
%thus, look for the events IN_T1_WIN_CD and IN_T2_WIN_CD.  GCD, 5/30/2000
num_trials = length(h_disp_p1);
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

% loop through each pedestal
for ped_idx = 1:length(unique_hdisp_p4)
    
    % generate data in a form that Stata can process to do logistic regression
    % use this for logit in matlab too 05/30/02 TU
    yy=[];
    count = 1;
    for i = 1:length(unique_ctr_disp)
        for k = 1:length(unique_mstim)
            yy(count,1) = unique_ctr_disp(i);
            yy(count,2) = unique_mstim(k);
            yy(count,3) = sum((choice == PREFERRED) & (ctr_disp == unique_ctr_disp(i)) & (h_disp_p4 == unique_hdisp_p4(ped_idx)) & (mstim == unique_mstim(k)) & select_trials);  % # preferred decisions
            yy(count,4) = sum((ctr_disp == unique_ctr_disp(i)) & (h_disp_p4 == unique_hdisp_p4(ped_idx)) & (mstim == unique_mstim(k)) & select_trials);		% # trials
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
        disp_ramp = linspace(min(unique_ctr_disp), max(unique_ctr_disp), 500);
        thresh = []; slope = []; crv = [];
        % first, the no-stim case
        crv(1,:) = 1./(1+exp(-1*(b(1)+b(2).*disp_ramp)));
        thresh(1) = -1*b(1)/b(2);		% 50 pct PD threshold
        slope(1) = b(2);
        se_thresh = (1/b(2))*sqrt(stats.se(1)^2 + (b(1)^2/b(2)^2)*stats.se(2)^2  );
        %[thresh(1) se_thresh slope(1) 2*se_thresh*slope(1)]
        % now, the stim case
        crv(2,:) = 1./(1+exp(-1*(b(1)+b(3)+b(2).*disp_ramp+b(4).*disp_ramp)));
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
    outfile = [BASE_PATH 'ProtocolSpecific\Stereoacuity\StereoacuityMstimThresholdErrs.dat'];
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
    %		% now create the Stata input file
    %		cols = 4;
    %		file_input = [PATH '..\analysis\stata\' FILE '.stata.inp'];
    %   		fid = fopen(file_input,'w');
    %   		fprintf(fid, 'disp\tmstim\tch\ttot\n');
    %		for i = 1:(count-1)
    %			for j = 1:cols
    %				fprintf(fid,'%4d\t', yy(i,j));
    %			end
    %			fprintf(fid,'\n');		
    %   		end
    %   		fclose(fid);
    %   
    %		disp_ramp = linspace(min(unique_ctr_disp), max(unique_ctr_disp), 500);
    %   
    %   		if (length(unique_mstim) == 1)		% there are no microstim trials, so don't do logit
    %      		LOGIST_REG = 0;
    %   		else LOGIST_REG = 1;
    %    	end
    %      
    %		if (LOGIST_REG == 1)
    %      		unix(['copy ' file_input ' Z:\LabTools\Stata_analysis\temp.dat']);
    %			unix('del Z:\LabTools\Stata_analysis\logistic_reg.log');  %just to make sure
    %			unix('D:\Stata\wstata do Z:\LabTools\Stata_analysis\logistic_reg');
    %			%unix('C:\Winapps\Stata\wstata do Z:\LabTools\Stata_analysis\logistic_reg');
    %			[coeffs] = Stata_getLogistRegCoefs('Z:\LabTools\Stata_analysis\logistic_reg.log');
    %			file_output = [PATH '..\analysis\stata\' FILE '.stata.log']
    %			unix(['copy Z:\LabTools\Stata_analysis\logistic_reg.log ' file_output]);
    %			unix('del Z:\LabTools\Stata_analysis\logistic_reg.log');
    %      
    %			% compute the best-fitting curves
    %			%thresh = []; slope = []; crv = [];
    %			% first, the no-stim case for the first pedestal
    %			crv(1,:) = 1./(1+exp(-1*((coeffs(4,1))+coeffs(1,1).*disp_ramp)));
    %			thresh(ped_idx,1) = -1*(coeffs(4,1))/coeffs(1,1);		% 50 pct PD threshold
    %			slope(ped_idx,1) = coeffs(1,1);
    %			% now, the stim case for the first pedestal
    %			crv(2,:) = 1./(1+exp(-1*((coeffs(4,1))+coeffs(2,1)+coeffs(1,1).*disp_ramp+coeffs(3,1).*disp_ramp)));
    %			thresh(ped_idx,2) = -1*(coeffs(4,1) + coeffs(2,1))/(coeffs(1,1) + coeffs(3,1));
    %      		slope(ped_idx,2) = coeffs(1,1) + coeffs(3,1);
    %     		P_thresh(ped_idx) = coeffs(2,4);  % P value for shift
    %      		P_slope(ped_idx) = coeffs(3,4);	% P value for slope
    %        end
    %----------------------------------------------------------------------------------------------------
    
    % now compute and plot the psychometric functions
    pdecs = []; y = [];
    for k = 1:length(unique_mstim)
        for i = 1:length(unique_ctr_disp)
            y(i,1) = unique_ctr_disp(i);
            y(i,2) = sum((choice == PREFERRED) & (ctr_disp == unique_ctr_disp(i)) & (h_disp_p4 == unique_hdisp_p4(ped_idx)) & (mstim == unique_mstim(k)) & select_trials) / sum((ctr_disp == unique_ctr_disp(i)) & (h_disp_p4 == unique_hdisp_p4(ped_idx)) & (mstim == unique_mstim(k)) & select_trials);  % fraction preferred decisions
            %y(i,3) = sum((ctr_disp == unique_ctr_disp(i)) & (mstim == unique_mstim(k)) & select_trials);		% # trials
            
            pdecs(i,k) = y(i,2);
        end
        %pdecs';
        
        hold on;
        if (ped_idx == 1)
            if (k == 1) 
                Handl(1) = plot(unique_ctr_disp, pdecs(:,k), 'ko', 'MarkerFaceColor', 'k');
                %plot([unique_unsigned_hdisp(1) unique_unsigned_hdisp(length(unique_unsigned_hdisp))], [0.5 0.5], 'k-.');
                plot(disp_ramp, crv(k,:), 'k-'); 
            elseif (k == 2)
                Handl(2) = plot(unique_ctr_disp, pdecs(:,k), 'ko');
                %plot([unique_unsigned_hdisp(1) unique_unsigned_hdisp(length(unique_unsigned_hdisp))], [0.5 0.5], 'k-.');
                plot(disp_ramp, crv(k,:), 'k--'); 
            end
        elseif (ped_idx == 2)
            if (k == 1) 
                Handl(1) = plot(unique_ctr_disp, pdecs(:,k), 'kv', 'MarkerFaceColor', 'k');
                %plot([unique_unsigned_hdisp(1) unique_unsigned_hdisp(length(unique_unsigned_hdisp))], [0.5 0.5], 'k-.');
                plot(disp_ramp, crv(k,:), 'k-'); 
            elseif (k == 2)
                Handl(2) = plot(unique_ctr_disp, pdecs(:,k), 'kv');
                %plot([unique_unsigned_hdisp(1) unique_unsigned_hdisp(length(unique_unsigned_hdisp))], [0.5 0.5], 'k-.');
                plot(disp_ramp, crv(k,:), 'k--'); 
            end                  
        end      
        hold off;
        
    end	%for k
end %for ped     

xlabel('Binocular Disparity (deg)');
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
for ped_idx = 1:length(unique_hdisp_p4)
    ypos = 10;
    line = sprintf('Non-stim: threshold = %6.3f deg, slope = %6.3f', thresh(ped_idx,1), slope(ped_idx,1) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Stim: threshold = %6.3f deg, slope = %6.3f', thresh(ped_idx,2), slope(ped_idx,2) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('P values: threshold = %6.3f, slope = %6.3f', P_thresh(ped_idx), P_slope(ped_idx) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Pedestal disparity: %6.3f deg', unique_hdisp_p4(ped_idx) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    xpos = xpos + 50; 
end


%get the average eye positions to calculate vergence
Leyex_positions = data.eye_positions(1, :);
Leyey_positions = data.eye_positions(2, :);
Reyex_positions = data.eye_positions(3, :);
Reyey_positions = data.eye_positions(4, :);

vergence_h = Leyex_positions - Reyex_positions;
vergence_v = Leyey_positions - Reyey_positions;

if (data.eye_calib_done == 1)
    Leyex_positions = data.eye_positions_calibrated(1, :);
    Leyey_positions = data.eye_positions_calibrated(2, :);
    Reyex_positions = data.eye_positions_calibrated(3, :);
    Reyey_positions = data.eye_positions_calibrated(4, :);
    
    vergence_h = Leyex_positions - Reyex_positions;
    vergence_v = Leyey_positions - Reyey_positions;
end

% select1 = logical(mstim == 0);
% select2 = logical(mstim == 1);
% figure;
% plot(ctr_disp(select1), vergence_h(select1), 'ro');
% hold on;
% plot(ctr_disp(select2), vergence_h(select2), 'bo');

[H, ATAB, CTAB, STATS] = aoctool(ctr_disp, vergence_h, mstim);

P_mstim = ATAB{2,6}
P_ctr_disp = ATAB{3,6}
P_interact = ATAB{4,6}
std_verg = std(vergence_h)
mean_verg = mean(vergence_h)
verg_ctr_disp_slope = CTAB{5,2}
num_pts = length(vergence_h)

close(4); close(5);

%---------------------------------------------------------------------------------------
%Write out some summary data to a cumulative summary file
buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %8.6f\t %8.6f\t %8.6f\t %7.3f\t %7.3f\t %7.3f\t %6d\t %6.3f\t', ...
    FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
    P_mstim, P_ctr_disp, P_interact, std_verg, mean_verg, verg_ctr_disp_slope, num_pts, unique_hdisp_p4(1));
outfile = [BASE_PATH 'ProtocolSpecific\Stereoacuity\StereoacuityMstimVergence_analysis.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Pmstim\t\t Pctrdsp\t Pinter\t\t SDverg\t\t AvVerg\t\t Slope\t\t Npts\t PedDsp\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------

%--------------------------------------------------------------------------
%now, print out data and fits to a file for external plotting purposes (e.g., in Origin)
i = size(PATH,2) - 1;
while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
    i = i - 1;
end   
PATHOUT = [PATH(1:i) 'Analysis\Microstim\'];
i = size(FILE,2) - 1;
while FILE(i) ~='.'
    i = i - 1;
end
FILEOUT = [FILE(1:i) 'curve_fits'];

fileid = [PATHOUT FILEOUT]
fwriteid = eval(['fopen(fileid, ''w'')']);

fprintf(fwriteid,'DispI\tS_fit\tNS_fit\tDispR\tS_raw\tNS_raw\n');
for i=1:length(disp_ramp)
    fprintf(fwriteid, '%8.5f\t%6.3f\t%6.3f\t', disp_ramp(i), crv(2,i), crv(1,i));
    if (i <= length(unique_ctr_disp))
        fprintf(fwriteid, '%8.5f\t%6.3f\t%6.3f\n', unique_ctr_disp(i),  pdecs(i,2), pdecs(i,1));
    else
        fprintf(fwriteid, '\n');
    end
end

fclose(fwriteid);


return;