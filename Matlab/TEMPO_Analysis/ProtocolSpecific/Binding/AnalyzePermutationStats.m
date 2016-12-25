% AnalyzeCorrelogramStats.m - Loads correlogram analysis 
%       11/5/01 - BJP

clear all; close all; clc;

batchfiledir = 'Z:\Data\Tempo\Batch Files\'
filename = input('Enter Batch File Name: ', 's');
%close_flag = input('Do you wish to close all windows (0=N, 1=Y)? ');
%print_flag = input('Do you wish to print figures (0=N, 1=Y)? ');

status = ['Loading batch file: ' batchfiledir filename];
disp(status);

filename = [batchfiledir filename];
fid = fopen(filename);
output = 1;
print_fig = 0;

line = fgetl(fid);
while (line ~= -1)   
    if (line(1) ~= '%')
	    spaces = isspace(line);
		space_index = find(spaces);

		%get path / file
		PATH = line(1:space_index(1) - 1);
		FILE = line(space_index(1) + 1:space_index(2) - 1);
        FILE(1,end-2:end) = 'prm'
        PATH = PATH(1,1:end-4);
        PATH = [PATH 'Analysis\Correlograms\'];
        if (exist([PATH FILE])) 
            sync_p_values = [];
            comp_p_values = [];
            eval (['load ' PATH FILE ' -MAT']);       
            num_conditions = length(cross_corr_peak_power);
            for cond = 1: num_conditions
                if (print_fig == 1)
                    figid = figure
                    set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', 'Correlogram Peak Heights');
                    subplot(3,1,1);
                    axis([0 100 0 100]);
                    
                    axis('off');
                    xpos = -10;
                    ypos = 110;
                    font_size = 9;
                    bump_size = 10;
                    
                    temp = strcat(PATH, FILE);
                    temp(temp == '\') = '/';
                    % this prevents a stupid error from appearing on the screen
                    line = sprintf('File: %s', temp);
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Condition: %2d', cond );
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Window for Averaging of Correlogram Peak: %2d ms', (2*peak_range + 1) );
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Range for Finding Correlogram Peak: %3d to %3d ms', -peak_window, peak_window );
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Number of Bootstraps: %d', num_bootstraps);
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Shuffle-Subtracted Cross Correlogram Peak Power: %4.3f', cross_corr_peak_power(cond) );
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Mean of Bootstrap Correlogram Peak Power: %4.3f', mean(bootstrap_peak_power(:,cond) ));
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Variance of Bootstrap Correlogram Peak Heights: %4.3f', var(bootstrap_peak_heights(:,cond) ));
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Shuffle-Subtracted Cross Correlogram Peak Lag Times: %2.0f ms', cross_corr_peak_time(cond) );
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Mean of Bootstrap Correlogram Peak Lag Times: %4.2f ms', mean(bootstrap_peak_time(:,cond) ));
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('Variance of Bootstrap Correlogram Peak Lag Times: %4.2f', var(bootstrap_peak_time(:,cond) ));
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    line = sprintf('P Value: %7.6f', power_p_values(cond) );
                    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
                    
                    histo = NaN*zeros(num_bootstraps,3);
                    subplot(3,1,2);
                    histo(:,2) = bootstrap_peak_heights(:,cond);
                    histo(1: fix(num_bootstraps*0.05), 3) = cross_corr_peak_power(cond);
                    hist(histo,100);
                    XLabel('Correlogram Peak Height');
                    YLabel('Count');
                    
                    histo = NaN*zeros(num_bootstraps,3);            
                    subplot(3,1,3);
                    histo(:,2) = bootstrap_peak_time(:,cond);
                    histo(1: fix(num_bootstraps*0.2), 3) = cross_corr_peak_time(cond);
                    hist(histo,(2*peak_window + 1) );
                    XLabel('Lag Time of Correlogram Peak (ms)');
                    YLabel('Count');
                    print;
                    close(figid);
                end    
            end    
            
            if (output == 1)
                    for comp = 1: length(corr_peak_diff)
                        comp_p_values(comp) = mean( perm_peak_diff(:,comp) > corr_peak_diff(comp) );
                    end
                PATH = PATH(1,1:end-4);
                PATH = [PATH 'Analysis\Correlograms\'];

                FILEOUT = 'perm_values.txt';
                fileid = [PATH FILEOUT];
                fwriteid = eval(['fopen(fileid, ''a'')']);
                
                fprintf(fwriteid, '%s', FILE);
                fprintf(fwriteid, ' %4.3f', sync_p_values  );
                fprintf(fwriteid, '\r\n');
                
                fprintf(fwriteid, '%s', FILE);
                fprintf(fwriteid, ' %0.6f', comp_p_values);
                fprintf(fwriteid, '\r\n');
                fprintf(fwriteid, '\r\n');
                
            end    
        end %if ccg file exists       
    end % if (line(1) ~=...
  	line = fgetl(fid);
end %while...
fclose all;

