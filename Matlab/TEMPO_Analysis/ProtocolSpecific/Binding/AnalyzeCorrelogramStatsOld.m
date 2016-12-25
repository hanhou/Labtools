% AnalyzeCorrelogramStats.m - Loads correlogram analysis 
%       11/5/01 - BJP

%clear all; close all; clc;

batchfiledir = 'Z:\Data\Tempo\Batch Files\Binding\'
filename = input('Enter Batch File Name: ', 's');
%close_flag = input('Do you wish to close all windows (0=N, 1=Y)? ');
%print_flag = input('Do you wish to print figures (0=N, 1=Y)? ');

status = ['Loading batch file: ' batchfiledir filename];
disp(status);

filename = [batchfiledir filename];
fid = fopen(filename);
output = 1;
print_fig = 1;
file_num = 1;

line = fgetl(fid);
while (line ~= -1)   
    if (line(1) ~= '%')
	    spaces = isspace(line);
		space_index = find(spaces);

		%get path / file
		PATH = line(1:space_index(1) - 1);
		FILE = line(space_index(1) + 1:space_index(2) - 1);
%        FILE = [FILE(1:end-4) '_s003.ccg']
        FILE = [FILE(1:end-4) '.ccg']
        PATH = PATH(1,1:end-4);
     %   PATH = [PATH 'Analysis\Simulated_Spike_Coherence\'];
        PATH = [PATH 'Analysis\Correlograms\'];
        if (exist([PATH FILE]))  
            eval (['load ' PATH FILE ' -MAT']);       
            num_conditions = length(NCC);
            
            %for cond = 1: num_conditions
                if (print_fig == 1)
                    PATH = PATH(1,1:end-4);
                    %   PATH = [PATH 'Analysis\Simulated_Spike_Coherence\'];
                    PATH = [PATH 'Analysis\Correlograms\'];
                    
                    FILE = [FILE(1:end-4) '.fig']
                    eval (['open ' PATH FILE]);       
                    print;
%                     set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', 'Correlogram Peak Heights');
%                     subplot(3,1,1);
%                     axis([0 100 0 100]);
%                     
%                     axis('off');
%                     xpos = -10;
%                     ypos = 110;
%                     font_size = 9;
%                     bump_size = 10;
%                     
%                     temp = strcat(PATH, FILE);
%                     temp(temp == '\') = '/';
%                     % this prevents a stupid error from appearing on the screen
%                     line = sprintf('File: %s', temp);
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Condition: %2d', cond );
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Window for Averaging of Correlogram Peak: %2d ms', (2*peak_range + 1) );
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Range for Finding Correlogram Peak: %3d to %3d ms', -peak_window, peak_window );
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Number of Bootstraps: %d', num_bootstraps);
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Shuffle-Subtracted Cross Correlogram Peak Height: %4.3f', cross_corr_peak_power(cond) );
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Mean of Bootstrap Correlogram Peak Heights: %4.3f', mean(bootstrap_peak_power(:,cond) ));
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Variance of Bootstrap Correlogram Peak Heights: %4.3f', var(bootstrap_peak_power(:,cond) ));
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Shuffle-Subtracted Cross Correlogram Peak Lag Times: %2.0f ms', cross_corr_peak_time(cond) );
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Mean of Bootstrap Correlogram Peak Lag Times: %4.2f ms', mean(bootstrap_peak_time(:,cond) ));
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('Variance of Bootstrap Correlogram Peak Lag Times: %4.2f', var(bootstrap_peak_time(:,cond) ));
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     line = sprintf('P Value: %7.6f', p_values(cond) );
%                     text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
%                     
%                     histo = NaN*zeros(num_bootstraps,3);
%                     subplot(3,1,2);
%                     histo(:,2) = bootstrap_peak_power(:,cond);
%                     histo(1: fix(num_bootstraps*0.05), 3) = cross_corr_peak_power(cond);
%                     hist(histo,100);
%                     XLabel('Correlogram Peak Height');
%                     YLabel('Count');
%                     
%                     histo = NaN*zeros(num_bootstraps,3);            
%                     subplot(3,1,3);
%                     histo(:,2) = bootstrap_peak_time(:,cond);
%                     histo(1: fix(num_bootstraps*0.2), 3) = cross_corr_peak_time(cond);
%                     hist(histo,(2*peak_window + 1) );
%                     XLabel('Lag Time of Correlogram Peak (ms)');
%                     YLabel('Count');
%                     print;
%                     close(figid);
                end    
                %end    
       
%            power(file_num, :) = power_single_object;
%            file_num = file_num + 1;
            if (output == 1)
                FILEOUT = 'c:\correlogram_stats.txt';
                fileid = [FILEOUT];
                fwriteid = eval(['fopen(fileid, ''a'')']);
                
                fprintf(fwriteid, '%s', FILE);
%                 fprintf(fwriteid, ' %3d', cross_corr_peak_time);
%                 fprintf(fwriteid, ' %6.4f', cross_corr_peak_height);
%                 fprintf(fwriteid, ' %6.4f', cross_corr_peak_area);
%                 fprintf(fwriteid, ' %6.4f', cross_corr_peak_power);
%                 fprintf(fwriteid, ' %6.4f', peak_height_p_values);
%                 fprintf(fwriteid, ' %6.4f', peak_area_p_values);
%                 fprintf(fwriteid, ' %6.4f', power_p_values);
%                fprintf(fwriteid, ' %6.4f', cross_corr_peak_power(1,:)  );
%                fprintf(fwriteid, ' %6.4f', cross_corr_peak_power(2,:)  );
%                fprintf(fwriteid, ' %6.4f', cross_corr_peak_power(3,:)  );
                
                fprintf(fwriteid, ' %6.4f', NC);
                fprintf(fwriteid, ' %6.4f', NCC);
                
%                  fprintf(fwriteid, ' %3d', cross_corr_peak_time);
%                  fprintf(fwriteid, ' %6.4f', cross_corr_peak_height);
%                  fprintf(fwriteid, ' %6.4f', cross_corr_peak_area);
%                  fprintf(fwriteid, ' %6.1f', selected_FR1);
%                  fprintf(fwriteid, ' %6.1f', selected_FR2);
%                  fprintf(fwriteid, ' %6.1f', selected_spikes1);
%                  fprintf(fwriteid, ' %6.1f', selected_spikes2);
%                  fprintf(fwriteid, ' %3d', num_selected_segments);
%                  fprintf(fwriteid, ' %3d', num_reps);
%                  fprintf(fwriteid, ' %3d', minimum_segment_length);
                   
                
                
                
                fprintf(fwriteid, '\r\n');
                                
%                 fprintf(fwriteid, '%s Peak_Time', FILE);
%                 fprintf(fwriteid, ' %3d', cross_corr_peak_time);
%                 fprintf(fwriteid, '\r\n');
%                 
%                 fprintf(fwriteid, '%s Peak_Height', FILE);
%                 fprintf(fwriteid, ' %6.4f', cross_corr_peak_height);
%                 fprintf(fwriteid, '\r\n');
% 
%                 fprintf(fwriteid, '%s Peak_Area', FILE);
%                 fprintf(fwriteid, ' %6.4f', cross_corr_peak_area);
%                 fprintf(fwriteid, '\r\n');
% 
%                 fprintf(fwriteid, '%s Peak_Power', FILE);
%                 fprintf(fwriteid, ' %6.4f', cross_corr_peak_power);
%                 fprintf(fwriteid, '\r\n');
%  
%                 fprintf(fwriteid, '%s Peak_height_p', FILE);
%                 fprintf(fwriteid, ' %6.4f', peak_height_p_values);
%                 fprintf(fwriteid, '\r\n');
% 
%                 fprintf(fwriteid, '%s Peak_area_p', FILE);
%                 fprintf(fwriteid, ' %6.4f', peak_area_p_values);
%                 fprintf(fwriteid, '\r\n');
%                 
%                 fprintf(fwriteid, '%s Power_p', FILE);
%                 fprintf(fwriteid, ' %6.4f', power_p_values);
%                 fprintf(fwriteid, '\r\n');
%                 fprintf(fwriteid, '\r\n');
                
            end    
        else
            sprintf('File %s not exist', FILE);    
        end %if ccg file exists       
    end % if (line(1) ~=...
  	line = fgetl(fid);
end %while...
fclose all;

%figure
%plot(freq2, power);
%xlim([0 50]);