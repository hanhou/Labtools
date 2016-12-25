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
print_fig = 0;
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
            
            if (print_fig == 1)
                PATH = PATH(1,1:end-4);
                %   PATH = [PATH 'Analysis\Simulated_Spike_Coherence\'];
                PATH = [PATH 'Analysis\Correlograms\'];
                
                FILE = [FILE(1:end-4) '.fig']
                eval (['open ' PATH FILE]);       
                print;
            end    
            
            %            power(file_num, :) = power_single_object;
%            file_num = file_num + 1;
            if (output == 1)
                FILEOUT = 'c:\correlogram_stats.txt';
                fileid = [FILEOUT];
                fwriteid = eval(['fopen(fileid, ''a'')']);
                
                %baseline of fit for on bar condition = pars{2}(1)
               
                               
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
                
%                fprintf(fwriteid, ' %6.4f', NCC);
%                 fprintf(fwriteid, ' %6.4f', NCC2);
%                 fprintf(fwriteid, ' %6.4f', NCC3);
%                 fprintf(fwriteid, ' %6.4f', NCC4);
%                 fprintf(fwriteid, ' %6.4f', NCC5);
                 fprintf(fwriteid, ' %6.4f', length( find(fit_cross_correlogram(2,:) >  0.5*(max(fit_cross_correlogram(2,:) ) -  pars{2}(1) ) + pars{2}(1)  ) ) );


                 fprintf(fwriteid, ' %6.4f', permutation_NCC_p_values(2) );
%                 fprintf(fwriteid, ' %6.4f', NCC_diff_p_value );
                
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