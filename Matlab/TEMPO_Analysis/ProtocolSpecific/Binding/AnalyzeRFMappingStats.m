% AnalyzeRFMappingStats.m - Loads correlogram analysis 
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
        FILE(1,end-2:end) = 'map'
        PATH = PATH(1,1:end-4);
        PATH = [PATH 'Analysis\RF_Mapping\'];
        if (exist([PATH FILE]))  
            eval (['load ' PATH FILE ' -MAT']);       
%            power(file_num, :) = power_single_object;
%            file_num = file_num + 1;
            if (output == 1)
                FILEOUT = 'c:\rf_mapping_stats.txt';
                fileid = [FILEOUT];
                fwriteid = eval(['fopen(fileid, ''a'')']);
                
                fprintf(fwriteid, '%s', FILE);

                fprintf(fwriteid, ' %2d', SpikeChan);
                 fprintf(fwriteid, ' %5.1f', fit_xctr);
                 fprintf(fwriteid, ' %5.1f', fit_yctr);
                 fprintf(fwriteid, ' %5.1f', x_diam);
                 fprintf(fwriteid, ' %5.1f', y_diam);
                 fprintf(fwriteid, ' %5.1f', rf_xctr);
                 fprintf(fwriteid, ' %5.1f', rf_yctr);
                 fprintf(fwriteid, '\r\n');                
            end    
        else
            sprintf('File %s not exist', FILE);    
        end %if dir file exists       
    end % if (line(1) ~=...
  	line = fgetl(fid);
end %while...
fclose all;