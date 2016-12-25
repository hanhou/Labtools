% AnalyzeCorrelogramStats.m - Loads correlogram analysis 
%       11/5/01 - BJP

%clear all; close all; clc;

ProtocolDefs

batchfiledir = 'z:\Data\Tempo\Batch Files\Binding\'
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
        FILE(1,end-2:end) = 'log'
        if (exist([PATH FILE]))  
            clear dots_params; 
            clear targ_params;
            clear misc_params; 
            clear one_time_params;
            clear eye_calib_params; 
            clear obj_params;
            clear bar_params;
            clear bkgnd_params;
            clear neuron_params;
            [dots_params, moog_params, targ_params, misc_params, one_time_params, eye_calib_params, obj_params, bar_params, bkgnd_params, neuron_params] = ReadTEMPOLog([PATH FILE], 300);
            
            if (output == 1)
                FILEOUT = 'c:\tempo_log_stats.txt';;
                fileid = [FILEOUT];
                fwriteid = eval(['fopen(fileid, ''a'')']);
                
                fprintf(fwriteid, '%s', FILE);
                fprintf(fwriteid, ' %6.2f ', neuron_params(RF_XCTR,1) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(RF_YCTR,1) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(RF_DIAMETER,1) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(PREFERRED_DIRECTION,1) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(PREFERRED_SPEED,1) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(RF_XCTR,2) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(RF_YCTR,2) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(RF_DIAMETER,2) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(PREFERRED_DIRECTION,2) );
                fprintf(fwriteid, ' %6.2f ', neuron_params(PREFERRED_SPEED,2) );
                
           
                fprintf(fwriteid, '\r\n');
                                                
            end    
        else
            sprintf('File %s not exist', FILE);    
        end %if ccg file exists       
    end % if (line(1) ~=...
  	line = fgetl(fid);
end %while...
fclose all;