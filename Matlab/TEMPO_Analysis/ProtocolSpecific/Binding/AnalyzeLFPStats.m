% AnalyzeLFPStats.m - Loads LFP analyses
%       11/5/01 - BJP

%clear all; close all; clc;

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
file_num = 1;

line = fgetl(fid);
while (line ~= -1)   
    if (line(1) ~= '%')
	    spaces = isspace(line);
		space_index = find(spaces);

		%get path / file
		PATH = line(1:space_index(1) - 1);
		FILE = line(space_index(1) + 1:space_index(2) - 1);
        FILE(1,end-2:end) = 'lfp'
        PATH = PATH(1,1:end-4);
        PATH = [PATH 'Analysis\LFP\'];
        if (exist([PATH FILE]))  
            eval (['load ' PATH FILE ' -MAT']);       
            num_conditions = size(low_band, 1);
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
                end    
            end    
       
            if (output == 1)
                FILEOUT = 'c:\LFP_stats.txt';
                fileid = [FILEOUT];
                fwriteid = eval(['fopen(fileid, ''a'')']);
                
                fprintf(fwriteid, '%s', FILE);
                fprintf(fwriteid, ' %8.4f', rms(1, :)  );
                fprintf(fwriteid, ' %8.4f', rms(2, :)  );

                fprintf(fwriteid, ' %8.4f', low_band(:, 1, 1)  );
                fprintf(fwriteid, ' %8.4f', low_band(:, 1, 2)  );
          
                fprintf(fwriteid, ' %8.4f', band1(:, 1, 1)  );
                fprintf(fwriteid, ' %8.4f', band1(:, 1, 2)  );
                
                fprintf(fwriteid, ' %8.4f', band2(:, 1, 1)  );
                fprintf(fwriteid, ' %8.4f', band2(:, 1, 2)  );

 %               fprintf(fwriteid, ' %8.4f', band3(:, 1, 1)  );
 %               fprintf(fwriteid, ' %8.4f', band3(:, 1, 2)  );

  %              fprintf(fwriteid, ' %8.4f', high_band(:, 1, 1)  );
  %              fprintf(fwriteid, ' %8.4f', high_band(:, 1, 2)  );

                fprintf(fwriteid, '\r\n');
                                

                
            end    
        else
            sprintf('File %s not exist', FILE);    
        end %if ccg file exists       
    end % if (line(1) ~=...
  	line = fgetl(fid);
end %while...
fclose all;
