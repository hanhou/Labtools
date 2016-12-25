function batch_cmd(batchfiledir, filename)

if nargin < 1 | strcmp(batchfiledir, '');
    batchfiledir = 'Z:\data\Tempo\Batch Files';
end
if nargin < 2 | strcmp(filename, '');
    filename = 'MOOG_vary_fix.m';
end

filename = fullfile(batchfiledir, filename);
fid = fopen(filename);

% copied from BATCH_GUI_Tempo_Analysis
line = fgetl(fid);
cnt = 0;
while (line ~= -1)
    
    %pause
    % format for batch files
    % PATH  FILE    
    
    if (line(1) ~= '%')
        
        % first remove any comment text at the end of a line (following a %), GCD, added 9/25/01
        comment_start = find(line == '%');
        if ~isempty(comment_start)
            line = line(1:(comment_start(1)-1));
        end
        
        spaces = isspace(line);
        space_index = find(spaces);
        
        %get path / file
        PATH = line(1:space_index(1) - 1);
        FILE = line(space_index(1) + 1:space_index(2) - 1)
        l = length(FILE);
        if (FILE(l-3:l) == '.htb')	% .htb extension already there
            filename = [PATH FILE];   %the HTB data file
            logfile = [PATH FILE(1:l-4) '.log'];   %the TEMPO log file
        else	%no extension in FILE, add extensions
            filename = [PATH FILE '.htb'];   %the HTB data file
            logfile = [PATH FILE '.log'];   %the TEMPO log file
        end
        OPT_FILE = fullfile('C:\MATLAB6p5\work\tempo_backdoor', 'Curvefit_default_options.m');

        % go through the backdoor
        cnt = cnt + 1;
        DirectionTuningPlot_Curvefit([],[],[],[],[],[],[],[],[],[],PATH,FILE,OPT_FILE);
        %DirectionTuningPlot_Curvefit([],[],[],[],[],[],[],[],[],[],PATH,FILE);
        %Heading_CurveFit([],[],[],[],[],[],[],[],[],[],PATH,FILE,OPT_FILE);

        pause;
    end
    line = fgetl(fid);
end

fclose(fid);
