function [return_value, info_text] = GetBatchInfo(BatchPATH, FILEName)

TEMPO_Defs;
ProtocolDefs;
Path_Defs;

BatchFile = [BatchPATH FILEName];
fid = fopen(BatchFile);
line = fgetl(fid);

%get the filename to be processed, then determine what protocol it is
%then repeat for each line
try

info_text = [];
while (line ~= -1)
    if (line(1) ~= '%')
        
        spaces = isspace(line);
        space_index = find(spaces);
        %get path / file
        PATH = line(1:space_index(1) - 1);
        FILE = line(space_index(1) + 1:space_index(2) - 1);
        
        l = length(FILE);
        if (FILE(l-3:l) == '.htb')	% .htb extension already there
            filename = [PATH FILE];   %the HTB data file
            logfile = [PATH FILE(1:l-4) '.log'];   %the TEMPO log file
        else	%no extension in FILE, add extensions
            filename = [PATH FILE '.htb'];   %the HTB data file
            logfile = [PATH FILE '.log'];   %the TEMPO log file
        end
        
        %read in analysis type
        i = 3;
        done = 0;
        while ~done
            endanalysis = line(space_index(i) - 1);
            if endanalysis == ''''
                done = 1;
            end
            i = i + 1;
        end
        i = i-1;
        analysis = line(space_index(2)+1:space_index(i)-1);
        
        % Add some more to check the batch file (From BATCH_GUI_Tempo_Analysis). HH20150413
        %get beginning trial
        BegTrial = line(space_index(i) + 1:space_index(i+1) - 1);
        
        %get ending trial
        EndTrial = line(space_index(i+1) + 1:space_index(i+2) - 1);
        
        if (i + 9 - 1 < length(space_index) )
            %if there are sync pulse value
            SpikeChan = line(space_index(i+7) + 1:space_index(i+8) - 1);
            SpikeChan = str2num(SpikeChan);
            SpikeChan2 = line(space_index(i+8) + 1:space_index(i+9) - 1);
            SpikeChan2 = str2num(SpikeChan2);
            UseSyncPulses = line(space_index(i+9) +  1:length(line) );
            UseSyncPulses = str2num(UseSyncPulses);
        else
            SpikeChan = line(space_index(i+7) + 1:length(line));
            SpikeChan = str2num(SpikeChan);
            SpikeChan2 = SpikeChan;
            %if sync pulses flag not set, assume we want to use TEMPO START_CD & STOP_CD
            UseSyncPulses = 0;
        end
        
        
        %read in protocol type and add to cell array
        if exist(logfile)
            protocol_info = textread(logfile, '%s', 3);
            protocol_name = protocol_info{3};
            protocol_name = protocol_name(2:length(protocol_name)-1);
        else
            protocol_name = '** Rescued From CED **   ';
        end
        
        %t = sprintf('%s \t %s \t %s \t %s',PATH, FILE, protocol_name, analysis);
        t = sprintf('%-25s %-15s %-20s %-32s %4g',PATH, FILE, protocol_name, analysis, SpikeChan);
        
        info_text{length(info_text)+1} = t;
    end
    line = fgetl(fid);
end

fclose(fid);
return_value = length(info_text);

catch err
    err
    line
    keyboard;
    fclose(fid);
end
return


