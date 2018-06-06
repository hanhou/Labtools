function BATCH_GUI_Tempo_Analysis(batchfiledir, batch_filename, close_flag, print_flag, protocol_filter, analysis_filter)

TEMPO_Defs;
Path_Defs;

ori_filename = batch_filename;
batch_filename = [batchfiledir batch_filename];
fid = fopen(batch_filename);

% Count total number of lines for progressbar. HH20140526
lineNumTotal = 0;
line = fgetl(fid);
while (line~=-1)
    if (line(1) ~= '%')
        lineNumTotal = lineNumTotal + 1;
    end
    line = fgetl(fid);
end

% Reopen the file
fclose(fid);
fid = fopen(batch_filename);
line = fgetl(fid);

% Clear persistent variable "XlsData" in SaveResult.m in case I have changed xls file between calls of BATCH_GUI. HH20150724
clear SaveResult;

% Preparation for speed-up version of "xlswrite1.m"
global Excel;
File='Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm'; 
try
    Excel = actxGetRunningServer('Excel.Application');  % Use the current server
catch
    Excel = actxserver('Excel.Application');  % Start a new server
end

isopen = xls_check_if_open(File,'');
if ~isopen && print_flag ~= -999 % Open file
    winopen('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm');
end

% invoke(Excel.Workbooks,'Open',File); % There is not need for this because "isopen" has ensured that it has been opened.

% Error tolerance is important in batch process! HH20140510
errors = 0;

progressbar('Batch Numer'); lineNum = 0;

while (line ~= -1)
    
    %pause;
    % format for batch files
    % PATH  FILE
    
    if (line(1) ~= '%')
        lineNum = lineNum + 1;
        
        % first remove any comment text at the end of a line (following a %), GCD, added 9/25/01
        comment_start = find(line == '%');
        if ~isempty(comment_start)
            line = line(1:(comment_start(1)-1));
        end
        
        spaces = isspace(line);
        space_index = find(spaces);
        
        %get path / file
        PATH = line(1:space_index(1) - 1);
        FILE = line(space_index(1) + 1:space_index(2) - 1);
        mat_PATH = [PATH(1:end-4) 'Analysis\SortedSpikes2\'];
        
        
        % I use this to indicate that we only need to export the files (for data sharing). HH20141103
        if print_flag == -999
            
            exportPath = [batchfiledir ori_filename(1:strfind(ori_filename,'.')-1) '_Exported\'];
            if ~exist(exportPath,'dir')
                mkdir(exportPath);
            end
            
            htbOK =copyfile([PATH FILE '.htb'],[exportPath FILE '.htb']);
            logOK=copyfile([PATH FILE '.log'],[exportPath FILE '.log']);
            matOK=copyfile([mat_PATH FILE '.mat'],[exportPath FILE '.mat']);
            
            if ~(htbOK && logOK)
                errors = errors + 1;
                errorFiles(errors).fileName = [PATH FILE];
                errorFiles(errors).line = lineNum;
            end
            
            progressbar(lineNum/lineNumTotal);
            line = fgetl(fid);
            continue; % Jump to next file directly
        end
        
        if protocol_filter ~= -1
            
            l = length(FILE);
            if (FILE(l-3:l) == '.htb')	% .htb extension already there
%                 data_filename = [PATH FILE];   %the HTB data file
                logfile = [PATH FILE(1:l-4) '.log'];   %the TEMPO log file
            else	%no extension in FILE, add extensions
%                 data_filename = [PATH FILE '.htb'];   %the HTB data file
                logfile = [PATH FILE '.log'];   %the TEMPO log file
            end
            
            % ------ Read in protocol type and add to cell array ------
            
            % Add protocol names manually to files that have been rescued from CED. HH20150723
            switch logfile 
                case { 'Z:\Data\MOOG\Polo\raw\m5c77r2.log',...
                       'Z:\Data\MOOG\Polo\raw\m5c90r2.log',...
                       'Z:\Data\MOOG\Polo\raw\m5c118r3.log',...
                       'Z:\Data\MOOG\Messi\raw\m10c50r3.log',...
                       'Z:\Data\MOOG\Messi\raw\m10c70r3.log',...
                       'Z:\Data\MOOG\Messi\raw\m10c104r9.log',...
                       'Z:\Data\MOOG\Messi\raw\m10c168r3.log'}   % HD tasks rescued from CED
                    beep;
                    protocol_name = 101;
                case {'Z:\Data\MOOG\Polo\raw\m5c91r1.log'}   % MemSac tasks rescued from CED
                    beep;
                    protocol_name = 125; 
                case {'Z:\Data\MOOG\Messi\raw\m10c102r5.log'} % DelSac task
                    beep;
                    protocol_name = 111;
                otherwise
                    % Read the log file normally
                    protocol_info = textread(logfile, '%s', 3);
                    protocol_name = str2num(protocol_info{2});
            end
            
            if protocol_name == protocol_filter
                %this protocol matches the ones that we want to analyze.
                %proceed with the analysis.
                do_file = 1;
            else %otherwise goto next file.
                do_file = 0;
            end
            
        else
            do_file = 1;
        end % if protocol_filter ~= 1 ...
        
        if do_file == 1
            
            try   % Error tolerance is important in batch process! HH20140510
                
                %get analysis type
                if analysis_filter ~= -1
                    TempAnalysis1 = analysis_strings{protocol_filter + 1}(analysis_filter);
                    TempAnalysis = TempAnalysis1{1};
                    i = 2;
                    while(line(space_index(i)-1) ~= '''')
                        i = i+1;
                    end
                else
                    %special commands for analysis read
                    i = 2;
                    TempAnalysis = '';
                    while(line(space_index(i)-1) ~= '''')
                        i = i+1;
                    end
                    
                    TempAnalysis = strcat(TempAnalysis, line(space_index(2) + 2:space_index(i)-2));
                end
                
                Analysis = {};
                Analysis{1} = TempAnalysis;
                
                %get beginning trial
                BegTrial = line(space_index(i) + 1:space_index(i+1) - 1);
                BegTrial = str2num(BegTrial);
                
                %get ending trial
                EndTrial = line(space_index(i+1) + 1:space_index(i+2) - 1);
                EndTrial = str2num(EndTrial);
                
                %get startcode
                StartCode = line(space_index(i+2) + 1:space_index(i+3) - 1);
                StartCode = str2num(StartCode);
                
                %get startoffset
                StartOffset = line(space_index(i+3) + 1:space_index(i+4) - 1);
                StartOffset = str2num(StartOffset);
                
                %get stopcode
                StopCode = line(space_index(i+4) + 1:space_index(i+5) - 1);
                StopCode = str2num(StopCode);
                
                %get stopoffset
                StopOffset = line(space_index(i+5) + 1:space_index(i+6) - 1);
                StopOffset = str2num(StopOffset);
                
                %get select_data - good trials or bad trials
                good_flag = line(space_index(i+6) + 1:space_index(i+7) - 1);
                good_flag = str2num(good_flag);
                
                %get spikechannels & sync code flags
                
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
                
                %load data file
                [return_val, g_data, b_data] = LoadTEMPOData(PATH,FILE);
                
                if return_val == -999 && SpikeChan >= 5
                    % I throw an error message only when we are in batch mode and we DO use the offline spike2 data. 
                    error('Htb & spike2 files not matched ...');  
                end
                
                if (good_flag)
                    select_data = g_data;
                else
                    select_data = b_data;
                end
                
                if BegTrial == -1
                    BegTrial = 1;
                end
                
                if EndTrial == -1
                    EndTrial = size(g_data.event_data,3);
                end
                
                %get protocol type
                Protocol = g_data.one_time_params(PROTOCOL);		%get the protocol string descriptio
                
                %analyze data
                batch_flag = ori_filename;
                
                if SpikeChan == 999
                    chan_proc = size(select_data.spike_data);
                    if chan_proc > 2
                        SpikeChan = chan_proc(1)
                        
                        Analysis_Switchyard(select_data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag, UseSyncPulses);
                        
                        if print_flag
                            printhandle = figure;
                            close(printhandle);
                            for printindex = 2:1:printhandle - 1
                                print(printindex);
                            end
                        end
                        if close_flag
                            closehandle = figure;
                            for closeindex = 2:1:closehandle
                                close(closeindex);
                            end
                        end
                    else
                        t_string = sprintf('No MU Channels for %s.', FILE);
                        disp(t_string);
                    end
                else
                    
                    if (SpikeChan == 8888)  %tells us to use the last spike channel (MU)
                        %%!!!TEMP KLUGE, GCD 10/8/03 !!!
                        SpikeChan = size(select_data.spike_data,1); %take the last channel for MU data
                        SpikeChan2 = SpikeChan;
                    end
                    
                    Analysis_Switchyard(select_data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag, UseSyncPulses);
                    
                    if print_flag
                        printhandle = figure;
                        close(printhandle);
                        for printindex = 2:1:printhandle - 1
                            print(printindex);
                        end
                    end
                    if close_flag
                        closehandle = figure;
                        for closeindex = 2:1:closehandle
                            close(closeindex);
                        end
                    end
                end
                
            catch exception
                errors = errors + 1;
                errorFiles(errors).fileName = [PATH FILE];
                errorFiles(errors).info = exception;
                errorFiles(errors).line = lineNum;
                
            end
            
        end %do_file
        
        progressbar(lineNum/lineNumTotal);
        
    end %if (line(1) ~= '%')...
    
    line = fgetl(fid);
    
end %while
fclose(fid);

% Close xls server. See "xlswrite1.m" for details
try
    invoke(Excel.ActiveWorkbook,'Save');
    % Excel.Quit
    Excel.delete
    clear global Excel
end

if print_flag == -999
    
    copyfile(batch_filename,[exportPath ori_filename]);
    
else
    outpath = ['Z:\Data\Tempo\Batch\' ori_filename(1:end-2) '\'];
    
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
    
    copyfile(batch_filename,[outpath ori_filename]);
end



home;
if errors > 0
    
    fprintf('\nOops, the following file(s) have error:\n\n');
    for i = 1:errors
        disp(['Batch file Line ' num2str(errorFiles(i).line) ':  ' errorFiles(i).fileName]);
        if isfield(errorFiles(i),'info');
            fprintf('       %s:  line %g , "%s"\n', errorFiles(i).info.stack(1).name,errorFiles(i).info.stack(1).line,errorFiles(i).info.message);
        end
    end
    assignin('base','BatchErrorInfo',errorFiles);
    
    delete([outpath 'BatchError.mat']);
    save([outpath 'BatchError.mat'], 'errorFiles');
else
    fprintf('All success!\n');
end

audio = audioplayer(audioread('Z:\LabTools\Matlab\TEMPO_Analysis\type.wav'),20000);
playblocking(audio);

% Window vibration, HaHa. HH20130829
figure(1001);
set(1001,'Unit','pixel');
currPos = get(1001,'Position');

for ii = 1:6
    set(1001,'Position',[currPos(1) + 10*mod(ii,2) currPos(2) currPos(3) currPos(4)]);
    drawnow;
    pause(0.02);
end
