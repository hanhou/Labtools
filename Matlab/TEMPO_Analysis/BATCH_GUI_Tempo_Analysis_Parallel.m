function BATCH_GUI_Tempo_Analysis_Parallel(batchfiledir, batch_filename, close_flag, print_flag, protocol_filter, analysis_filter)

hostname = char( getHostName( java.net.InetAddress.getLocalHost)); % Get host name
if contains(hostname,'clc') || contains(hostname,'node')
    ifCluster = 1;
else
    ifCluster = 0;
end
 
% For parallel processing
if ~verLessThan('matlab','R2014a')
    if isempty(gcp('nocreate'))
        parpool('local',17,'IdleTimeout',inf);
    end
else
    if matlabpool('size') == 0
        try
            matlabpool;
        catch
        end
    end
end

TEMPO_Defs;
Path_Defs;

ori_filename = batch_filename;
batch_filename = [batchfiledir batch_filename];

% fid = fopen(batch_filename);
% Count total number of lines for progressbar. HH20140526
% lineNumTotal = 0;
% line = fgetl(fid);
% while (line~=-1)
%     if (line(1) ~= '%')
%         lineNumTotal = lineNumTotal + 1;
%     end
%     line = fgetl(fid);
% end
% fclose(fid);

% Reopen the file
fid = fopen(batch_filename);
line = fgetl(fid);

% Clear persistent variable "XlsData" in SaveResult.m in case I have changed xls file between calls of BATCH_GUI. HH20150724
% clear SaveResult;
%
% % Preparation for speed-up version of "xlswrite1.m"
% global Excel;
% File='Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm';
% try
%     Excel = actxGetRunningServer('Excel.Application');  % Use the current server
% catch
%     Excel = actxserver('Excel.Application');  % Start a new server
% % end
%
% isopen = xls_check_if_open(File,'');
% if ~isopen && print_flag ~= -999 % Open file
%     winopen('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm');
% end

% invoke(Excel.Workbooks,'Open',File); % There is not need for this because "isopen" has ensured that it has been opened.

% Error tolerance is important in batch process! HH20140510
errors = 0;

% ================ First get all lines, then do parallel computing on the cluster. HH20180607 ============
nn = 0; % Use to count files that needed to be batched.
tic

while (line ~= -1) % Propocessing
    
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
        FILE = line(space_index(1) + 1:space_index(2) - 1);
        
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
            if ifCluster
                switch logfile
                    case { 'Z/Data/MOOG/Polo/raw/m5c77r2.log',...
                            'Z/Data/MOOG/Polo/raw/m5c90r2.log',...
                            'Z/Data/MOOG/Polo/raw/m5c118r3.log',...
                            'Z/Data/MOOG/Messi/raw/m10c50r3.log',...
                            'Z/Data/MOOG/Messi/raw/m10c70r3.log',...
                            'Z/Data/MOOG/Messi/raw/m10c104r9.log',...
                            'Z/Data/MOOG/Messi/raw/m10c168r3.log'}   % HD tasks rescued from CED
                        beep;
                        protocol_name = 101;
                    case {'Z/Data/MOOG/Polo/raw/m5c91r1.log'}   % MemSac tasks rescued from CED
                        beep;
                        protocol_name = 125;
                    case {'Z/Data/MOOG/Messi/raw/m10c102r5.log'} % DelSac task
                        beep;
                        protocol_name = 111;
                    otherwise
                        % Read the log file normally
                        protocol_info = textscan(logfile, '%s', 3);
                        protocol_name = str2double(protocol_info{2});
                        
                end
            else
                switch logfile
                    case {  'Z:\Data\MOOG\Polo\raw\m5c77r2.log',...
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
                        protocol_info = textscan(logfile, '%s', 3);
                        protocol_name = str2double(protocol_info{2});
                        
                end
            end
            
            do_file = protocol_name == protocol_filter;
            %             if protocol_name == protocol_filter
            %                 %this protocol matches the ones that we want to analyze.
            %                 %proceed with the analysis.
            %                 do_file = 1;
            %             else %otherwise goto next file.
            %                 do_file = 0;
            %             end
            
        else  % if no protocol_filter, do all
            do_file = 1;
        end 

        
        if do_file == 1
            
            
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
            BegTrial = str2num(BegTrial);  % Must be str2num!!!
            
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
            
            % ==== Cache files ==== % 
            
            nn = nn + 1;
            para_batch{nn} = {Analysis, SpikeChan , SpikeChan2, StartCode, StopCode,...
                BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, UseSyncPulses, good_flag};
        
        end
        
    end
    
    line = fgetl(fid);
end
% pause;
% format for batch files
% PATH  FILE


%                 if SpikeChan == 999
%                     chan_proc = size(select_data.spike_data);
%                     if chan_proc > 2
%                         SpikeChan = chan_proc(1)
%
%                         Analysis_Switchyard(select_data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag, UseSyncPulses);
%
%                         if print_flag
%                             printhandle = figure;
%                             close(printhandle);
%                             for printindex = 2:1:printhandle - 1
%                                 print(printindex);
%                             end
%                         end
%                         if close_flag
%                             closehandle = figure;
%                             for closeindex = 2:1:closehandle
%                                 close(closeindex);
%                             end
%                         end
%                     else
%                         t_string = sprintf('No MU Channels for %s.', FILE);
%                         disp(t_string);
%                     end
%                 else

%                     if (SpikeChan == 8888)  %tells us to use the last spike channel (MU)
%                         %%!!!TEMP KLUGE, GCD 10/8/03 !!!
%                         SpikeChan = size(select_data.spike_data,1); %take the last channel for MU data
%                         SpikeChan2 = SpikeChan;
%                     end

N = nn; % Total number of files
errorFiles(N).fileName = [];

% ======================== Do parallel processing ===================

parfor_progress(N);

PROTOCOLForParWorkers = PROTOCOL;

parfor nn = 1:N
    try   % Error tolerance is important in batch process! HH20140510
        
        %                       1         2          3           4         5
        % para_batch{nn} = {Analysis, SpikeChan, SpikeChan2, StartCode, StopCode,...
        %                  6         7           8          9         10    11        12            13
        %               BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, UseSyncPulses, good_flag};
        
        %load data file
        [return_val, g_data, b_data] = LoadTEMPOData(para_batch{nn}{10},para_batch{nn}{11});
        
        if return_val == -1
            fprintf('*** HTB not found: %s ***\n',[para_batch{nn}{10},para_batch{nn}{11}]);
        end
        
        if return_val == -999 && para_batch{nn}{2} >= 5
            % I throw an error message only when we are in batch mode and we DO use the offline spike2 data.
            error('Htb & spike2 files not matched ...\n');
        end
        
        if (para_batch{nn}{13})
            select_data = g_data;
        else
            select_data = b_data;
        end
        
        if para_batch{nn}{6} == -1
            para_batch{nn}{6} = 1;
        end
        
        
        if para_batch{nn}{7} == -1
            para_batch{nn}{7} = size(g_data.event_data,3);
        end
        
        %get protocol type
        Protocol = g_data.one_time_params(PROTOCOLForParWorkers);		%get the protocol string descriptio
        
        %analyze data
        batch_flag = ori_filename;
        
        % --------- Do analysis finally ---------
        %                       1         2          3           4         5
        % para_batch{nn} = {Analysis, SpikeChan, SpikeChan2, StartCode, StopCode,...
        %                  6         7           8          9         10    11        12            13
        %               BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, UseSyncPulses, good_flag};
        Analysis_Switchyard(select_data, Protocol, para_batch{nn}{1}, para_batch{nn}{2}, para_batch{nn}{3}, para_batch{nn}{4}, para_batch{nn}{5}, ...
            para_batch{nn}{6}, para_batch{nn}{7}, para_batch{nn}{8}, para_batch{nn}{9}, para_batch{nn}{10}, para_batch{nn}{11}, batch_flag, para_batch{nn}{12});
        
    catch exception
        errorFiles(nn).fileName = [para_batch{nn}{11}];
        errorFiles(nn).info = exception;
        fprintf('*** Got error at %s ***\n',para_batch{nn}{11});
        fprintf('\t%s\n\tLine%g\n\t%s\n', exception.message, exception.stack(1).line, exception.stack(1).file)
    end
    
    parfor_progress;
end


%                     if print_flag
%                         printhandle = figure;
%                         close(printhandle);
%                         for printindex = 2:1:printhandle - 1
%                             print(printindex);
%                         end
%                     end
%                     if close_flag
%                         closehandle = figure;
%                         for closeindex = 2:1:closehandle
%                             close(closeindex);
%                         end
%                     end
%                 end

%if (line(1) ~= '%')...


fclose(fid);

% Close xls server. See "xlswrite1.m" for details
% try
%     invoke(Excel.ActiveWorkbook,'Save');
%     % Excel.Quit
%     Excel.delete
%     clear global Excel
% end

if print_flag == -999
    
    copyfile(batch_filename,[exportPath ori_filename]);
    
else
    if if_cluster
        outpath = ['Z/Data/TEMPO/Batch/' ori_filename(1:end-2) '/'];
    else
        outpath = ['Z:\Data\TEMPO\Batch\' ori_filename(1:end-2) '\'];
    end
    
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
    
    copyfile(batch_filename,[outpath ori_filename]);
end



% home;
realError = [];
if ~isempty([errorFiles(:).fileName])
    
    fprintf('\nOops, the following file(s) have error:\n\n');
    for i = 1:N
        if ~isempty(errorFiles(i).fileName)
            disp([num2str(i) ':  ' errorFiles(i).fileName]);
            realError = [realError, errorFiles(i)];
        end
    end
    
    assignin('base','BatchErrorInfo',realError);
 
    save([outpath 'BatchError.mat'], 'realError');
else
    fprintf('All success!\n');
end

toc

% audio = audioplayer(audioread('LabTools/Matlab/TEMPO_Analysis/type.wav'),20000);
% playblocking(audio);

% Window vibration, HaHa. HH20130829
% figure(1001);
% set(1001,'Unit','pixel');
% currPos = get(1001,'Position');
%
% for ii = 1:6
%     set(1001,'Position',[currPos(1) + 10*mod(ii,2) currPos(2) currPos(3) currPos(4)]);
%     drawnow;
%     pause(0.02);
% end
