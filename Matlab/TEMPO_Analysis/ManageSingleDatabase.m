function recordDiff = ManageSingleDatabase(inputFiles, outputFilePath)
% Manage each database file by selecting trial number.

warning('off','MATLAB:dispatcher:InexactMatch');

recordDiff = [];
htbFiles = cast(inputFiles, 'char');

fileNumber = size(inputFiles, 1);

% open htb files
fidHtbFiles = [];
for i = 1:fileNumber
    fidHtbFiles(i) = fopen(htbFiles(i,:), 'r');
    if fidHtbFiles(i) == -1
        warndlg(['Cannot open htb file: ' htbFiles(i,:)],'!! Warning !!');
        fclose('all');
        return
    end  
end

% open log file
logFiles = [];
for i = 1:fileNumber
    tmpFile = strrep(htbFiles(i,:), '.htb', '.log');
    logFiles = [logFiles; tmpFile];
end

fidLogFiles = [];
for i = 1:fileNumber
    fidLogFiles(i) = fopen(logFiles(i,:), 'r');
    if fidLogFiles(i) == -1
        warndlg(['Cannot open log file: ' logFiles(i,:)],'!! Warning !!');
        fclose('all');
        return
    end
end

logFileNames = [];
for i = 1:fileNumber
    tmp = strread(logFiles(i,:), '%s', 'delimiter', '\\');
    logFileNames = [ logFileNames ; tmp(end)];
end

recordDiff = sprintf('\n================================================\n');
for i = 1:fileNumber
    recordDiff = sprintf('%s\nInput Files:%s\n',recordDiff, logFiles(i,:));
    outputLogFile = [outputFilePath '\' cast(logFileNames(i),'char')]; 
    recordDiff = sprintf('%s\nOutput File:%s\n',recordDiff, outputLogFile);
end 

epochNumber = GetEpochNumber();
recordDiff = sprintf('%s\nUser selected Epoch(Trial)#:\n',recordDiff);
for i = 1:fileNumber
    lfn = strtrim(cast(logFileNames(i),'char'));
    tmpAns = str2num(cast(epochNumber(i,1:2),'char'));
    recordDiff = sprintf('%s%s: Start Epoch(Trial)# = %d; End Epoch(Trial)# = %d\n',recordDiff,lfn,tmpAns(1),tmpAns(2));
end

for fn = 1:fileNumber
    CreateLogFile(fn);
    CreateHtbFile(fn);
end

recordDiff = sprintf('%s\nManage Database files successful!\n',recordDiff);
fclose('all');

%% create log file
% first create log file
    function CreateLogFile(fn)
        outputLogFile = [outputFilePath '\' cast(logFileNames(fn),'char')]; 
        if exist(outputLogFile, 'file') ~= 0
            tmp = sprintf('The following file has already exist!\n%s',outputLogFile);
            button = questdlg(tmp,'Output LOG file','Overwrite', 'Cancel', 'Cancel');
            if strcmp(button, 'Cancel')
                fclose('all');
                return
            end
        end

        % discard existing contents if file already exist
        fidApend = fopen(outputLogFile, 'w');

        trial = 0;
        s = textread(logFiles(fn,:), '%s', 'whitespace', '', 'delimiter', '\n', 'bufsize', 4095*2);

        % write log file header
        startEpoch = str2num(cast(epochNumber(fn,1),'char'));
        foundTRIAL = 0;
        for j = 1:length(s)
            tline = cast(s(j),'char');
            if findstr(tline, 'TRIAL#')
                foundTRIAL = 1;
                [a b] = strread(tline,'%s %d');
                if b == startEpoch
                    break;
                end
            elseif (foundTRIAL == 0)
                fprintf(fidApend, '%s\n', tline);
            end
        end

        % write log file each trial
        for i = j:length(s)
            tline = cast(s(i),'char');
            if findstr(tline, 'TRIAL#')
                [a b] = strread(tline,'%s %d');
                endEpoch = str2num(cast(epochNumber(fn,2),'char'));
                if b == (endEpoch+1)
                    break;
                end
                trial = trial + 1;
                fprintf(fidApend, 'TRIAL#\t%d\n', trial);
            else
                fprintf(fidApend, '%s\n', tline);
            end
        end

        fclose(fidApend);
    end


%% Then create htb files
    function CreateHtbFile(fn)
        outputLogFile = [outputFilePath '\' cast(logFileNames(fn),'char')];
        outputHtbFile = [outputLogFile(1:end-3) 'htb'];
        if exist(outputHtbFile, 'file') ~= 0
            tmp = sprintf('The following file has already exist!\n%s',outputHtbFile);
            button = questdlg(tmp,'Output HTB file','Overwrite', 'Cancel', 'Cancel');
            if strcmp(button, 'Cancel')
                fclose('all');
                return
            end
        end

        % discard existing contents if file already exist
        fidAppend = fopen(outputHtbFile, 'w');

        ndbs = HTBCOUNT(fidHtbFiles(fn));                  % Get number of databases in it
        for database = 1:ndbs
            oldHtbHeader = HTBGETHD(fidHtbFiles(fn), database);   % HTB header #1
            newHtbHeader = oldHtbHeader;
            
            startNum = str2num(cast(epochNumber(1,1),'char'));
            endNum = str2num(cast(epochNumber(1,2),'char'));
            epoch = [];
            for nEpoch = startNum : endNum
                oldEpoch = htbGetEp(oldHtbHeader, nEpoch);
                epoch = [epoch ; oldEpoch];
            end
            newHtbHeader.sweep = endNum - startNum + 1;
            newHtbHeader.alloc = (oldHtbHeader.alloc - 512)*(newHtbHeader.sweep/oldHtbHeader.sweep) + 512;
            % 512 bytes is the fixed header size
            
            newHtbHeader.fileoffset = ftell(fidAppend);
            err = HTBWRITE(fidAppend, newHtbHeader, epoch);
        end
    end

%% Get Epoch number from each database
    function epochNumber = GetEpochNumber()
        epochNumber = [];
        for i = 1:fileNumber
            tmp = strtrim(htbFiles(i,:));
            index = strfind(tmp, '\');
            fileName = tmp(index(end)+1:end);

            %if showGoodTrialNum
                listBoxTitle = sprintf('Good Trail Number: %s',fileName);
                str = goodTrialList(i);
                          
                listBoxFig = figure(...
                    'Name',listBoxTitle,...
                    'HandleVisibility','on',...
                    'IntegerHandle','off',...
                    'Menubar','none',...
                    'NumberTitle','off',...
                    'Units','characters',...
                    'Resize', 'off',...
                    'Position',[20 20 50 50]);
                listBox = uicontrol(listBoxFig,...
                'Style','listbox',...
                'FontSize',10,...
                'Units','characters',...
                'Max', 10, 'Min', 0,...
                'Position',[0 0 50 50],...
                'String', str);
            %end

            tmpstr = sprintf('File Name: %s\nStart Epoch number:',fileName);
            prompt = {tmpstr,'End Epoch number:'};
            dlg_title = 'Select Epoch';
            num_lines = 1;
            htbHeader = HTBGETHD(fidHtbFiles(i), 1);
            def = {'1', num2str(htbHeader.sweep)};

            s = struct('Resize', 'on', 'WindowStyle', 'normal');
            answer = inputdlg(prompt,dlg_title,num_lines,def, s);
            if isempty(answer)
                answer = def';
            end
           
            if(ishandle(listBoxFig))
                close(listBoxFig);
            end
            
            tmpAns = str2num(cast(answer,'char'));
            
            while tmpAns(1) > tmpAns(2) || tmpAns(1) < 1 || tmpAns(2) > htbHeader.sweep
                h = warndlg(['Input Epoch number wrong! The number must between 1 to ' num2str(htbHeader.sweep)],'!! Warning !!');
                uiwait(h);
                answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
                if isempty(answer)
                    answer = def';
                end
                tmpAns = str2num(cast(answer,'char'));
            end
            
            epochNumber = [epochNumber answer];
        end
        epochNumber = epochNumber';
    end

%% find good Trial number
    function str = goodTrialList(currFileNum)
        %SUCCESS_CD must be same as in TEMPO_Defs.m file
        SUCCESS_CD = 12;		%trial was successful
        
        str=[];
        ndbs = HTBCOUNT(fidHtbFiles(currFileNum));               % Get number of databases in it
        for database = 1:ndbs               % Loop though each one...
            %store the database headers in good_data and in bad_data
            temp_hd = HTBGETHD(fidHtbFiles(currFileNum), database);
            
            if strcmp(temp_hd.title, 'Events')
                %disp('reading event data...');
                temp3 = HTBGETDA(temp_hd);
                temp3 = reshape(temp3', [temp_hd.nchannels, temp_hd.period, temp_hd.sweep]);
            end
        end
        
        good_trials = find(temp3(:,:,:) == SUCCESS_CD);
        good_trials = ceil(good_trials/temp_hd.period);  %these are now trial indices
        
        for t=1:size(good_trials,1)
            tmp = sprintf('trial#: %d   good trial#: %d',good_trials(t),t);
            tmp = cellstr(tmp);
            str = [str; tmp];
        end
        str = cellstr(str);
        %str = num2str(good_trials, '%d\n');
    end
end