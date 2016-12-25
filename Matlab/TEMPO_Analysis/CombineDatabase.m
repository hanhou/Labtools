function recordDiff = CombineDatabase(inputFiles, outputFile)
% Really combine database after get input files and output file.

warning('off','MATLAB:dispatcher:InexactMatch');

SkipHandle = findobj(gcbf, 'Tag', 'SKIP FIXED_SEED');
skipSelected = get(SkipHandle, 'Value');

selectEpochHandle = findobj(gcbf, 'Tag', 'SELECT EPOCH');
selectEpoch = get(selectEpochHandle, 'Value');

ShowGoodTrialHandle = findobj(gcbf, 'Tag', 'SHOW GOOD TRIAL NUM');
% showGoodTrialNum = get(ShowGoodTrialHandle, 'Value');
showGoodTrialNum = selectEpoch;

if strcmp(outputFile(end-3:end), 'log') || strcmp(outputFile(end-3:end), 'htb')
    outputFile = outputFile(1:end-4);
end

recordDiff = [];
htbFiles = cast(inputFiles, 'char');

[fileNumber c] = size(inputFiles);
if fileNumber < 2
    warndlg('Please add more than one input file!','!! Warning !!');
    return;
end

good_trials_all = cell(1,fileNumber);

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

recordDiff = sprintf('\n================================================\nInput Files:\n');
for i = 1:fileNumber
    recordDiff = sprintf('%s%s\n',recordDiff, logFiles(i,:));
end 
recordDiff = sprintf('%s\nOutput File:\n%s.log\n',recordDiff, outputFile);
    

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

header = CheckLogFileHeader();
if isempty(header)
    return;
end

if selectEpoch
    epochNumber = GetEpochNumber();
    recordDiff = sprintf('%s\nUser selected Epoch(Trial)#:\n',recordDiff);
    for i = 1:fileNumber
        lfn = strtrim(cast(logFileNames(i),'char'));
        tmpAns = str2num(cast(epochNumber(i,1:2),'char'));
        recordDiff = sprintf('%s%s: Start Epoch(Trial)# = %d; End Epoch(Trial)# = %d\n',recordDiff,lfn,tmpAns(1),tmpAns(2));
    end
end

CombineLogFiles(header);% fclose('all'); return;
CombineHtbFiles();
CombineSpike2Files(); % Added by HH20150210

if skipSelected
    recordDiff = sprintf('%s\nSkip FIXED_SEED checking\n',recordDiff);
end

recordDiff = sprintf('%s\nCombine successful!\n',recordDiff);
fclose('all');

%% check log file
    function header = CheckLogFileHeader()
        % find max line in log files' head
        header = [];
        maxline = 0;
        line = 0;
        for i = 1:fileNumber
            while feof(fidLogFiles(i)) == 0
                tline = fgetl(fidLogFiles(i));
                line = line + 1;
                if findstr(tline, 'TRIAL#')
                    if line > maxline
                        maxline = line;
                    end
                    line = 0;
                    break;
                end
            end
        end

        line = maxline - 1;

        logHeader = [];
        for i = 1:fileNumber
            s = textread(logFiles(i,:), '%s', 'whitespace', '', 'delimiter', '\n');
            logHeader = [logHeader ; s(1:line)];
        end

        result = ~strcmp(logHeader(1:line), logHeader(line+1 : 2*line));
        for j = 0:fileNumber-2
            for i = j+1:fileNumber-1
                result = ~strcmp(logHeader(j*line+1:(j+1)*line), logHeader(i*line+1 : (i+1)*line)) + result;
            end
        end

%         SkipHandle = findobj(gcbf, 'Tag', 'SKIP FIXED_SEED');
%         skipSelected = get(SkipHandle, 'Value');
        diffHeader = [];
        for i = 1:length(result)
            if result(i) > 0
                for j = 0:fileNumber-1
                    tmpstr = cast(logHeader(j*line+i), 'char');
                    tmpresult = strfind(tmpstr, 'FIXED_SEED');
                    if skipSelected && ~isempty(tmpresult)
                        continue;
                    end
                    diffHeader = [ diffHeader ; logHeader(j*line+i)]; 
                end
            end
        end

        if ~isempty(diffHeader)

            tmp = strtrim(logFiles(1,:));
            index = strfind(tmp, '\');
            tmp = tmp(1:index(end)-1);

            recordDiff = sprintf('%s\nFollowing is different between log file headers:\n',recordDiff);

%             title = sprintf('Different in log files: %s',Path);
            prompt = [];
            def = [];
            tmpStr = '';
            j = 1;
            for i = 1:length(diffHeader)
                lfn = strtrim(cast(logFileNames(j),'char'));
                dh = strtrim(cast(diffHeader(i),'char'));
                recordDiff = sprintf('%s%s: \t%s\n',recordDiff,lfn,dh);
                tmpStr = sprintf('%s%s: \t%s\n',tmpStr,lfn,dh);

                j = j+1;
                if j>fileNumber
                    j=1;
                    recordDiff = sprintf('%s\n',recordDiff);
                    prompt = [prompt; cellstr(tmpStr)];
                    def = [def; diffHeader(i-fileNumber+1)];
                    tmpStr = [];
                end
            end

            title_ = 'Conflict in log file?';
            answer = inputdlg(prompt,title_,1, def,'on');
            if isempty(answer)
                recordDiff = sprintf('%sUser Cancel combining database!',recordDiff);
                fclose('all');
                return
            else
                recordDiff = sprintf('%sUser selected header parameters:',recordDiff);
                for i = 1:length(answer)
                    dh = strtrim(cast(answer(i),'char'));
                    recordDiff = sprintf('%s\n%s',recordDiff,dh);
                end
                recordDiff = sprintf('%s\n',recordDiff);
            end
        end

        header = s(1:line);
        j = 1;
        for i = 1:length(result)
            if result(i) > 0
                tmpstr = cast(logHeader(i), 'char');
                tmpresult = strfind(tmpstr, 'FIXED_SEED');
                if skipSelected && ~isempty(tmpresult)
                    continue;
                end
                header(i) = answer(j);        
                j = j+1;
            end
        end
    end
%% Continue Combine all files
% first combine log files
    function CombineLogFiles(header)
        outputLogFile = [outputFile '.log'];
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


        for i = 1:length(header)
            tline = cast(header(i),'char');
            fprintf(fidApend, '%s\n', tline);
        end

        trial = 0;
        for fn = 1:fileNumber
            s = textread(logFiles(fn,:), '%s', 'whitespace', '', 'delimiter', '\n');

            for j = 1:length(s)
                tline = cast(s(j),'char');
                if findstr(tline, 'TRIAL#')
                    if selectEpoch
                        [a b] = strread(tline,'%s %d');
                        startEpoch = str2num(cast(epochNumber(fn,1),'char'));
                        if b == startEpoch
                            break;
                        end
                    else
                        break;
                    end
                end
            end

            for i = j:length(s)
                tline = cast(s(i),'char');
                if findstr(tline, 'TRIAL#')
                    if selectEpoch
                        [a b] = strread(tline,'%s %d');
                        endEpoch = str2num(cast(epochNumber(fn,2),'char'));
                        if b == (endEpoch+1)
                            break;
                        end
                    end
                    trial = trial + 1;
                    fprintf(fidApend, 'TRIAL#\t%d\n', trial);
                else
                    fprintf(fidApend, '%s\n', tline);
                end
            end
        end
    end


%% Then combine htb files
    function CombineHtbFiles()
        outputHtbFile = [outputFile '.htb'];
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

        ndbs = HTBCOUNT(fidHtbFiles(1));                  % Get number of databases in it
        for i = 2:fileNumber
            if ndbs ~= HTBCOUNT(fidHtbFiles(i));          % Get number of databases in it
                warndlg('Database numbers are different from different htb files. Cannot combine database!','!! Warning !!')
                fclose('all');
                return
            end
        end

        for database = 1:ndbs
            oldHtbHeader = HTBGETHD(fidHtbFiles(1), database);   % HTB header #1
            newHtbHeader = oldHtbHeader;
            
            if selectEpoch % user selected epoch number
                startNum = str2num(cast(epochNumber(1,1),'char'));
                endNum = str2num(cast(epochNumber(1,2),'char'));
                epoch = [];
                for nEpoch = startNum : endNum
                    oldEpoch = HTBGETEP(oldHtbHeader, nEpoch);
                    epoch = [epoch ; oldEpoch];
                end
                newHtbHeader.sweep = endNum - startNum + 1;
                newHtbHeader.alloc = (oldHtbHeader.alloc - 512)*(newHtbHeader.sweep/oldHtbHeader.sweep) + 512;
            else
                epoch = HTBGETDA(oldHtbHeader);
            end
            
            for i = 2:fileNumber
                oldHtbHeader = HTBGETHD(fidHtbFiles(i), database);   % HTB header
                tmpHeader = oldHtbHeader;
                
                if selectEpoch % user selected epoch number
                    startNum = str2num(cast(epochNumber(i,1),'char'));
                    endNum = str2num(cast(epochNumber(i,2),'char'));
                    for nEpoch = startNum : endNum
                        oldEpoch = HTBGETEP(oldHtbHeader, nEpoch);
                        epoch = [epoch ; oldEpoch];
                    end
                    tmpHeader.sweep = endNum - startNum + 1;
                    tmpHeader.alloc = (oldHtbHeader.alloc - 512)*(tmpHeader.sweep/oldHtbHeader.sweep) + 512;
                    
                    newHtbHeader.sweep = newHtbHeader.sweep + tmpHeader.sweep;
                    newHtbHeader.alloc = newHtbHeader.alloc + tmpHeader.alloc - 512; % reduce one header file
                else
                    newHtbHeader.sweep = newHtbHeader.sweep + oldHtbHeader.sweep;
                    newHtbHeader.alloc = newHtbHeader.alloc + oldHtbHeader.alloc - 512; % reduce one header file 

                    oldEpoch = HTBGETDA(oldHtbHeader);              % Get all epoch data for database 1
                    epoch = [epoch ; oldEpoch];
                end
            end
            
            newHtbHeader.fileoffset = ftell(fidAppend);
            err = HTBWRITE(fidAppend, newHtbHeader, epoch);
        end
    end

    function CombineSpike2Files()  % HH20150210

        fa = find(inputFiles{1} == '\');
        spike2_path = ['Z:\Data\MOOG', inputFiles{1}(fa(3):fa(4)), 'Analysis\SortedSpikes2\'];
        
        % Load data
        commonUnits = 0:1000; % A crazy buffer
        
        for ss = 1:fileNumber
            spike2_mat = [spike2_path inputFiles{ss}(fa(5)+1:end-4) '.mat'];
            
            if (exist(spike2_mat,'file'))
                origin_spsData2(ss) = load(spike2_mat);
                commonUnits = intersect(commonUnits,[origin_spsData2(ss).spsData2.UnitId]);
            else
                disp('At least one file does not have Spike2 data...');
                return;
            end
        end
        
        % Merge spike sorted data to new .mat file (only merge the units with the same ID)
        for cc = 1:length(commonUnits)
            thisUnit = find([origin_spsData2(1).spsData2.UnitId] == commonUnits(cc));
            spsData2(cc) = origin_spsData2(1).spsData2(thisUnit); % Copy other stuffs (only first file)
            
            % In case we only need part of trials for file #1
            startNum = str2double(epochNumber(1,1));
            endNum = str2double(epochNumber(1,2));
            spsData2(cc).spikeInfo = spsData2(cc).spikeInfo(find(good_trials_all{1}>=startNum,1):find(good_trials_all{1}<=endNum,1,'last'));
            
            for ss = 2:fileNumber % Merge SpikeTimes for file #>1
                startNum = str2double(epochNumber(ss,1));
                endNum = str2double(epochNumber(ss,2));

                thisUnit = find([origin_spsData2(ss).spsData2.UnitId] == commonUnits(cc));
                spsData2(cc).spikeInfo = [spsData2(cc).spikeInfo ,...
                    origin_spsData2(ss).spsData2(thisUnit).spikeInfo(find(good_trials_all{ss}>=startNum,1):find(good_trials_all{ss}<=endNum,1,'last'))];
            end
        end
        
        save([spike2_path outputFile(find(outputFile=='\',1,'last')+1:end) '.mat'],'spsData2');
        disp('Spike2 data merged...');
    end

%% Get Epoch number from each database
    function epochNumber = GetEpochNumber()
        epochNumber = [];
        for i = 1:fileNumber
            tmp = strtrim(htbFiles(i,:));
            index = strfind(tmp, '\');
            fileName = tmp(index(end)+1:end);

            if showGoodTrialNum
                listBoxTitle = sprintf('Good Trail Number: %s',fileName);
                %goodTrialList();
                %d = dir;
                %str = {d.name};
                
                [str, good_trials_all{i}] = goodTrialList(i);
                          
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
            end

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
           
%             if(ishandle(listBoxFig))
%                 close(listBoxFig);
%             end
            
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
            close(listBoxFig);
        end
        epochNumber = epochNumber';
    end

%% find good Trial number
    function [str,good_trials] = goodTrialList(currFileNum)
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