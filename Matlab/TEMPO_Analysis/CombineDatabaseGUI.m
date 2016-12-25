function fig = CombineDatabaseGUI()
% Combine two or more database GUI
% Merge file1.htb, file2.htb and so on ...
% Merge file1.log, file2.log and so on ...

global STOPCOMBINE
STOPCOMBINE = false;

set(0,'Units','characters')
Color = get(0,'DefaultUicontrolBackgroundcolor');
Position = [20 20 180 40];

% Create the GUI
gui.main = figure(...
    'HandleVisibility','on',...
    'IntegerHandle','off',...
    'Menubar','none',...
    'NumberTitle','off',...
    'Name','Combine Database',...
    'Tag','CombineDatabase',...
    'Color',Color,...
    'Units','characters',...
    'Position',Position); 

uicontrol(gui.main,'Style','text','FontSize',10,'HorizontalAlign','left','Units','characters','String','Input Data Files:','Position',[5 38 25 1.7]);

gui.listBox = uicontrol( gui.main,...
    'Style','listbox',...
    'FontSize',10,...
    'Units','characters',...
    'Max', 10, 'Min', 0,...
    'Position',[5 23 150 15]);

gui.addFile = uicontrol( gui.main,...
    'Style','push',...
    'FontSize',10,...
    'Units','characters',...
    'String','Add Data File',...
    'Position',[158 26 20 6],...
    'BackgroundColor','g',...
    'Callback',{@AddFile, gui.listBox});

gui.removeFile = uicontrol( gui.main,...
    'Style','push',...
    'FontSize',10,...
    'Units','characters',...
    'String','Remove File',...
    'Position',[158 23 20 1.7],...
    'Callback',{@RemoveFile, gui.listBox});

uicontrol(gui.main,'Style','text','FontSize',10,'HorizontalAlign','left','Units','characters','String','Input Data Files from Batch File:','Position',[5 19.5 60 1.7]);

gui.inputBatchFile = uicontrol(gui.main,...
    'Style','edit',...
    'FontSize',10,...
    'HorizontalAlign','left',...
    'Units','characters',...
    'String','',...
    'Position',[5 18 150 1.7]);

gui.browse1 = uicontrol( gui.main,...
    'Style','push',...
    'FontSize',10,...
    'Units','characters',...
    'String','Browse',...
    'Position',[158 18 20 1.7],...
    'Callback',{@Browse1});


uicontrol(gui.main,'Style','text','FontSize',10,'FontWeight','bold','HorizontalAlign','left','Units','characters','Position',[5 16 150 1.7],...
    'String','Notice: Output File Name will assign automatically! You only need select Output File Path.');


uicontrol(gui.main,'Style','text','FontSize',10,'HorizontalAlign','left','Units','characters','String','Output File Path:','Position',[5 11.5 25 1.7]);

gui.outputFilePath = uicontrol(gui.main,...
    'Style','edit',...
    'FontSize',10,...
    'HorizontalAlign','left',...
    'Units','characters',...
    'String',pwd,...
    'Position',[5 10 150 1.7]);

gui.browse2 = uicontrol( gui.main,...
    'Style','push',...
    'FontSize',10,...
    'Units','characters',...
    'String','Browse',...
    'Position',[158 10 20 1.7],...
    'Callback',{@Browse2});

uicontrol(gui.main,'Style','text','FontSize',10,'HorizontalAlign','left','Units','characters','String','Output File Name:','Position',[5 6.8 25 1.7]);

gui.outputFileName = uicontrol(gui.main,...
    'Style','edit',...
    'FontSize',10,...
    'HorizontalAlign','left',...
    'Units','characters',...
    'String','',...
    'Position',[27 7 50 1.7]);

gui.combine = uicontrol( gui.main,...
    'Style','push',...
    'FontSize',10,...
    'Units','characters',...
    'String','Combine',...
    'Position',[158 3 20 6],...
    'BackgroundColor','g',...
    'Callback',{@Combine});

gui.combine = uicontrol( gui.main,...
    'Style','push',...
    'FontSize',10,...
    'Units','characters',...
    'String','Stop',...
    'Position',[128 5 20 1.7],...
    'Callback',{@Stop});

gui.combine = uicontrol( gui.main,...
    'Style','push',...
    'FontSize',10,...
    'Units','characters',...
    'String','Manual',...
    'Position',[158 36 20 1.7],...
    'Callback',{@Manual});

gui.skipCheckbox = uicontrol( gui.main,...
    'Style','checkbox',...
    'Value', 1,...
    'FontSize',10,...
    'Units','characters',...
    'String','Skip FIXED_SEED checking',...
    'Tag','SKIP FIXED_SEED', ...
    'Position',[5 3 50 1.7]);

gui.skipCheckbox = uicontrol( gui.main,...
    'Style','checkbox',...
    'Value', 1,...
    'FontSize',10,...
    'Units','characters',...
    'String','Select Epoch(or Trial)',...
    'Tag','SELECT EPOCH', ...
    'Position',[50 3 50 1.7]);

gui.skipCheckbox = uicontrol( gui.main,...
    'Style','checkbox',...
    'Value', 1,...
    'FontSize',10,...
    'Units','characters',...
    'String','Show good trial number',...
    'Tag','SHOW GOOD TRIAL NUM', ...
    'Position',[87 3 50 1.7]);

gui.manageData = uicontrol( gui.main,...
    'Style','push',...
    'FontSize',10,...
    'Units','characters',...
    'String','Manage Data',...
    'Position',[158 33 20 1.7],...
    'Callback',{@manageData});

%% add mult file
    function AddFile(Object, eventdata, handles)
        s = get(gui.listBox, 'String');
        
        if ~isempty(s)
            file = cast(s(end),'char');
            [pathstr, name, ext, versn] = fileparts(file);
            filterSpec = sprintf('%s\\*.htb', pathstr);
        else
            filterSpec = 'Z:\Data\*.htb';
        end
 
        [filename, pathname, FilterIndex] = uigetfile(filterSpec,'Select .htb file', 'MultiSelect', 'on');
        if FilterIndex == 0 
            return
        end
        
        filename = strrep(filename,'.log','.htb');
        if iscell(filename)
            for i = 1:length(filename)
               tmp = strcat(pathname, filename(i));
               s = [s ; tmp];
            end
        else
            tmp = cellstr(strcat(pathname, filename));
            s = [s ; tmp];
        end

        set(gui.listBox, 'String', s);
        
        fileNames = [];
        [r c] = size(s);
        s = cast(s,'char');
        
        for i=1:r
            tmp = strread(s(i,:), '%s', 'delimiter', '\\');
            fileNames = [ fileNames ; tmp(end)];
        end

        if r>1
            outputFileName = findOutputFileName(fileNames);
        else
            outputFileName = cast(fileNames,'char');
            outputFileName = outputFileName(1:end-4);
        end
        set(gui.outputFileName, 'String', outputFileName);
        set(gui.outputFilePath, 'String', pathname); % HH20150210
        
        
    end

%% remove mult file
    function RemoveFile(Object, eventdata, handles)
        s = get(gui.listBox, 'String');

        if isempty(s)
            return;
        end

        v = get(gui.listBox, 'Value');
        if iscell(s)
            for i = 1:length(v)
                s(v(i)) = [];
                v = v - 1;
            end
        elseif v == 1
            s = [];
        end
        
        set(gui.listBox, 'String', s, 'Value', 1);
        
        fileNames = [];
        [r c] = size(s);
        s = cast(s,'char');
        for i=1:r
            tmp = strread(s(i,:), '%s', 'delimiter', '\\');
            fileNames = [ fileNames ; tmp(end)];
        end

        if r>1
            outputFileName = findOutputFileName(fileNames);
        else
            outputFileName = cast(fileNames,'char');
            outputFileName = outputFileName(1:end-4);
        end
        set(gui.outputFileName, 'String', outputFileName);
        
    end

%%  get input batch file
    function Browse1(Object, eventdata, handles)
        [filename, pathname] = uigetfile();
        if filename~=0
            file = [pathname filename];
            set(gui.inputBatchFile, 'String', file);
        end
    end

%%  get output file path
    function Browse2(Object, eventdata, handles)
        directory_name = uigetdir;
        if directory_name~=0
            set(gui.outputFilePath, 'String', directory_name);
        end
    end

%%  combine all input files to output file
    function Combine(Object, eventdata, handles)
        recordDiff = [];
        listFiles = get(gui.listBox, 'String');
        inputBatchFile = get(gui.inputBatchFile, 'String');
        outputFilePath = get(gui.outputFilePath, 'String');
        outputFileName = get(gui.outputFileName, 'String');
        
        if isempty(listFiles) && isempty(inputBatchFile)
            warndlg('Please add htb files or batch file!', '!! Warning !!');
        elseif ~isempty(listFiles) && ~isempty(inputBatchFile)
            button = questdlg('Which way to combine files?','Question','Files in List box','Batch File','Cancel', 'Cancel');
            if strcmp(button, 'Files in List box')
                recordDiff = CombineDatabase(listFiles, [outputFilePath '\' outputFileName]);
            elseif strcmp(button, 'Batch File')
                ReadBatchFile();
            end
        elseif ~isempty(listFiles)
            recordDiff = CombineDatabase(listFiles, [outputFilePath '\' outputFileName]);
        elseif ~isempty(inputBatchFile)
            ReadBatchFile();
        end
        
        if ~isempty(recordDiff)
            fid = fopen([outputFilePath '\CombineDatabase.log'], 'a');
            fprintf(fid, '\n\nCombine Database record: %s\n\n%s', datestr(now),recordDiff);
            fclose(fid);
            msgbox('Finished combine database in List box!','Confirm');
        end
    end

%% combine database according to Batch file
    function ReadBatchFile()
        inputBatchFile = get(gui.inputBatchFile, 'String');
        outputFilePath = get(gui.outputFilePath, 'String');
        fid = fopen([outputFilePath '\CombineDatabase.log'], 'a');
        fprintf(fid, '\n\nCombine Database record: %s\n\nBatch File:\n%s\n', datestr(now),inputBatchFile);
        fclose(fid);
        
        s = textread(inputBatchFile, '%s', 'whitespace', '', 'delimiter', '\n');
        c = cast(s,'char');
        for i = 1:length(s)
            if STOPCOMBINE
                STOPCOMBINE = false;
                return;
            end
            
            if c(i) == 'Z' || c(i) == 'z'
                for j = i:length(s)
                    if c(j) ~= 'Z' && c(j) ~= 'z'
                        break;
                    end
                end

                % find database that need combine
                if j-1 > i
                    listFiles = [];
                    fileNames = [];
                    % setup input list files
                    for k = i:j-1
                        cfile = cast(s(k),'char');
                        file = strread(cfile, '%s', 2);
                        fileNames = [fileNames; file(2)]; 
                        tmp = cellstr( [cast(file(1),'char') cast(file(2),'char')] );
                        listFiles = [listFiles; tmp];
                    end
                    
                    outputFileName = findOutputFileName(fileNames);
                                        
                    recordDiff = CombineDatabase(listFiles, [outputFilePath '\' outputFileName]);
                    if ~isempty(recordDiff)
                        fid = fopen([outputFilePath '\CombineDatabase.log'], 'a');
                        fprintf(fid, '%s', recordDiff);
                        fclose(fid);
                    end
                end
                
                i = j;
            end
        end
        
        msgbox('Finished combine all database in batch file!','Confirm');
    end

%% setup output file name, i.e. 'm4 c39 r4.htb'
    function outputFileName = findOutputFileName(fileNames)
        % monkey number
        filenum = length(fileNames);
        [token, remain] = strtok(fileNames, 'm');
        [token, remain] = strtok(token, 'c');
        diff = sum(~strcmp(token(1), token(2:filenum)));
        if diff
            outputFileName = ['m' cast(token(1),'char')];
            for k = 2:filenum
                outputFileName = [outputFileName '_' cast(token(k),'char')];
            end
        else
            outputFileName = ['m' cast(token(1),'char')];
        end
        % cell number
        [token, remain] = strtok(remain, 'c');
        [token, remain] = strtok(token, 'r');
        diff = sum(~strcmp(token(1), token(2:filenum)));
        if diff
            outputFileName = [outputFileName 'c' cast(token(1),'char')];
            for k = 2:filenum
                outputFileName = [outputFileName '_' cast(token(k),'char')];
            end
        else
            outputFileName = [outputFileName 'c' cast(token(1),'char')];
        end
        %r number
        [token, remain] = strtok(remain, 'r');
        [token, remain] = strtok(token, '.');
        diff = sum(~strcmp(token(1), token(2:filenum)));
        if diff
            outputFileName = [outputFileName 'r' cast(token(1),'char')];
            for k = 2:filenum
                outputFileName = [outputFileName '_' cast(token(k),'char')];
            end
        else
            outputFileName = [outputFileName 'r' cast(token(1),'char')];
        end
    end

%% load manual
    function Manual(Object, eventdata, handles)
        file = which('manual');
        command = sprintf('notepad ''%s'' &',file);
        dos(command);
    end

%%  stop combine input files from Batch file.
    function Stop(Object, eventdata, handles)
        button = questdlg('Do you want to stop combining database?','Confirm','Ok','Cancel','Cancel');
        if strcmp(button,'Ok')
            STOPCOMBINE = true;
        end
    end

%% manage each single data file by select the trial number
    function manageData(Object, eventdata, handles)
        listFiles = get(gui.listBox, 'String');
        outputFilePath = get(gui.outputFilePath, 'String');
        
        if isempty(listFiles)
            warndlg('Please add htb files!', '!! Warning !!');
        elseif ~isempty(listFiles)
            recordDiff = ManageSingleDatabase(listFiles, outputFilePath);
        end
        
        if ~isempty(recordDiff)
            fid = fopen([outputFilePath '\CombineDatabase.log'], 'a');
            fprintf(fid, '\n\nManage Database record: %s\n\n%s', datestr(now),recordDiff);
            fclose(fid);
            msgbox('Finished managing database in List box!','Confirm');
        end

    end
end
%directory_name = uigetdir
