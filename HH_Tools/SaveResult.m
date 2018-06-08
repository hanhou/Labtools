function SaveResult(config, result)
% Output text for xls; save figures, .mat and .dat file (BATCH only)
% HH20141124

persistent XlsData; % Hold some xls data to prevent read xls file repeatly during batch processing. HH20150724

%%%%%%%%%%%%%%%%%%%%%  Output   HH20140510 / HH20140621 / HH20141003 /HH20141124 %%%%%%%%%%%%%%%%%

if ~isempty(config.batch_flag)  % Figures and raw data (always in "result" structure)
    
    outpath = ['Z:\Data\Tempo\Batch\' config.batch_flag(1:end-2) '\'];
    
    % Check directory
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
    
    dot_pos = strfind(result.FILE,'.');
    if ~isempty(dot_pos)
        result.FILE = result.FILE(1:dot_pos-1);
    end
    
    savefilename = [outpath [result.FILE '_' num2str(result.SpikeChan)] '_' config.suffix];
    
    
    % Overwrite or Append
    if exist([savefilename '.mat'],'file')
        if isfield(config,'append') && config.append  % Load already existing fields before deletion
            temp = load([savefilename '.mat']);
            result_old = temp.result;
            result_new = result;
            
            % Merge two structures
            M = [fieldnames(result_new)' fieldnames(result_old)' ; struct2cell(result_new)' struct2cell(result_old)' ];
            [~, rows] = unique(M(1,:), 'stable'); % Handle Conflicting ('stable' selects the first occurrence, i.e., result_new)
            M = M(:, rows);
            %              result = struct(M{:}); % This is wrong because we have cells in some fields...
            for i = 1:size(M,2)
                result.(M{i*2-1})= M{i*2};
            end
            disp('Appending fields in .mat file...');
        else
            delete([savefilename '*.*']);
            disp('Overwrite the .mat file...');
        end
    end
    
    % Save raw data
    save(savefilename,'result');
    disp('Saved to .mat ...');
    
    % Save figures
    for ff = 1:length(config.save_figures)
        orient landscape;
        set(config.save_figures(ff),'Visible','on','PaperPositionMode','auto','PaperOrientation','landscape');
        %         print(config.save_figures(ff),'-dbitmap',[savefilename '_fig_' num2str(config.save_figures(ff)) '.bmp']);
        
        saveas(config.save_figures(ff),[savefilename '_fig_' num2str(config.save_figures(ff))],'jpeg');
        saveas(config.save_figures(ff),[savefilename '_fig_' num2str(config.save_figures(ff))],'fig');
        
        if ~strcmp(config.batch_flag,'test.m')  % HH20160415
            close(config.save_figures(ff));
        end
    end
    disp('Saved to .fig...');
end


% Print part of data to texts (clipboard / .dat file / Data hub "Result.xlsm")

% try
    
    if ~isempty(config.batch_flag)
        
        % Print to file
        outfile = [outpath config.suffix '.dat'];
        printHead = 0;
        if (exist(outfile, 'file') == 0)   % file does not yet exist
            printHead = 1;
        end
        
        fid = fopen(outfile, 'a');
        % This line controls the output format
        
        if (printHead)
            fprintf(fid, ['FILE\t ' config.sprint_once_contents '|\t']);
            
            for ll = 1:length(config.sprint_loop_contents)
                fprintf(fid,[config.sprint_loop_contents{ll} '|\t']);
            end
            fprintf(fid, '\r\n');
        end
        
        fprintf(fid,'%s\t',[result.FILE '_' num2str(result.SpikeChan)]);
        
    else  % Print to screen
        fid = 1;
    end
    
    toClip = [];
%     toXls = {};
    
    
    sprint_once_marker = [];
    for i = 1:length(config.sprint_once_marker)
        sprint_once_marker = [sprint_once_marker '%' config.sprint_once_marker(i) '\t'];
    end
    
    % Print once
    if ~isempty(config.sprint_once_marker)
        eval(['buff = sprintf(sprint_once_marker,' config.sprint_once_contents ');']);
%         eval(['buff_toXls = {' config.sprint_once_contents '};' ]); % HH20150724

        fprintf(fid, '%s', buff);
        toClip = [toClip sprintf('%s', buff)];
%         toXls = [toXls, buff_toXls];
    end
    
    % Print loops
    for ll = 1:length(config.sprint_loop_contents)
        
        sprint_loop_marker = [];
        for i = 1:length(config.sprint_loop_marker{ll})
            sprint_loop_marker = [sprint_loop_marker '%' config.sprint_loop_marker{ll}(i) '\t '];
        end
        
        for condition = 1:3 % Always output 3 conditions (if not exist, fill with NaNs)
            if sum(result.unique_stim_type == condition)==0
                buff = sprintf(sprint_loop_marker,nan(1,sum(sprint_loop_marker=='%')));
%                 buff_toXls = num2cell(nan(1,sum(sprint_loop_marker=='%')));
            else
                k = find(result.unique_stim_type == condition);
                eval(['buff = sprintf(sprint_loop_marker,' config.sprint_loop_contents{ll} ');']);
%                 eval(['buff_toXls = {' config.sprint_loop_contents{1} '};' ]); % HH20150724
            end
            fprintf(fid, '%s', buff);
            toClip = [toClip sprintf('%s', buff)];
%             toXls = [toXls, buff_toXls]; % HH20150724
        end
        
    end
    
    fprintf(fid, '\r\n');
    
    if ~isempty(config.batch_flag)  % Close .dat file
        fclose(fid);
    end
    
    toClip = [toClip(1:end-1) sprintf('\r\n')]; % Delete the last '\t' for clipboard
    %     clipboard('copy',toClip);
    
    
    % Turn [NaN]s (number) into 'NaN's (string)
%     for ll = 1:length(toXls)
%         if isnan(toXls{ll})
%             toXls{ll} = 'NaN';
%         end
%     end

    % --- Save back to xls ---
    % (finally I decide to save somethings back to xls for easier and better visualization in xls). HH20150724
    
    if ~isempty(config.batch_flag) && isfield(config,'xls_column_begin')
        
        % Turn toClip into toXls (Separate strings by TAB and reoranize into cells)
        toXls = textscan(toClip,'%s','Delimiter','\t');
        toXls = toXls{1}';
        
        % Read xls if needed. (only for the first file in BATCH mode)
        if isempty(XlsData) || strcmp(config.batch_flag,'test.m')  % If we are in test mode, we reload xls each time. HH20160415
            XlsData = ReadXls('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm',2,3);
        end
        
        % Locate where to paste "toClip"
        row = find(strcmp(XlsData.txt(:,XlsData.header.FileNo),result.FILE)) + 3; 
        
        if ~isempty(row)
            if numel(row) > 1 % More than one SpikeChans
                row = intersect(row,find(XlsData.num(:,XlsData.header.Chan1) == result.SpikeChan)+3);
            else % Only one SpikeChan or no SpikeChan (Training sessions)
            end

            column_begin = XlsData.header.(config.xls_column_begin);
            column_end = XlsData.header.(config.xls_column_end);

    %         keyboard;

            % Write "toClip" into excel file
            if length(toXls) == column_end - column_begin + 1;
                column_begin_name = num2ExcelName(column_begin);
                range_name = [column_begin_name num2str(row)];
                try
                    xlswrite1('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm',toXls,2,range_name);  % Speed-up of xlswrite
                    disp('Writing to .xls finished...');
                catch
                    disp('Writing to .xls failed :<');
                end
            else
                disp('Size not match when write back to xls...');
                keyboard
            end
        else
            disp('No file entry found in .xls...');
        end
    end
    
% catch error
%     beep;
%     keyboard;
%     disp('Warning: Print error (it''s fine if you are only running part of the batch code)...');
% end

    if strcmp(config.batch_flag,'test.m')
        assignin('base','result',result);    
    end

end

function [col_str] = num2ExcelName(num_loc) % Convert xls column number to column name (such as "A", "AB", ...)
test = 2;
old = 0;
x = 0;
while test >= 1
    old = 26^x + old;
    test = num_loc/old;
    x = x + 1;
end
num_letters = x - 1;
str_array = zeros(1,num_letters);
for i = 1:num_letters
    loc = floor(num_loc/(26^(num_letters-i)));
    num_loc = num_loc - (loc*26^(num_letters-i));
    str_array(i) = char(65 + (loc - 1));
end
col_str = strcat(str_array(1:length(str_array)));
end