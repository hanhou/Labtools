%----------------------------------------------------------------------------------------------------------
%-- BATCH_GUI_Switchyard.m: This function contains the callback routines for all of the BATCH_GUI
%--	interface objects that require callback actions.  I have done it this way to keep all of the callback
%-- 	code in one centralized location.  JDN, 7/6/2000
%----------------------------------------------------------------------------------------------------------
function BATCH_GUI_Switchyard(action)

TEMPO_Defs;
ProtocolDefs;
Path_Defs;


switch(action)
    case 'file open'
        PathHandle = findobj(gcbf, 'Tag', 'Pathname');
        PATH = get(PathHandle, 'String');
        t = sprintf('%s*.m', PATH);
        [fname, pname] = uigetfile(t,'Select batch file');
        if (fname ~= 0)
            set(PathHandle, 'String', pname);
            FileHandle = findobj(gcbf, 'Tag', 'Filename');
            set(FileHandle, 'String', fname);
        end
        
        %set protocol names even if they aren't going to be used
        ProtocolHandle = findobj(gcbf, 'Tag', 'ProtocolList');
        set(ProtocolHandle, 'String', protocol_names);
        PROTOCOL = get(ProtocolHandle, 'Value');
        
        %set analysis options based upon protocol
        AnalysisHandle = findobj(gcbf, 'Tag', 'AnalysisList');
        set(AnalysisHandle, 'String', analysis_strings{PROTOCOL});
        
    case 'file info'
        PathHandle = findobj(gcbf, 'Tag', 'Pathname');
        PATH = get(PathHandle, 'String');
        FileHandle = findobj(gcbf, 'Tag', 'Filename');
        FILE = get(FileHandle, 'String');
        [return_val, info_text] = GetBatchInfo(PATH,FILE);
        ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
        if (return_val == 0)  %problem opening file
            listtext{length(listtext)+1} = 'ERROR: File could not be opened';
            set(ListHandle, 'String', listtext);
        else
            set(ListHandle, 'String', info_text);
        end
        set(ListHandle, 'Value',1);
        set(findobj(gcbf,'Tag','ExportFiles'),'enable','on');
        
        
    case 'export files'
        PathHandle = findobj(gcbf, 'Tag', 'Pathname');
        PATH = get(PathHandle, 'String');
        FileHandle = findobj(gcbf, 'Tag', 'Filename');
        FILE = get(FileHandle, 'String');
        
        print_flag = -999; % I use this to indicate that we only need to export the files (for data sharing). HH20141103
        BATCH_GUI_Tempo_Analysis(PATH, FILE, [], print_flag, [], []);
        
    case 'open batch file'
        PathHandle = findobj(gcbf, 'Tag', 'Pathname');
        PATH = get(PathHandle, 'String');
        FileHandle = findobj(gcbf, 'Tag', 'Filename');
        FILE = get(FileHandle, 'String');
        eval(sprintf('edit %s', [PATH FILE]));
        
    case  'Show Errors'
        PathHandle = findobj(gcbf, 'Tag', 'Pathname');
        PATH = get(PathHandle, 'String');
        FileHandle = findobj(gcbf, 'Tag', 'Filename');
        FILE = get(FileHandle, 'String');
        
        try
            load([PATH FILE(1:end-2) '\BatchError']);
            errors = length(errorFiles);
            home;
            for i = 1:errors
                disp(['Batch file Line ' num2str(errorFiles(i).line) ':  ' errorFiles(i).fileName]);
                if isfield(errorFiles(i),'info');
                    fprintf('       %s:  line %g , "%s"\n', errorFiles(i).info.stack(1).name,errorFiles(i).info.stack(1).line,errorFiles(i).info.message);
                end
            end
        catch
            disp('Congratulations! There seems no error for your last running...');
        end
        
        
    case 'load_mu_data_checkbox'
        global LOAD_MU_DATA;
        if get(findobj(gcbf, 'Tag', 'LoadMUDataCheckBox'), 'Value')
            LOAD_MU_DATA = 1;
        else
            LOAD_MU_DATA = 0;
        end
        
    case 'filter off'
        FilterOnHandle = findobj(gcbf, 'Tag', 'FilterOn');
        set(FilterOnHandle, 'Value', 0);
        
        %deacivate filter settings
        ProtocolHandle = findobj(gcbf, 'Tag', 'ProtocolList');
        set(ProtocolHandle, 'Enable', 'off');
        AnalysisBoxHandle = findobj(gcbf, 'Tag', 'AnalysisBox');
        set(AnalysisBoxHandle, 'Value', 0);
        set(AnalysisBoxHandle, 'Enable', 'off');
        
    case 'filter on'
        FilterOffHandle = findobj(gcbf, 'Tag', 'FilterOff');
        set(FilterOffHandle, 'Value', 0);
        
        %acivate filter settings
        ProtocolHandle = findobj(gcbf, 'Tag', 'ProtocolList');
        set(ProtocolHandle, 'Enable', 'on');
        AnalysisBoxHandle = findobj(gcbf, 'Tag', 'AnalysisBox');
        set(AnalysisBoxHandle, 'Enable', 'on');
        
        %setup analysis options based upon protocol
        AnalysisHandle = findobj(gcbf, 'Tag', 'AnalysisList');
        PROTOCOL = get(ProtocolHandle, 'Value');
        set(AnalysisHandle, 'String', analysis_strings{PROTOCOL});
        
    case 'change protocol'
        ProtocolHandle = findobj(gcbf, 'Tag', 'ProtocolList');
        AnalysisHandle = findobj(gcbf, 'Tag', 'AnalysisList');
        
        %setup analysis options based upon protocol
        PROTOCOL = get(ProtocolHandle, 'Value');
        set(AnalysisHandle, 'String', analysis_strings{PROTOCOL});
        
    case 'override analysis'
        AnalysisHandle = findobj(gcbf, 'Tag', 'AnalysisList');
        AnalysisBoxHandle = findobj(gcbf, 'Tag', 'AnalysisBox');
        override_flag = get(AnalysisBoxHandle, 'Value');
        if override_flag == 1
            set(AnalysisHandle, 'Enable', 'on');
        else
            set(AnalysisHandle, 'Enable', 'off');
        end
        
        
    case 'load data'
        %get file info
        PathHandle = findobj(gcbf, 'Tag', 'Pathname');
        PATH = get(PathHandle, 'String');
        FileHandle = findobj(gcbf, 'Tag', 'Filename');
        FILE = get(FileHandle, 'String');
        
        %get print and close info
        CloseHandle = findobj(gcbf, 'Tag', 'CloseBox');
        PrintHandle = findobj(gcbf, 'Tag', 'PrintBox');
        
        close_flag = get(CloseHandle, 'Value');
        print_flag = get(PrintHandle, 'Value');
        
        %print status
%         ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
%         set(ListHandle, 'String', 'Loading Data...');
        
        %get filter options
        FilterOffHandle = findobj(gcbf, 'Tag', 'FilterOff');
        filter_flag = get(FilterOffHandle, 'Value');
        
        if filter_flag ~=1
            ProtocolHandle = findobj(gcbf, 'Tag', 'ProtocolList');
            protocol_type = get(ProtocolHandle, 'Value') - 1;
            
            AnalysisBoxHandle = findobj(gcbf, 'Tag', 'AnalysisBox');
            analysis_flag = get(AnalysisBoxHandle, 'Value');
            if analysis_flag == 1
                AnalysisHandle = findobj(gcbf, 'Tag', 'AnalysisList');
                analysis_type = get(AnalysisHandle, 'Value');
            else
                analysis_type = -1;
            end
            
        else
            protocol_type = -1;
            analysis_type = -1;
        end
        
        %do analysis
        a=tic;
        
        BATCH_GUI_Tempo_Analysis(PATH, FILE, close_flag, print_flag, protocol_type, analysis_type);
        
        toc(a);
        
        %if (return_val == -1)  %problem opening file
        %    listtext{length(listtext)+1} = 'ERROR: File could not be opened';
        %    set(ListHandle, 'String', listtext);
        %end
end