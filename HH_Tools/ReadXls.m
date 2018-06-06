function XlsData = ReadXls(FILE,Sheet,HEADS_N)

%% Read xls. HH20140624
% global num txt raw header;
try
    [num,txt,raw] = xlsread(FILE,Sheet);
catch err
    disp(err);
    warning('Xls read error... use basic mode');
    [num,txt,raw] = xlsread(FILE,Sheet,'','basic');
end
% [num,txt,raw] = xlsread('Z:\Data\MOOG\Results\Result.xlsm',2);

% Get Header infomation
% HEADS_N = 3;

if HEADS_N == 3 
    header_all = txt(HEADS_N-1,:);
    header_all(strcmp(header_all,'')) = txt(HEADS_N-2,strcmp(header_all,''));
else
    header_all = txt(1,:);
end

for i = 1:length(header_all)
    try
        if i == num(1,i)
            eval(['header.' header_all{i} '=' num2str(i) ';']);
            
            if sum(~isnan(num(HEADS_N+1:end,i))) > 0 % Number type
                eval(['types.' header_all{i} '= 1;']);
            else
                eval(['types.' header_all{i} '= 2;']);
            end
        else
            disp('Header info error...');
            keyboard;
        end
    catch
    end
end

% Delete headers

% end_line = find(~isnan(num(:,1)),1,'last'); 
if isfield(header,'Monkey')
    end_line = find(~isnan(num(:,header.Monkey)),1,'last'); 
else
    end_line = find(~isnan(num(:,1)),1,'last'); 
end

XlsData.num = num(HEADS_N+1:end_line,:);
XlsData.txt = txt(HEADS_N : end_line - 1,:); % Note here
XlsData.raw = raw(HEADS_N+1:end_line,:);


XlsData.header = header;
XlsData.type = types;
XlsData.hName = fieldnames(header);

