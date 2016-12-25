%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND_PROT.M -- Takes in a batch file (from Z:\Data\Tempo\Batch Files), reads each  
% cell number, and looks for log files containing the desired protocol number.  
% Returns a list of filenames for which that protocol was run.   CRF 10/11/05
%
% EXAMPLE:   files = find_prot('vary_fix_horiz.m', 106)

function files_list = find_prot(batchfile, prot_num)

files = textread(['Z:\Data\Tempo\Batch Files\' batchfile],'%q','delimiter','\n');
files_orig = files;
index = 0;

for n = 1:length(files)
    
    spaces = find(files{n}==' ');
    files{n} = files{n}(spaces(1)+1:spaces(1)+12);
    files{n}(end-3:end) = '.log';
    files_orig{n} = files{n};
    
    for m = 1:9
        files{n}(8) = num2str(m);
        if exist(files{n}) == 2
            clear temp_prot;
            temp_prot = textread(files{n}, '%12c', 1);
            if temp_prot == ['PROTOCOL ' num2str(prot_num)]
                index = index+1;
                files_orig{n}(8) = num2str(m);
                files_orig{n}(end-3:end) = '.htb';
                files_list{index} = files_orig{n};  %  OR if you don't want a cell array:
%                files_list(index,:) = files_orig{n};
            end
        end
    end
    
end

files_list = files_list';

return;
