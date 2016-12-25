%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND_PROT_DIR.M -- Searches the supplied directory for .log files 
% that contain the desired protocol number.  Returns a list of 
% filenames for which that protocol was run.   CRF 10/14/05
% 
% EXAMPLE:   files = find_prot_dir('z:\data\moog\que\raw\', 106)

function files_list = find_prot(dirpath, prot_num)

files = dir(dirpath);
index = 0;

for n = 1:length(files)
    if n/10 == floor(n/10)
        n
        out_of = length(files)
    end
    files_orig{n} = files(n).name;
    if length(files(n).name) > 9
        if files(n).name(end-3:end) == '.log'
            temp_prot = textread([dirpath files(n).name], '%12c', 1);
            if temp_prot == ['PROTOCOL ' num2str(prot_num)]
                index = index+1;
                files_orig{n}(end-3:end) = '.htb';
                files_list{index} = files_orig{n};
            end
        end
    end
end

files_list = files_list';
    
return;
