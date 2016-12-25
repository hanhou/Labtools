%-----------------------------------------------------------------------------------------------------------------------
%-- MPGetNote.m -- Pulls a custom note from the batch note file. 
%-- Started by JWN, 3/21/05
%-----------------------------------------------------------------------------------------------------------------------
function note = MPGetNote(FILE,notefilepath);

disp(sprintf('(MPGetNote) Started at %s.',datestr(now,14)));

cid = strtok(FILE,'.');  % Remove .htb extension.
fid = fopen(notefilepath, 'r');  % Open text file.
note = 'Unknown cell id.';  % Default note.
y = 0;
while feof(fid) == 0
    tline = fgetl(fid);
    if(findstr(tline,cid))
        note = tline;
    end
end
    
disp('(MPGetNote) Done.');