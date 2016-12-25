%-----------------------------------------------------------------------------------------------------------------------
%-- MP_ISI.m -- Gets some specific inter-spike-interval information from each file
%-- Exists to identify poorly isolated neurons
%-- Started by JWN, 07/09/07
%-- Last by JWN, 07/09/07
%-----------------------------------------------------------------------------------------------------------------------

function MP_ISI(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

TEMPO_Defs;
ProtocolDefs;
Path_Defs;

disp(sprintf('(MP_ISI) Started at %s.',datestr(now,14)));

% Get monkey and cell numbers
[monkid, cellid, runstr]=strread(FILE,'m%dc%dr%s.htb');
% Calculate number of trials
trials = size(data.spike_data,3);

spikes = squeeze(data.spike_data(1,:,:));
catspikes = reshape(spikes, 1, 5000*trials);

counts = zeros(1,10);
counts(1) = sum(catspikes>1);
catspikes = cast(catspikes>0,'double');
counts(2) = sum(conv(catspikes,[1 1])>1);
counts(3) = sum(conv(catspikes,[1 -1 1])>1);
counts(4) = sum(conv(catspikes,[1 -1 -1 1])>1);
counts(5) = sum(conv(catspikes,[1 -1 -1 -1 1])>1);
counts(6) = sum(conv(catspikes,[1 -1 -1 -1 -1 1])>1);
counts(7) = sum(conv(catspikes,[1 -1 -1 -1 -1 -1 1])>1);
counts(8) = sum(conv(catspikes,[1 -1 -1 -1 -1 -1 -1 1])>1);
counts(9) = sum(conv(catspikes,[1 -1 -1 -1 -1 -1 -1 -1 1])>1);
counts(10) = sum(conv(catspikes,[1 -1 -1 -1 -1 -1 -1 -1 -1 1])>1);
counts = counts/trials;

% Write out to a file
PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
outfile = cell2mat(strcat(PATHOUT,{'ISI'},'.txt'));
headerflag = 0;
if (exist(outfile) == 0) % File does not yet exist, so print a header
    headerflag = 1;
end
fid = fopen(outfile, 'a');  % Open text file.
if (headerflag)
    fprintf(fid, 'FILE ');
    fprintf(fid, 'monkid cellid ');
    fprintf(fid, 'corrupts 1ms 2ms 3ms 4ms 5ms 6ms 7ms 8ms 9ms');
    fprintf(fid, '\r\n');
end
fprintf(fid,'%10s', strtok(FILE,'.'));
fprintf(fid,' %+2.5f', monkid, cellid, counts);
fprintf(fid,'\r\n');
fclose(fid);

disp('(MP_ISI) Done.');
return;