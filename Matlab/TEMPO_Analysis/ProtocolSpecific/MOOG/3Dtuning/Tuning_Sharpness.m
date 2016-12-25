%----------------------------------------------------------------------------------------------------------------------
%-- Tuning_Sharpness.m -- Computes different metrics of tuning sharpness (steepness),
%-- to compare cells pre- and post- discrim training.  - CRF, 4/2007
%-----------------------------------------------------------------------------------------------------------------------

function Tuning_Sharpness(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_fix_x = data.moog_params(FIX_X,:,MOOG);
temp_fix_x(isnan(temp_fix_x)) = 0;
temp_spike_rates = data.spike_rates(SpikeChan, :);    

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_trials ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
fix_x = temp_fix_x(~null_trials & select_trials);
spike_rates= temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_fix_x = munique(fix_x');

if Protocol == 100
    trials_per_rep = 79;
elseif Protocol == 104 
    trials_per_rep = 81;
elseif Protocol == 107
    trials_per_rep = length(unique_azimuth)*length(unique_stim_type);
elseif Protocol == 108
    trials_per_rep = 93;
else
    disp('ERROR -- Protocol not recognized');
    return
end

num_reps = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);

relevant_azimuths = find(unique_azimuth >= 45 & unique_azimuth <= 135)';

resp_mat = []; resp_mean = []; resp_std = [];
for j = 1:length(unique_stim_type)
    for i = 1:length(relevant_azimuths)
        select = logical( stim_type==unique_stim_type(j) & azimuth==unique_azimuth(relevant_azimuths(i)) & elevation==0 & fix_x==0 );
        if (sum(select) > 0)                
            resp_mat{i,j} = spike_rates(select)';
        else
            resp_mat{i,j} = zeros(num_reps,1);
        end
        resp_mat{i,j} = resp_mat{i,j}(1:num_reps);  % remove the incomplete reps
        resp_mean(i,j) = mean(resp_mat{i,j});
        resp_std(i,j) = std(resp_mat{i,j});
    end
end

for j = 1:length(unique_stim_type)
    y = horzcat(resp_mat{:,j});
    temp = size(y);
    for l = 1:temp(1)
        x(l,:) = unique_azimuth(relevant_azimuths)';
    end
    p = polyfit(x,y,1);
    tuning_slope(j)= p(1);
    f = polyval(p,x);
    [c,P] = corrcoef(f,y);
    tuning_r(j) = c(1,2);
    tuning_p(j) = P(1,2);
end

%---------------------------------------------------------------------------------------
outfile = ['C:\MATLAB6p5\work\Tuning_Sharpness\sharpness.dat'];
for j = 1:length(unique_stim_type)
	buff = sprintf('%s\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t', FILE, unique_stim_type(j), tuning_slope(j), tuning_r(j), tuning_p(j));
	printflag = 0;
	if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
	end
	fid = fopen(outfile, 'a');
	if (printflag)
        fprintf(fid, 'FILE\t stim_type\t slope\t r\t p\t');
        fprintf(fid, '\r\n');
	end
	
	fprintf(fid, '%s', buff);
	fprintf(fid, '\r\n');
	fclose(fid);
end

return;