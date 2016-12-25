%-----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	GCD, 6/27/03
%-- Modified for Vary_Fixation protocol  CRF, 1/06/04
%-----------------------------------------------------------------------------------------------------------------------

function MOOG_PSTH_Fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
temp_spike_data = data.spike_data(1,:);
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
fix_x     = temp_fix_x(~null_trials & select_trials);
fix_y     = temp_fix_y(~null_trials & select_trials);
spike_rates= temp_spike_rates(~null_trials & (trials >= BegTrial) & (trials <= EndTrial));
% notice that this bad_trials is the number without spon trials 
bad_trials = find(spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response

unique_azimuth  = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_fix_x    =  munique(fix_x');
unique_fix_y    =  munique(fix_y');

if length(unique_fix_y) == 1
   condition_num = fix_x;
   temp_condition_num = temp_fix_x;
else
   condition_num = fix_y; 
   temp_condition_num = temp_fix_y;
end

unique_condition_num = munique(condition_num');

% Assign titles
if length(unique_stim_type) == 1
    if unique_stim_type == 1
        h_title{1,1} = 'Vestib, EyePos = -20';
        h_title{1,2} = 'Vestib, EyePos = 0';
        h_title{1,3} = 'Vestib, EyePos = +20';
    else
        h_title{1,1} = 'Visual, EyePos = -20';
        h_title{1,2} = 'Visual, EyePos = 0';
        h_title{1,3} = 'Visual, EyePos = +20';
    end
else
    h_title{1,1} = 'Vestib, EyePos = -20';
    h_title{1,2} = 'Vestib, EyePos = 0';
    h_title{1,3} = 'Vestib, EyePos = +20';
    h_title{2,1} = 'Visual, EyePos = -20';
    h_title{2,2} = 'Visual, EyePos = 0';
    h_title{2,3} = 'Visual, EyePos = +20';
end

% add parameters here
% timebin for plot PSTH
timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% find spontaneous trials which azimuth,elevation,stim_type=-9999
spon_found = find(null_trials==1);     

% creat a new matrix without spontaneous trails from spike_data
hist_spike_data = [];                     

% take spontaneous activity out of the whole spike_data
for z=1:(length(spon_found)-1)
   hist_spike_data=[hist_spike_data, temp_spike_data( (frequency*spon_found(z)+1):(spon_found(z+1)-1)*frequency )];
end

% add the first part and the last part
hist_spike_data = [temp_spike_data( 1: frequency*(spon_found(1)-1) ), hist_spike_data, temp_spike_data( (frequency*spon_found(end)+1):end )]; 

% temp to get rid off bad trials, only temperally!!!!!!!
if ( bad_trials ~= NaN) 
    hist_spike_data = [hist_spike_data(1:(bad_trials-1)*frequency), hist_spike_data( (frequency*bad_trials+1):end )]; 
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % BEGIN TEMPORARY CODE (CRF 10/2006)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Temporarily compute only PSTH from preferred (max response) direction
% % (requires response matrix first)
% select = [];
% resp_mat = [];
% for n=1:length(unique_stim_type)
%     for i=1:length(unique_azimuth)
%         for j=1:length(unique_elevation)
%             for k=1:length(unique_condition_num)
%                 select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) & (stim_type==unique_stim_type(n)) ); 
%                 if (sum(select) > 0)
%                     resp_mat{n}(i,j,k) = mean(spike_rates(select));
%                 else
%                     resp_mat{n}(i,j,k) = resp_mat{n}(1,j,k);
%                 end
%             end
%         end
%     end
% end
% 
% for n = 1:length(unique_stim_type)
%     clear x y;
%     [x,y] = find(resp_mat{n}(:,:,2) == max(max(resp_mat{n}(:,:,2))));
%     ii{n} = x(1); jj{n} = y(1);
% end
% 
% % count spikes from raster data (spike_data)
% time_step=1;
% for n = 1:length(unique_stim_type)
%     select = logical( (azimuth==unique_azimuth(ii{n})) & (elevation==unique_elevation(jj{n})) & (condition_num==unique_condition_num(2)) & (stim_type==unique_stim_type(n)) ); 
%     % get rid off -90 and 90 cases
%     if (sum(select) > 0)
%        act_found = find( select==1 );
%        % count spikes per timebin on every same condition trials
%        for repeat=1:length(act_found) 
%            for m=1:(x_length)
%                temp_count(repeat,m)=sum(hist_spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+m*timebin)));
%                time_step=time_step+timebin;
%            end
%            time_step=1;
%         end
%         % get the average of the total same conditions if repetion is > 1
%         if (length(act_found) > 1);
%            count_y{n} = mean(temp_count);
%         else
%            count_y{n} = temp_count;     % for only one repetition cases
%         end
%        
%     else
%         disp('something wrong here');
%         dbstop;        
%     end
%     % change to spikes per second
%     rate_y{n} = count_y{n} * (1000/timebin);
% end
% 
% 
% % plot PSTH now
% % for n = 1:length(unique_stim_type)
% %     figure(n+10);  bar(rate_y{n});
% % end
% 
% 
% % output data
% if unique_stim_type == 1
%     stim_text = 'Vestibular';
% elseif unique_stim_type == 2
%     stim_text = 'Visual';
% else
%     stim_text = 'Vestibular'; % (both, but ves first)
% end
% 
% buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %s', ...
%        FILE, rate_y{1}, stim_text);
% outfile = ['C:\MATLAB6p5\work\new_PSTH_fix\PSTH_fix_3D.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% 
% % And another line if two stim types
% if length(unique_stim_type) > 1
%     stim_text = 'Visual';
%     buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %s', ...
%     FILE, rate_y{2}, stim_text);
%     outfile = ['C:\MATLAB6p5\work\new_PSTH_fix\PSTH_fix_3D.dat'];
% 	fid = fopen(outfile, 'a');
% 	fprintf(fid, '%s', buff);
% 	fprintf(fid, '\r\n');
% 	fclose(fid);
% end
% 
% return;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % END TEMP CODE, RESUME NORMAL PSTH CODE (CRF 10/2006)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% count spikes from raster data (spike_data)
time_step=1;
for n=1:length(unique_stim_type)
	for i=1:length(unique_azimuth)
        for j=1:length(unique_elevation)
            for k=1: length(unique_condition_num)
	%            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
                select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) & (stim_type==unique_stim_type(n)) ); 
                % get rid off -90 and 90 cases
                if (sum(select) > 0)
                    
                   act_found = find( select==1 );
                   % count spikes per timebin on every same condition trials
                   for repeat=1:length(act_found) 
                       for m=1:(x_length)
                           temp_count(repeat,m)=sum(hist_spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+m*timebin)));
                           time_step=time_step+timebin;
                       end
                       time_step=1;
                    end
                    % get the average of the total same conditions if repetion is > 1
                    if (length(act_found) > 1);
                       count_y{i,j,k,n} = mean(temp_count);
                    else
                       count_y{i,j,k,n}= temp_count;     % for only one repetition cases
                    end
                   
                 else
                    count_y{i,j,k,n}=count_y{1,j,k,n};
                 end
                 
            end
        end        
	end
end


% plot PSTH now

% get the largest count_y so that make the scale in each figures equal
max_count = max( cat(2, count_y{:}) ); 
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker =  [0,  max_count];

% define figure
for n=1:length(unique_stim_type)    
    figure(n+1);
    set(n+1,'Position', [5+(n-1)*100,25 1200,900], 'Name', '3D Direction Tuning');
    title(FILE);
    axis off;

    xoffset=0;
    yoffset=0;

    % now plot
    for k=1: length(unique_condition_num) 
        
        if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
            yoffset = yoffset-0.42;
            xoffset = 0;
        end
        % output some text 
        axes('position',[0 0 1 0.9]); 
        xlim([-50,50]);
        ylim([-50,50]);
        text(-30+xoffset*100,52+yoffset*110, h_title{n,k} );
        text(-43,-40, 'Azim: 0          45          90          135          180           225           270          315');    
        axis off;
        hold on;
        for i=1:length(unique_azimuth)
            for j=1:length(unique_elevation)
                axes('position',[0.05*i+0.02+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]); 
                bar( x_time,count_y{i,j,k,n}(1,:) ); 
                hold on;
                plot( x_start, y_marker, 'r-');
                plot( x_stop,  y_marker, 'r-');
                %set( gca, 'xticklabel', '0|on|20|off|40|' );
                % set the same scale for all plot
                xlim([0,x_length]);
                ylim([0,max_count]);
            end    
        end 

        xoffset=xoffset+0.46;
        
    end
end

return;

