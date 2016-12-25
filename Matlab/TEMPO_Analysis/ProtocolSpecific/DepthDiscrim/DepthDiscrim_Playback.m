%-----------------------------------------------------------------------------------------------------------------------
%-- DepthDiscrim_Playback.m -- Play back data from a discrimination run
%--	GCD, 3/5/01
%-----------------------------------------------------------------------------------------------------------------------
function DepthDiscrim_Playback(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

        DDiscrim_Playback_Loop(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
return;


%----------------------------------------------------------------------------------------------------------------------------------------------------
function DDiscrim_Playback_Loop(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

PREFERRED = 1;
NULL = 2;

SORT_CHOICE = 1;
SORT_DISPARITY = 2;
%sort_method = SORT_DISPARITY;
sort_method = SORT_CHOICE;
SELECT_CORR = 0;

%get the column of values of horiz. disparities in the dots_params matrix
h_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
unique_hdisp = munique(h_disp');

%get the binocular correlations
binoc_corr = data.dots_params(DOTS_BIN_CORR, :, PATCH1);
unique_bin_corr = munique(binoc_corr');

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);
unique_direc = munique(direction');

%get the column of values of speeds in the dots_params matrix
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

Pref_HDisp = data.one_time_params(PREFERRED_HDISP);

figure(2);
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [2 55 900 600], 'Name', 'Depth Discrimination');

rf_x = data.one_time_params(RF_XCTR);
rf_y = data.one_time_params(RF_YCTR);
rf_diam = data.one_time_params(RF_DIAMETER);

%first, we'll need the bin_width (in sec) of our spike raster + events + eye data
h = data.htb_header{SPIKE_DB};
spike_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
h = data.htb_header{EVENT_DB};
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
h = data.htb_header{EYE_DB};
eye_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

pref_dist = [NaN]; null_dist = [NaN];

%pull out trials for the selected correlation level
select_trials = find(binoc_corr == SELECT_CORR);

%loop through all the selected trials
for i = 1:length(select_trials)

    %find event bins corresponding to STIM_START and STIM_STOP
    start_eventbin = find(data.event_data(1,:,select_trials(i)) == VSTIM_ON_CD);
    stop_eventbin = find(data.event_data(1,:,select_trials(i)) == VSTIM_OFF_CD);
    stim_duration = (stop_eventbin - start_eventbin)*event_bin_width;
    
    %[binoc_corr(select_trials(i))    h_disp(select_trials(i))]
    spikes = [];
    spikes = data.spike_data(SpikeChan,:,select_trials(i));
    bad_ones = (spikes > 1); %shouldn't happen, but just in case
    spikes(logical(bad_ones)) = 1;
    
    %plot the spike train
    h_spikes = subplot('position', [.07 .77 .86 .2]);
    subplot(h_spikes);
    xtim = 1:length(spikes);
    h_raster = plot((xtim-data.htb_header{SPIKE_DB}.offset)*spike_bin_width, spikes, 'w-', 'LineWidth', 1.3);
    set(h_raster, 'EraseMode', 'none');
    hold on;
    plot([0 stim_duration], [1.2 1.2], 'k-', 'LineWidth', 4);
    hold off;
    xlim([-1 3]);
    ylim([0 1.5])
    xlabel('Time (seconds)', 'FontSize', 16);
    pause(0.5);
    
    %determine the monkey's choice on this trial
    temp = data.event_data(1,:,select_trials(i));
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice = PREFERRED;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice = NULL;
    end
    
    %widen the width of each spike pulse, to make them more audible
    ii = find(spikes == 1);
    spikes(ii + 1) = 1;
    spikes(ii + 2) = 1;
    spikes(ii + 3) = 1;
        
    %convert to audio sampling rate, 8192 Hz
    spikes_interp = interp(spikes, 8); %this is an approximation
    sound(8*spikes_interp);

    h_eyes = subplot('position', [.07 .08 .4 .6]);
    subplot(h_eyes);

    count = 0;
    t0 = clock;
    elap = etime(clock, t0);
    targ_x =   [data.targ_params(TARG_XCTR,select_trials(i),FP) data.targ_params(TARG_XCTR,select_trials(i),T1) data.targ_params(TARG_XCTR,select_trials(i),T2)];
    targ_y =   [data.targ_params(TARG_YCTR,select_trials(i),FP) data.targ_params(TARG_YCTR,select_trials(i),T1) data.targ_params(TARG_YCTR,select_trials(i),T2)];
    targ_wid = [data.targ_params(TARG_WIN_FWIDTH,select_trials(i),FP) data.targ_params(TARG_WIN_FWIDTH,select_trials(i),T1) data.targ_params(TARG_WIN_FWIDTH,select_trials(i),T2)];
    targ_hgt = [data.targ_params(TARG_WIN_FHEIGHT,select_trials(i),FP) data.targ_params(TARG_WIN_FHEIGHT,select_trials(i),T1) data.targ_params(TARG_WIN_FHEIGHT,select_trials(i),T2)];

    dots_x_start = []; dots_y_start = []; dots_x = []; dots_y = [];
    % compute some random dot locations
    num_dots = 50;
    dots_x_start = rand(num_dots, 1)*rf_diam + rf_x - rf_diam/2;
    dots_y_start = rand(num_dots, 1)*rf_diam + rf_y - rf_diam/2;
     
    hold on;
    h_cursor = plot(data.eye_data(LEYE_H,1,select_trials(i)), data.eye_data(LEYE_V,1,select_trials(i)), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    set(h_cursor, 'EraseMode', 'xor');

    h_dots = plot(dots_x_start, dots_y_start, 'wo', 'MarkerSize', 3, 'MarkerFaceColor', 'w');
    set(h_dots, 'EraseMode', 'xor');
    
    %plot a circle to show the RF location
    rectangle('Curvature', [1 1], 'Position',[(rf_x-rf_diam/2) (rf_y-rf_diam/2) rf_diam rf_diam]);
    %plot boxes to show the FP and targ windows
    for j=FP:T2
        rectangle('Position',[(targ_x(j)-targ_wid(j)/2) (targ_y(j)-targ_hgt(j)/2) targ_wid(j) targ_hgt(j)]);
    end
    %plot a marker for the fixation point
    plot(targ_x(FP), targ_y(FP), 'ks', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    %plot markers for choice targets
    h_T1 = plot(targ_x(T1), targ_y(T1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    h_T2 = plot(targ_x(T2), targ_y(T2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    
    %set the axis limits
    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('Horizontal Position (deg)', 'FontSize', 16);
    ylabel('Vertical Position (deg)', 'FontSize', 16);

    count = 1;
    
    while (elap < 3.9)

        %compute an index into the spike array
        spk_index = floor(elap/spike_bin_width + 1);
        if ( mod(count,4) == 0 )
            set(h_raster, 'XData', (xtim(1:spk_index)-data.htb_header{SPIKE_DB}.offset)*spike_bin_width, 'YData', spikes(1:spk_index), 'Color', 'b' );
        end

        %compute an index into the eye position array
        index = floor(elap/eye_bin_width + 1);
        %plot the eye cursor
        if (data.eye_data(LEYE_H,index,select_trials(i)) ~= 0)
            set(h_cursor, 'XData', data.eye_data(LEYE_H,index,select_trials(i)), 'YData', data.eye_data(LEYE_V,index,select_trials(i)));
        end
        %plot markers for the choice targets, and make the correct choice blink on and off
        if (h_disp(select_trials(i)) == Pref_HDisp)
            if ( (rem(elap,1) > 0.5) & (binoc_corr(select_trials(i)) ~= 0) )
                set(h_T1, 'MarkerEdgeColor', 'w', 'MarkerSize', 8, 'MarkerFaceColor', 'w' );
            else
                set(h_T1, 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'k' );
            end
        else
            if ( (rem(elap,1) > 0.5) & (binoc_corr(select_trials(i)) ~= 0) )
                set(h_T2, 'MarkerEdgeColor', 'w', 'MarkerSize', 8, 'MarkerFaceColor', 'w' );
            else
                set(h_T2, 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'k' );
            end
        end

        %during the visual stimulus period, plot text in RF to show stimulus polarity
        if ( (elap > start_eventbin*event_bin_width) & (elap < stop_eventbin*event_bin_width) )  %the visual stimulus period
            dots_x = dots_x_start + speed(select_trials(i))*(elap-start_eventbin*event_bin_width)*cos(direction(select_trials(i))*pi/180);
            while ( sum((dots_x>(rf_x+rf_diam/2))) > 0)
                dots_x(dots_x>(rf_x+rf_diam/2)) = dots_x(dots_x>(rf_x+rf_diam/2)) - rf_diam;
            end
            while ( sum((dots_x<(rf_x-rf_diam/2))) > 0)
                dots_x(dots_x<(rf_x-rf_diam/2)) = dots_x(dots_x<(rf_x-rf_diam/2)) + rf_diam;
            end

            dots_y = dots_y_start + speed(select_trials(i))*(elap-start_eventbin*event_bin_width)*sin(direction(select_trials(i))*pi/180);
            while ( sum((dots_y>(rf_y+rf_diam/2))) > 0)
                dots_y(dots_y>(rf_y+rf_diam/2)) = dots_y(dots_y>(rf_y+rf_diam/2)) - rf_diam;
            end
            while ( sum((dots_y<(rf_y-rf_diam/2))) > 0)
                dots_y(dots_y<(rf_y-rf_diam/2)) = dots_y(dots_y<(rf_y-rf_diam/2)) + rf_diam;
            end

            radius = sqrt((dots_x - rf_x).^2 + (dots_y - rf_y).^2);
            dots_x(radius > rf_diam/2) = NaN;
            dots_y(radius > rf_diam/2) = NaN;
            set(h_dots, 'XData', dots_x, 'YData', dots_y, 'MarkerEdgeColor', 'k', 'MarkerSize', 3, 'MarkerFaceColor', 'k' );

            if (binoc_corr(select_trials(i)) ~= 0)
                if (h_disp(select_trials(i)) < 0)
                    text(rf_x, rf_y+4.5, 'Near', 'FontSize', 18, 'HorizontalAlignment', 'Center');
                else
                    text(rf_x, rf_y+4.5, 'Far', 'FontSize', 18, 'HorizontalAlignment', 'Center');
                end            
            else
                text(rf_x, rf_y+4.5, '0% Corr.', 'FontSize', 18, 'HorizontalAlignment', 'Center');
            end
                
        end

        hold off;
        %force the figure to redraw
        drawnow;
        %get the current elapsed time
        elap = etime(clock, t0);
        count = count + 1;
    end
    cla;    
    
    %sort by disparity
    if (sort_method == SORT_DISPARITY)
        if (h_disp(select_trials(i)) == Pref_HDisp)
            pref_dist = [pref_dist spike_rates(select_trials(i))];
        else
            null_dist = [null_dist spike_rates(select_trials(i))];
        end
    end
    %sort by choice
    if (sort_method == SORT_CHOICE)
        if (choice == PREFERRED)
            pref_dist = [pref_dist spike_rates(select_trials(i))];
        else
            null_dist = [null_dist spike_rates(select_trials(i))];
        end
    end
    
    %update the response histogram
    h_curve = subplot('position', [.55 .08 .4 .6]);
    subplot(h_curve);    
    PlotTwoHists(pref_dist, null_dist);    
    xlabel('Mean Firing Rate (spikes/s)', 'FontSize', 16)
    ylabel('# of trials', 'FontSize', 16);
    
    subplot(h_spikes);
    cla;
    pause(1);

end


return;