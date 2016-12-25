%-----------------------------------------------------------------------------------------------------------------------
%-- MovTargFit.m -- Plots the response during the movingtarget task versus spatial position 
%--                 and estimates a center of the response field.
%--	VR, 7/14/06 
%-----------------------------------------------------------------------------------------------------------------------

function MovTargFit(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

TEMPO_defs;
ProtocolDefs;

%get angle values (deg)
angle = data.dots_params(AP_OFF_ANG,:,PATCH1); 
angle = squeeze_angle(angle);
unique_angle = munique(angle');

%get eccentricity (radius) values - stored in AP_OFF_RAD
rad = data.dots_params(AP_OFF_RAD,:,PATCH1);
unique_rad = munique(rad');

%make list of trials
trials = 1:length(angle);
select_trials = logical (data.misc_params(OUTCOME,BegTrial:EndTrial) == CORRECT);

%get firing rates for delay period (VSTIM_ON:VSTIM_OFF) and saccade (VSTIM_OFF:IN_T1_WIN)
delay_rates = data.spike_rates(SpikeChan,:);
for i = trials
    targ_on(i) = find(data.event_data(1,:,i) == TARGS_ON_CD);
    fix_off(i) = find(data.event_data(1,:,i) == VSTIM_OFF_CD);
    in_T1(i) = find(data.event_data(1,:,i) == IN_T1_WIN_CD,1,'last');
    sacc_rates(i) = sum(data.spike_data(SpikeChan, fix_off(i):in_T1(i), i)) / length(fix_off(i):in_T1(i)) * 1000;
end

%note some times
delay_period = mean(fix_off-targ_on);

keyboard
%now plot the delay_rates
figh(1) = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', sprintf('%s: Delay Period RF Map',FILE));
subplot(411);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
subplot(4, 1, 2); hold on;

% errorbar for varing radius
if (length(unique_angle) == 1) %includes null condition
    %first compute rates
    for i = 1:length(unique_rad)
        py = delay_rates( (rad == unique_rad(i)) & select_trials );
        px = unique_rad(i).*ones(length(py),1);
        mean_delay_rate(i) = mean(py);
        std_delay_rate(i) = std(py);
        plot(px,py,'b*'); %change size of dots?
    end
    %now plot
    hold on;
    errorbar(unique_rad, mean_delay_rate, std_delay_rate, 'b-');
    xlabel('Radius (deg)');ylabel('Firing Rate (Hz)');

% errorbar for varying angle
elseif (length(unique_rad) == 1) %includes null condition
    %first compute rates
    for i = 1:length(unique_angle)
        py = delay_rates( (angle == unique_angle(i)) & select_trials );
        mean_delay_rate(i) = mean(py);
        std_delay_rate(i) = std(py);
    end
    plot_x = angle(select_trials);
    plot_y = delay_rates(select_trials);
    
    %now first shift the data so that the angle with the highest response
    %is in the middle of the distribution
    [m, max_ind] = max(mean_delay_rate);
    adjust_angle = unique_angle - 360.*(unique_angle-unique_angle(max_ind) > 180) + 360.*(unique_angle-unique_angle(max_ind) < -180);
    adjust_plot_x = plot_x - 360.*(plot_x-unique_angle(max_ind) > 180) + 360.*(plot_x-unique_angle(max_ind) < -180);
    
    %now plot
    hold on;
    plot(adjust_plot_x, plot_y, 'b*');
    errorbar(adjust_angle, mean_delay_rate, std_delay_rate, 'bo');
    xlabel('Angle (deg)'); ylabel('Firing Rate (Hz)');
    axis('tight');
    
    %now fit the data to a gaussian and plot the result
    means = [adjust_angle mean_delay_rate'];
    raw = [adjust_plot_x' plot_y'];
    [pars] = gaussfit(means, raw, 0);   %last arg: allow positive going fit only
    x_interp = (getval(xlim,1): 0.5 : getval(xlim,2));
    y_interp = gaussfunc(x_interp - 360.*(x_interp-unique_angle(max_ind) > 180) + 360.*(x_interp-unique_angle(max_ind) < -180),pars);
    plot(x_interp, y_interp, 'r-');
    
    %now make some psths
    binwidth = 20;
    for i=1:length(unique_angle)
        select = trials(select_trials & (angle == unique_angle(i)));
        for j = 1:length(select)
            infix = find(data.event_data(1,:,select(j))==IN_FIX_WIN_CD);
            fixoff = find(data.event_data(1,:,select(j))==VSTIM_OFF_CD);
            [bins binned_rasters{i}(j,:)] = spikebinner(data.spike_data(SpikeChan, infix:infix+1525, select(j)), 1, binwidth, 100);
        end
        psth(i,:) = sum(binned_rasters{i},1)./length(select)./binwidth.*1000;
    end
    %plot these psths on one curve
    subplot(413); hold on;
    cm = colormap(hsv);
    spacer = floor(size(cm,1)/length(unique_angle));
    linecolors = cm(1:spacer:size(cm,1),:);
    legstr = 'legend(legh,';
    for i = 1:length(unique_angle)
        legh(i) = plot(bins,psth(i,:),'Color',linecolors(i,:));
        legstr = strcat(legstr,sprintf('''%d'',',unique_angle(i)));
    end
    legstr = strcat(legstr,'''Location'',''EastOutside'');');
    %eval(legstr);
    xlabel('Time About Target Onset');
    ylabel('F.R.(Hz)');
    axis('tight');
    plot([0 0],ylim,'k');
    plot([800 800],ylim,'k:')
    subplot(412); hold on;
    for i = 1:length(unique_angle)
        plot([adjust_angle(i)-5 adjust_angle(i)+5],[max(ylim) max(ylim)],'Color',linecolors(i,:),'LineWidth',4);
    end
    
    %now compute and plot the late delay data
    %first compute rates
    for i = 1:length(trials)
        fixoff = find(data.event_data(1,:,i)==VSTIM_OFF_CD);
        late_delay_rates(i) = sum(data.spike_data(SpikeChan, fixoff-300:fixoff, i)) / 301 * 1000;
    end
    for i = 1:length(unique_angle)
        py = late_delay_rates( (angle == unique_angle(i)) & select_trials );
        mean_late_delay_rate(i) = mean(py);
        std_late_delay_rate(i) = std(py);
    end
    plot_y = late_delay_rates(select_trials);
%     keyboard
    subplot(414);hold on;
    plot(adjust_plot_x, plot_y, 'b*');
    errorbar(adjust_angle, mean_late_delay_rate, std_late_delay_rate, 'bo');
    axis('tight')
    xlabel('Angle (deg)'); ylabel('Firing Rate (Hz)');
        
    %now fit the data to a gaussian and plot the result
    means = [adjust_angle mean_late_delay_rate'];
    raw = [adjust_plot_x' plot_y'];
    [pars] = gaussfit(means, raw, 0);   %last arg: allow positive going fit only
    x_interp = (getval(xlim,1): 0.5 : getval(xlim,2));
    y_interp = gaussfunc(x_interp - 360.*(x_interp-unique_angle(max_ind) > 180) + 360.*(x_interp-unique_angle(max_ind) < -180),pars);
    plot(x_interp, y_interp, 'r-');
    
    
    
% 2d colored scatter plot
else                        
    subplot(312);hold on;
    %first, convert polar coordinates into cartesian
    [xpos ypos] = pol2cart(angle.*pi./180, rad);
    xpos = round(1e6.*xpos)./1e6;  
    ypos = round(1e6.*ypos)./1e6;  
    unique_pos = munique([xpos' ypos']); %combining since xpos and ypos are not independent
    
    %next compute mean firing rate at each location
    for i = 1:size(unique_pos,1)
        tempx = unique_pos(i,1);   
        tempy = unique_pos(i,2);
        mean_delay_rate(i) = mean(delay_rates( (xpos==tempx) & (ypos==tempy) & select_trials ));
        std_delay_rate(i) = std(delay_rates( (xpos==tempx) & (ypos==tempy) & select_trials ));
        mean_sacc_rate(i) = mean(sacc_rates( (xpos==tempx) & (ypos==tempy) & select_trials ));
        std_sacc_rate(i) = std(sacc_rates( (xpos==tempx) & (ypos==tempy) & select_trials ));
    end
    
    %now plot a colored scatter plot
    temph = scatter(unique_pos(:,1), unique_pos(:,2), 80, mean_delay_rate); %3rd param is the marker size
    grid on;
    set(temph,'MarkerFaceColor','flat');
    mxspk = max(mean_delay_rate);
    cb = colorbar; set(cb,'YTickLabel',1e-2.*round(1e2.*[mxspk/10:mxspk/10:mxspk]));
    xl = xlim; yl = ylim; %store the plot dimensions
    
    %now fit the firing rates in angle/rad space (NOT cartesian) space with a 2d gaussian 
%     figh(2) = figure;
%     subplot(211);  
    hold on;
    
    fitrates = [];
    for i = 1:length(unique_angle)
        for j = 1:length(unique_rad)
            fitrates(i,j) = mean(delay_rates( (angle==unique_angle(i)) & (rad == unique_rad(j)) & select_trials));
        end
    end
    %contourf(unique_angle, unique_rad, fitrates);
    %organize delay_rates by 
    [ang_list rad_list] = cart2pol(unique_pos(:,1),unique_pos(:,2));
    ang_list = ang_list./pi.*180;
    ang_list = squeeze_angle(ang_list);
    ang_list = 1e-4.*round(1e4.*ang_list);
    rad_list = 1e-4.*round(1e4.*rad_list);
    
    raw = [angle(select_trials)' rad(select_trials)' delay_rates(select_trials)'];
    means = [ang_list rad_list mean_delay_rate'];
    pars = gauss2Dfit(means,raw);
    
    %create interpolated arrays for data display
    x_interp = [xl(1):0.5:xl(2)];
    y_interp = [yl(1):0.5:yl(2)];
    z_gauss = zeros(length(x_interp), length(y_interp));
    
    %obtain fitted data for interpolated arrays
    for i=1:length(x_interp)
        for j = 1:length(y_interp)
            [ang_temp rad_temp] = cart2pol(x_interp(i), y_interp(j));
            ang_temp = ang_temp/pi*180;
            z_gauss(i,j) =  gauss2Dfunc(ang_temp, rad_temp, pars);
        end
    end
    subplot(313);
    contourf(x_interp, y_interp, z_gauss')
    colorbar
%     axis image

    figure;
    set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', sprintf('%s: Delayed Saccade Fit',FILE));
    subplot(211)
    contourf(unique_angle, unique_rad, fitrates');
    %contourf(ang_list, rad_list, mean_delay_rate');
    colorbar;
    title(sprintf('%s: Data',FILE));
    subplot(212)
    ang_interp = [min(ang_list):1:max(ang_list)];
    rad_interp = [min(rad_list):1:max(rad_list)];
    z_gauss2 = zeros(ang_interp, rad_interp);
    for i=1:length(ang_interp)
        for j=1:length(rad_interp)
            z_gauss2(i,j) = gauss2dfunc(ang_interp(i), rad_interp(j), pars);
        end
    end
    contourf(ang_interp, rad_interp, z_gauss2');
    colorbar;
    title(sprintf('Peak at (%5.1f, %5.1f)',pars(3), pars(5)));
    
%     keyboard
end



% keyboard;