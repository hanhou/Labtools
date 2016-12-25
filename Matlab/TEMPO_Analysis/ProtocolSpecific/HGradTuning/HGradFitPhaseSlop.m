function HGradFitPhaseSlop(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE)

TEMPO_Defs;

symbols = {'bo' 'ro' 'go' 'ko' 'c*' 'm*' 'b*' 'r*' 'g*'};
lines = {'b-' 'r-' 'g-' 'k-' 'c-' 'm-' 'b-.' 'r-.' 'g-.'};
lines2 = {'b--' 'r--' 'g--' 'k--' 'c--' 'm--' 'b:' 'r:' 'g:'};
colors = {[0 0 1] [1 0 0] [0 1 0] [0 0 0] [0 1 1] [1 0 1] [0 0 1] [1 0 0] [0 1 0]};


%Start Data Retrieval Routines---------------------------------------------------------------------------------------------------------
%get the column of values of horiz. disparity magnitude in the dots_params matrix
mag_disp = data.dots_params(DOTS_HGRAD_MAG,BegTrial:EndTrial,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (mag_disp == data.one_time_params(NULL_VALUE)) );

unique_mag_disp = munique(mag_disp(~null_trials)');	

%get the column of values of horiz. disparity angle of orientation in the dots_params matrix
disp_ang = data.dots_params(DOTS_HGRAD_ANGLE,BegTrial:EndTrial,PATCH1);
unique_disp_ang = munique(disp_ang(~null_trials)');


%get the column of mean disparity values
mean_disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (mean_disp == LEYE_CONTROL) | (mean_disp == REYE_CONTROL) | (mean_disp == UNCORR_CONTROL) );

unique_mean_disp = munique(mean_disp(~null_trials & ~control_trials)');

%get the column of different aperture sizes
ap_size = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH1);

%do all sizes
all_sizes = 0
unique_ap_size = munique(ap_size(~null_trials)');
if all_sizes ~= 1
    unique_ap_size = unique_ap_size(length(unique_ap_size));
    num_ap_size = length(unique_ap_size);
else
    num_ap_size = length(unique_ap_size);
end

%unique_ap_size = munique(ap_size(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get the average horizontal eye positions to calculate vergence
Leyex_positions = data.eye_positions(1, :);
Reyex_positions = data.eye_positions(3, :);

vergence = Leyex_positions - Reyex_positions;


%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(mag_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%End Data Retrieval Routines---------------------------------------------------------------------------------------------------------
num_ap_size = length(unique_ap_size);
num_mag_disp = length(unique_mag_disp);
graph = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [500 50 500 773], 'Name', 'Fitted Tilt Tuning Curves');

stat_out = '';
f_out = '';
checkfile = ['Z:\Users\jerry\GradAnalysis\figure_data\simul_r_sqared_mdisp.dat'];
if (exist(checkfile, 'file') == 0)    %file does not yet exist
    stat_out{1} = sprintf('File\tMDisp\tR\tTDI\tChiP\tRall\tIndErr\tConstErr\tPseq\n');
end
font_size = 8;
bump_size = 10;

num_free_params = 5;
p_val = zeros(length(unique_ap_size), length(unique_mean_disp));
pref_tilt = zeros(length(unique_ap_size), length(unique_mean_disp));
TDI_save = zeros(1,(length(unique_ap_size)));
for i=1:length(unique_ap_size)
    TDIdata = [];
    if length(unique_mag_disp) < length(unique_mean_disp)
        for j=1:length(unique_mag_disp)
            start = zeros(length(unique_mean_disp), 1);
            stop = zeros(length(unique_mean_disp), 1);
            for k=1:length(unique_mean_disp)
                figure(graph);
                hold on;
                subplot(num_ap_size*2, num_mag_disp,  ((j-1)*(num_mag_disp) + i)*2);
                disp_select = logical((ap_size == unique_ap_size(i)) & (mag_disp == unique_mag_disp(j)) & (mean_disp == unique_mean_disp(k)) );
                
                plot_x = disp_ang(disp_select & ~null_trials & ~control_trials & select_trials);
                plot_y = spike_rates(disp_select & ~null_trials & ~control_trials & select_trials);
                
                %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
                [px, py, perr, spk_max, spk_min] = PlotTuningCurve(plot_x', plot_y', symbols{k}, '', 1, 1);
                [TDI(k), var_term] = Compute_DDI(plot_x, plot_y);
                
                %store data to calculate adjusted TDI later
                start(k) = length(TDIdata)+1;
                stop(k) = length(plot_x)+start(k)-1;
                TDIdata(start(k):stop(k), 1) = plot_x';
                TDIdata(start(k):stop(k), 2) = plot_y';
                
                x_interp = (px(1)): 1 : (px(length(px)));
                px = (px * pi)/180;
                plot_x = (plot_x * pi)/180;
                means{k} = [px py];
                raw{k} = [plot_x' plot_y'];
                
                ind_means = [px py];
                ind_raw = [plot_x' plot_y'];
                
                %fit with a distorted sin wave
                ind_pars{k} = sin_exp_fit(ind_means,ind_raw);
                ind_sinerr(k) = sin_exp_err(ind_pars{k})
                
                %grab pref_tilt of indp fits
                x_interp = (px(1)): .01 : (px(length(px)));
                x_deg = (x_interp * 180)/pi;
                y_sin = sin_exp_func(x_interp, ind_pars{k});
                y_err = sin_exp_err(ind_pars{k});
                y_sin(y_sin < 0) = 0;
                hold on
                plot(x_deg, y_sin, lines2{k});
                
                %store p-values of each curve
                temp_x =(plot_x *180)/pi;
                p_val(i,k) = calc_mdisp_anovap(disp_select, temp_x, plot_y, unique_disp_ang);
                [value, index_max] = max(y_sin);
                pref_tilt(i,k) = x_deg(index_max);
                
                %run chi^2 test on fit
                [chi2(k), chiP(k)] = Chi2_Test(ind_raw(:,1), ind_raw(:,2), 'sin_exp_func', ind_pars{k}, num_free_params);
                
                x_raw{k} = plot_x;
                y_raw{k} = plot_y;
                
                x_means{k} = px;
                y_means{k} = py;
            end %end mdisp
            %readjust mean disparity responses to fall on the same mean
            %then calc avg TDI
            %shifted_graphs = figure;
            %set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [750 50 500 773], 'Name', 'Mean Adjusted Tilt Tuning Curves');
            total_mean = mean(TDIdata(:,2));
            for count_meandisp = 1:length(unique_mean_disp)
                disp_mean = mean(TDIdata(start(count_meandisp):stop(count_meandisp),2));
                difference = total_mean - disp_mean;
                TDIdata(start(count_meandisp):stop(count_meandisp),2) = TDIdata(start(count_meandisp):stop(count_meandisp),2) + difference;
            end
            [TDI_adj(i), var_term] = compute_DDI(TDIdata(:,1)', TDIdata(:,2)');
        end %end mag disp
    else
        for j=1:length(unique_mean_disp)
            start = zeros(length(unique_mag_disp), 1);
            stop = zeros(length(unique_mag_disp), 1);
            for k=1:length(unique_mag_disp)
                figure(graph);
                hold on;
                subplot(num_mag_disp*2, num_ap_size,  ((k-1)*(num_ap_size) + i)*2);

                disp_select = logical((ap_size == unique_ap_size(i)) & (mag_disp == unique_mag_disp(k)) & (mean_disp == unique_mean_disp(j)) );
                
                plot_x = disp_ang(disp_select & ~null_trials & ~control_trials & select_trials);
                plot_y = spike_rates(disp_select & ~null_trials & ~control_trials & select_trials);
                
                %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
                [px, py, perr, spk_max, spk_min] = PlotTuningCurve(plot_x', plot_y', symbols{k}, '', 1, 1);
                [TDI(k), var_term] = Compute_DDI(plot_x, plot_y);
                
                %store data to calculate adjusted TDI later
                start(k) = length(TDIdata)+1;
                stop(k) = length(plot_x)+start(k)-1;
                TDIdata(start(k):stop(k), 1) = plot_x';
                TDIdata(start(k):stop(k), 2) = plot_y';
                
                x_interp = (px(1)): 1 : (px(length(px)));
                px = (px * pi)/180;
                plot_x = (plot_x * pi)/180;
                means{k} = [px py];
                raw{k} = [plot_x' plot_y'];
                
                ind_means = [px py];
                ind_raw = [plot_x' plot_y'];
                
                %fit with a distorted sin wave
                ind_pars{k} = sin_exp_fit(ind_means,ind_raw);
                ind_sinerr(k) = sin_exp_err(ind_pars{k});
                
                %grab pref_tilt of indp fits
                x_interp = (px(1)): .01 : (px(length(px)));
                x_deg= x_interp*180/pi;
                y_sin = sin_exp_func(x_interp, ind_pars{k});
                y_err = sin_exp_err(ind_pars{k});
                y_sin(y_sin < 0) = 0;
                %hold on
                %plot(x_deg, y_sin, lines2{k});
                
                %store p-values of each curve
                temp_x =(plot_x *180)/pi;
                p_val(i,k) = calc_mdisp_anovap(disp_select, temp_x, plot_y, unique_disp_ang);
                [value, index_max] = max(y_sin);
                pref_tilt(i,k) = x_deg(index_max);
                
                %run chi^2 test on fit
                [chi2(k), chiP(k)] = Chi2_Test(ind_raw(:,1), ind_raw(:,2), 'sin_exp_func', ind_pars{k}, num_free_params);
                
                x_raw{k} = plot_x;
                y_raw{k} = plot_y;
                
                x_means{k} = px;
                y_means{k} = py;
            end %end mag disp
            %readjust slant responses to fall on the same mean
            %then calc avg TDI
            %shifted_graphs = figure;
            %set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [750 50 500 773], 'Name', 'Mean Adjusted Tilt Tuning Curves');
            total_mean = mean(TDIdata(:,2));
            for count_meandisp = 1:length(unique_mean_disp)
                disp_mean = mean(TDIdata(start(count_meandisp):stop(count_meandisp),2));
                difference = total_mean - disp_mean;
                TDIdata(start(count_meandisp):stop(count_meandisp),2) = TDIdata(start(count_meandisp):stop(count_meandisp),2) + difference;
            end
            [TDI_adj(i), var_term] = compute_DDI(TDIdata(:,1)', TDIdata(:,2)');
        end %end mdisp
    end %end if
    
    
    %fill in starting pars
    avg_freq = 0;
    avg_phase = 0;
    for avg_count = 1:length(ind_pars)
        avg_freq = avg_freq + ind_pars{avg_count}(2);
        avg_phase = avg_phase + ind_pars{avg_count}(3);
    end
    avg_freq = avg_freq/length(ind_pars);
    avg_phase = avg_phase/length(ind_pars);
    
    %amplitude for data1
    q(1) = ind_pars{1}(1);
    
    %frequency
    q(2) = avg_freq;
    
    %phase
    q(3) = avg_phase;
    
    %baseline for data1
    q(4) = ind_pars{1}(4);
    
    %exponential
    q(5) = ind_pars{1}(5);
    
    %fill parameters so they are appropriate for phase slop fxn
    index = 6;
    for temp_counter=2:length(ind_pars)
        q(index) = ind_pars{temp_counter}(1) %amp
        q(index+1) = ind_pars{temp_counter}(4); %baseline
        q(index+2) = ind_pars{temp_counter}(5); %exp
        q(index+3) = ind_pars{temp_counter}(3); %phase slop
        index = index + 4;
    end
    
    %fit with a distorted sin wave with freq shared and phase shared +/- 45 slop
    pars_PhaseSlop = multi_sinfit_with_phase_slop(means, raw, q);    
    
    if(num_ap_size >= num_mag_disp)
        subplot(num_ap_size*2, num_mag_disp,  ((j-1)*(num_mag_disp) + i)*2);
        %plot PhaseSlop data
        hold on;
        pars_slop = pars_PhaseSlop(1:5);
        x_deg = (x_interp * 180)/pi;
        sin_slop{1} = sin_exp_func(x_interp, pars_slop);
        sin_slop{1}(sin_slop{1} < 0) = 0;
        plot(x_deg, sin_slop{1}, lines{1});
        
        %calculate R^2 of mean response
        %add a column of ones to yfit to make regress happy
        y_fit_mean = sin_exp_func(x_means{1}, pars_slop);
        
        %check to see if values are identical
        check = y_fit_mean(1);
        check2 = find(y_fit_mean==check);
        if length(check2)==length(y_fit_mean)
            y_fit_mean(1) = y_fit_mean(1) + (rand * .00001);
        end
        
        y_fit_cell{1} = [ones(length(y_fit_mean),1) y_fit_mean];
        [b_mean, bint_mean, r_mean, rint_mean, stats_mean] = regress(y_means{1}, y_fit_cell{1});
        
        r(1) = stats_mean(1);
        
        %run chi^2 test on fit
        [chi2_simul(1), chiP_simul(1)] = Chi2_Test(raw{1}(:,1), raw{1}(:,2), 'sin_exp_func', pars_slop, num_free_params);
        
        index = 6;
        index_slop = 6;
        for mdisp_ind = 2:length(unique_mean_disp)
            pars_slop(1) = pars_PhaseSlop(index_slop);
            pars_slop(4) = pars_PhaseSlop(index_slop+1);
            pars_slop(5) = pars_PhaseSlop(index_slop+2);
            pars_slop(3) = pars_PhaseSlop(3) + pars_PhaseSlop(index_slop+3);
            index_slop = index_slop + 4;
            hold on;
            sin_slop{mdisp_ind} = sin_exp_func(x_interp, pars_slop);
            sin_slop{mdisp_ind}(sin_slop{mdisp_ind} < 0) = 0;
            plot(x_deg, sin_slop{mdisp_ind}, lines{mdisp_ind});
            
            %calculate R^2 of mean response
            %add a column of ones to yfit to make regress happy
            y_fit_mean = sin_exp_func(x_means{mdisp_ind}, pars_slop);
            
            %check to see if values are identical
            check = y_fit_mean(1);
            check2 = find(y_fit_mean==check);
            if length(check2)==length(y_fit_mean)
                y_fit_mean(1) = y_fit_mean(1) + (rand * .00001);
            end
            
            y_fit_cell{mdisp_ind} = [ones(length(y_fit_mean),1) y_fit_mean];
            [b_mean, bint_mean, r_mean, rint_mean, stats_mean] = regress(y_means{mdisp_ind}, y_fit_cell{mdisp_ind});
            
            r(mdisp_ind) = stats_mean(1);
            
            %run chi^2 test on fit
            [chi2_simul(mdisp_ind), chiP_simul(mdisp_ind)] = Chi2_Test(raw{mdisp_ind}(:,1), raw{mdisp_ind}(:,2), 'sin_exp_func', pars_slop, num_free_params);
            
        end
        
        t=[];
        s=[];
        for lock=1:length(unique_mean_disp)
            s=[s;y_means{lock}];
            t=[t;y_fit_cell{lock}];
        end
        
        [b_all, bint_all, r_all, rint_all, stats_all] = regress(s, t);
        r_all = stats_all(1);
        
        z = sin_exp_func(raw{1}(:,1),pars_PhaseSlop(1:5));
        z(z < 0) = 0;
        error_temp(1) = norm(sqrt(z)-sqrt(raw{1}(:,2)))^2;
        index = 6;
        
        for mdisp_ind = 1:length(unique_mean_disp)
            q_temp = pars_PhaseSlop(1:5);
            
            if mdisp_ind > 1
                q_temp(1) = pars_PhaseSlop(index);
                q_temp(4) = pars_PhaseSlop(index+1);
                q_temp(5) = pars_PhaseSlop(index+2);
                q_temp(3) = pars_PhaseSlop(3) + pars_PhaseSlop(index+3);
                
                index = index + 4;
                z = sin_exp_func(raw{mdisp_ind}(:,1),q_temp);
                z(z < 0) = 0;
                error_temp(mdisp_ind) = norm(sqrt(z)-sqrt(raw{mdisp_ind}(:,2)))^2;
            end
            
            stat_string = sprintf(' %1.3f %1.4f %1.4f %1.4f %1.4f', unique_mean_disp(mdisp_ind), r(mdisp_ind),  TDI(mdisp_ind), chiP_simul(mdisp_ind), r_all, ind_sinerr(mdisp_ind), error_temp(mdisp_ind)); 
            stat_out{mdisp_ind+1} = stat_string;
        end
    elseif(num_ap_size < num_mag_disp)
        subplot(num_mag_disp*2, num_ap_size,  ((j-1)*(num_ap_size) + i)*2);
        hold on;
        pars_short = pars_PhaseSlop(1:5);
        x_deg = (x_interp * 180)/pi;
        sin{1} = sin_exp_func(x_interp, pars_short);
        sin{1}(sin{1} < 0) = 0;
        plot(x_deg, sin{1}, lines{1});
        
        %calculate R^2 of mean response
        %add a column of ones to yfit to make regress happy
        y_fit_mean = sin_exp_func(x_means{1}, pars_short);
        
        %check to see if values are identical
        check = y_fit_mean(1);
        check2 = find(y_fit_mean==check);
        if length(check2)==length(y_fit_mean)
            y_fit_mean(1) = y_fit_mean(1) + (rand * .00001);
        end
        
        y_fit_cell{1} = [ones(length(y_fit_mean),1) y_fit_mean];
        [b_mean, bint_mean, r_mean, rint_mean, stats_mean] = regress(y_means{1}, y_fit_cell{1});
        
        r(1) = stats_mean(1);
        
        %run chi^2 test on fit
        [chi2_simul(1), chiP_simul(1)] = Chi2_Test(raw{1}(:,1), raw{1}(:,2), 'sin_exp_func', pars_short, num_free_params);
        
        index = 6;
        for mag_ind = 2:length(unique_mag_disp)
            subplot(num_mag_disp*2, num_ap_size,  ((mag_ind-1)*(num_ap_size) + i)*2);
            pars_short(1) = pars_PhaseSlop(index);
            pars_short(4) = pars_PhaseSlop(index+1);
            pars_short(5) = pars_PhaseSlop(index+2);
            pars_short(3) = pars_PhaseSlop(3) + pars_PhaseSlop(index+3);
            index = index + 4;
            hold on;
            sin{mag_ind} = sin_exp_func(x_interp, pars_short);
            sin{mag_ind}(sin{mag_ind} < 0) = 0;
            plot(x_deg, sin{mag_ind}, lines{mag_ind});
            
            %calculate R^2 of mean response
            %add a column of ones to yfit to make regress happy
            y_fit_mean = sin_exp_func(x_means{mag_ind}, pars_short);
            
            %check to see if values are identical
            check = y_fit_mean(1);
            check2 = find(y_fit_mean==check);
            if length(check2)==length(y_fit_mean)
                y_fit_mean(1) = y_fit_mean(1) + (rand * .00001);
            end
            
            y_fit_cell{mag_ind} = [ones(length(y_fit_mean),1) y_fit_mean];
            [b_mean, bint_mean, r_mean, rint_mean, stats_mean] = regress(y_means{mag_ind}, y_fit_cell{mag_ind});
            
            r(mag_ind) = stats_mean(1);
    
            %run chi^2 test on fit
            [chi2_simul(mag_ind), chiP_simul(mag_ind)] = Chi2_Test(raw{mag_ind}(:,1), raw{mag_ind}(:,2), 'sin_exp_func', pars_short, num_free_params);            
        end
        
        t=[];
        s=[];
        for lock=1:length(unique_mag_disp)
            s=[s;y_means{lock}];
            t=[t;y_fit_cell{lock}];
        end
        
        [b_all, bint_all, r_all, rint_all, stats_all] = regress(s, t);
        r_all = stats_all(1);

        z = sin_exp_func(raw{1}(:,1),pars_PhaseSlop(1:5));
        z(z < 0) = 0;
        error_temp(1) = norm(sqrt(z)-sqrt(raw{1}(:,2)))^2;
        index = 6;
        for mag_ind = 1:length(unique_mag_disp)
            q_temp = pars_PhaseSlop(1:5);
            
            if mag_ind >= 2
                q_temp(1) = pars_PhaseSlop(index);
                q_temp(4) = pars_PhaseSlop(index+1);
                q_temp(5) = pars_PhaseSlop(index+2);
                q_temp(3) = pars_PhaseSlop(3) + pars_PhaseSlop(index+3);
            
                index = index + 4;
                z = sin_exp_func(raw{mag_ind}(:,1),q_temp);
                z(z < 0) = 0;
                error_temp(mag_ind) = norm(sqrt(z)-sqrt(raw{mag_ind}(:,2)))^2;
            end
            
            stat_string = sprintf(' %1.3f %1.4f %1.4f %1.4f %1.4f', unique_mag_disp(mag_ind), r(mag_ind),  TDI(mag_ind), chiP_simul(mag_ind), r_all, ind_sinerr(mag_ind), error_temp(mag_ind)); 
            stat_out{mag_ind+1} = stat_string;
        end
    end %end big old if statement
    

    
    %print out the parameters for the current mean disparity in the correct subplot
    if(num_ap_size >= num_mag_disp)
        subplot(num_ap_size*2, num_mag_disp,  (((j-1)*(num_mag_disp) + i)*2)-1);
    elseif(num_ap_size < num_mag_disp)
        subplot(num_mag_disp*2, num_ap_size,  (((j-1)*(num_ap_size) + i)*2)-1);
    end
    
    axis([0 100 0 100]);
    axis('off');
    xpos = -10;
    ypos = 110-(20*(i-1));         
    
    line = sprintf('M Disp = %3.2f', unique_mean_disp(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
    line = sprintf('Disp Mag = %3.2f', unique_mag_disp(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
    line = sprintf('Amp(1) = %3.2f', pars_PhaseSlop(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Base(1) = %3.2f', pars_PhaseSlop(4));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Exp(1) = %1.2f', pars_PhaseSlop(5));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Freq = %1.4f', pars_PhaseSlop(2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Phase = %1.4f', pars_PhaseSlop(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    
    line = sprintf('Chi^2 = %1.4f', chi2_simul(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Chi^2 P = %1.4f', chiP_simul(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Ind Err = %3.4f', ind_sinerr(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
  
    line = sprintf('Const Err = %3.4f', error_temp(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;    
    line = sprintf('Overall R = %1.4f', r_all);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;    
    
    index = 6;
    %print out the parameters for the current mean disparity in the correct subplot
    if(num_ap_size >= num_mag_disp)
        subplot(num_ap_size*2, num_mag_disp,  (((j-1)*(num_mag_disp) + i)*2)-1);
        for mdisp_ind = 2:length(unique_mean_disp)      
            if mdisp_ind == 2
                xpos = 15;
            else
                xpos = xpos+18;
            end
            ypos = 110-(20*(i-1));
            line = sprintf('M Disp = %3.2f', unique_mean_disp(mdisp_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
            
            line = sprintf('Amp(%1.0d) = %3.2f', mdisp_ind, pars_PhaseSlop(index));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Base(%1.0d) = %3.2f', mdisp_ind, pars_PhaseSlop(index+1));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
            line = sprintf('Exp(%1.0d) = %3.2f', mdisp_ind, pars_PhaseSlop(index+2));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
            line = sprintf('Phase(%1.0d) = %1.4f', mdisp_ind, pars_PhaseSlop(3) + pars_PhaseSlop(index+3));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            index = index + 4;
            
            line = sprintf('Chi^2 = %1.4f', chi2_simul(mdisp_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Chi^2 P = %1.4f', chiP_simul(mdisp_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Ind Err = %3.4f', ind_sinerr(mdisp_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

            line = sprintf('Const Err = %3.4f', error_temp(mdisp_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;    
        end
    elseif(num_ap_size < num_mag_disp)
        subplot(num_mag_disp*2, num_ap_size,  (((j-1)*(num_ap_size) + i)*2)-1);
        for mag_ind = 2:length(unique_mag_disp)      
            if mag_ind == 2
                xpos = 15;
            else
                xpos = xpos+18;
            end
            ypos = 110-(20*(i-1));
            line = sprintf('M Disp= %3.2f', unique_mean_disp(1));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
            line = sprintf('Disp Mag= %3.2f', unique_mag_disp(mag_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
            line = sprintf('Amp(%1.0d) = %3.2f', mag_ind, pars(index));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Base(%1.0d) = %3.2f', mag_ind, pars(index+1));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
            line = sprintf('Exp(%1.0d) = %3.2f', mag_ind, pars(index+2));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
            
            line = sprintf('Chi^2 = %1.4f', chi2_simul(mag_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Chi^2 P = %1.4f', chiP_simul(mag_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Ind Err = %3.4f', ind_sinerr(mag_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

            line = sprintf('Const Err = %3.4f', error_temp(mag_ind));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;    
        end
    end
    
    
    
       
    %do Sequential F-test
    slop_err = multi_sinerr_PhaseSlop(pars_PhaseSlop);
    independent_err = sum(ind_sinerr);
    
    if(num_ap_size >= num_mag_disp)
        indep_params = length(ind_pars{1})*length(unique_mean_disp);
        paired_slop = indep_params-(length(unique_mean_disp)-1)*2;
    elseif(num_ap_size < num_mag_disp)    
        indep_params = length(ind_pars{1})*length(unique_mag_disp);
        paired_slop = indep_params-(length(unique_mag_disp)-1)*2;        
    end


    Npts_slop = length(spike_rates(~null_trials & ~control_trials & select_trials));
    Fseq_slop = ( (slop_err - independent_err )/(indep_params-paired_slop) ) / ( independent_err/(Npts_slop - indep_params) );
    Pseq_slop = 1 - fcdf(Fseq_slop, (indep_params-paired_slop), (Npts_slop-indep_params) );      

    stat_string = sprintf(' %1.4f', Pseq_slop);     
    for loopcat = 2:length(stat_out)
        stat_out{loopcat} = strcat(stat_out{loopcat}, stat_string);
    end
    
    f_string = sprintf(' %1.4f %1.4f', TDI_adj(i), Pseq_slop); 
    f_out = strcat(f_out, f_string);
    
    line = sprintf('Fseq = %1.4f', Fseq_slop);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Pseq = %1.4f', Pseq_slop);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;   
    line = sprintf('Overall Err = %3.4f', slop_err);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
    
    if(num_ap_size >= num_mag_disp)
        subplot(num_ap_size*2, num_mag_disp,  ((j-1)*(num_mag_disp) + i)*2);
        height = axis;
        yheight = height(4);
        string = sprintf('File = %s', FILE);
        text(height(1)+2, .95*yheight, string, 'FontSize', 8);
        for counter =1:length(unique_mean_disp)
            string = sprintf('r^2 = %1.4f, TDI = %1.4f', r(counter), TDI(counter));
            text_handle = text(height(1)+2, (1-.05*counter-.05)*yheight, string, 'FontSize', 8);
            set(text_handle, 'Color', colors{counter});
        end 
    elseif(num_ap_size < num_mag_disp)
        subplot(num_mag_disp*2, num_ap_size,  ((j-1)*(num_ap_size) + i)*2);
        height = axis;
        yheight = height(4);
        string = sprintf('File = %s', FILE);
        text(height(1)+2, .95*yheight, string, 'FontSize', 8);
        for counter =1:length(unique_mag_disp)
            subplot(num_mag_disp*2, num_ap_size,  ((counter-1)*(num_ap_size) + i)*2);
            string = sprintf('r^2 = %1.4f, TDI = %1.4f', r(counter), TDI(counter));
            text_handle = text(height(1)+2, (1-.05*counter-.05)*yheight, string, 'FontSize', 8);
            set(text_handle, 'Color', colors{counter});
        end         
    end
    

end %end ap size

printme = 0;
if (printme==1)
    PATHOUT = 'Z:\Users\jerry\GradAnalysis\figure_data\';
    
    line = sprintf('%s', FILE);
    for i=2:length(stat_out)
        stat_out{i} = strcat(line, stat_out{i});
    end
        
    %print statistics for each mean disparity
    outfile = [PATHOUT 'simul_r_sqared_mdisp_slop_81902.dat'];
    print_label = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        print_label = 1;
    end
    
    fid = fopen(outfile, 'a');
    if print_label == 1
        fprintf(fid, '%s', [stat_out{1}]);
    end
    for i=2:length(stat_out)
        fprintf(fid, '%s', [stat_out{i}]);
        fprintf(fid, '\r\n');
    end
    fclose(fid);

    f_out = strcat(line, f_out);
    %print F statistics for single cell
    outfile = [PATHOUT 'simul_r_sqared_cell_slop_81902.dat'];
    print_label = 0;
    if (exist(outfile, 'file') == 0)
        print_label = 1;
    end
    fid = fopen(outfile, 'a');
    if print_label == 1    %file does not yet exist
        fprintf(fid, 'File\tTDI\tSeqF_pval\n');
    end
    fprintf(fid, '%s', [f_out]);
    fprintf(fid, '\r\n');
    fclose(fid);
end

