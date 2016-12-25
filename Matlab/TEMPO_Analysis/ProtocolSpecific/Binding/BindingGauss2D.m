function BindingGauss2D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

    output = 1;
    ProtocolDefs;
    symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
    lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};
    
    [conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

    %now, remove trials that do not fall between BegTrial and EndTrial
    trials = 1:size(data.obj_params, 2);		% a vector of trial indices
    select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
    
    %get null trials for spontaneous activity
    null_trials = logical( (data.misc_params(CONDITION,:) == 0) );  		%null condition = condition 0

    auto_channel_select = 0;
    
    if auto_channel_select == 1
%     %select spike chan# for automated batch
        switch SpikeChan
        case 1
            select_trials = select_trials & (data.obj_params(OBJ_STATUS,:,1) == 1) | null_trials;
        case 2
            select_trials = select_trials & (data.obj_params(OBJ_STATUS,:,1) == 0) | null_trials;
        end

        if sum(select_trials == 0)
            %find which channel we are mapping and override current channel
            %being analyzed.  Leave alone if channel correct or mapping both
            select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
            SpikeChan = ~(SpikeChan - 1) + 1;   % switch spike channels
            select_trials = select_trials & (data.obj_params(OBJ_STATUS,:,SpikeChan) == 1) | null_trials;
        end
    end

        %fixation point = 0,0 obj coordinates automatically relative to fix pt.
    fix_x_pos = 0;
    fix_y_pos = 0;

    % get the position values for each condition in the condition_list[]
    x_pos = data.obj_params(OBJ_TRANS_X,:,SpikeChan);
    
    % get the position values for each condition in the condition_list[]
    y_pos = data.obj_params(OBJ_TRANS_Y,:,SpikeChan);
    
    % subtract the fixation coordinates from the x and y position
    x_pos = x_pos - fix_x_pos;
    y_pos = y_pos - fix_y_pos;

    rf_xctr = data.neuron_params(RF_XCTR, SpikeChan);
    rf_yctr = data.neuron_params(RF_YCTR, SpikeChan);
  

    %now get unique values for both x and y position
    unique_x_pos = munique(x_pos(~null_trials & select_trials)');
    unique_y_pos = munique(y_pos(~null_trials & select_trials)');
    
    %now, get the firing rates for all the trials 
    spike_rates = data.spike_rates(SpikeChan, :);

    
    % Calculate spontaneous rates before looping through so can calculate DTI
    null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
    
    Z = zeros(length(unique_x_pos)* length(unique_y_pos), 3);
    display_contours = zeros(length(unique_x_pos), length(unique_y_pos));
            
    for i = 1:length(unique_x_pos)
        for j = 1:length(unique_y_pos)
            indices = logical((x_pos == unique_x_pos(i)) & (y_pos == unique_y_pos(j))  );            
            Z((i-1)*(length(unique_y_pos)) + j, 1) = unique_x_pos(i);
            Z((i-1)*(length(unique_y_pos)) + j, 2) = unique_y_pos(j);
            Z((i-1)*(length(unique_y_pos)) + j, 3) = mean(spike_rates(indices));
        end
    end
    
    for i=1:length(unique_x_pos)
        display_contours(:, i) = Z((i-1)*(length(unique_y_pos))+1:(i-1)*(length(unique_y_pos))+(length(unique_y_pos)), 3)
    end
    
    %display_contours = flipud(display_contours);
    
    curr_figure1 = figure;
    set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', 'Raw Receptive Field Map');
    
    %now, print out some useful information in the upper subplot
    subplot(2, 1, 1);
    PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

    subplot (2,1,2);
    contourf(unique_x_pos, unique_y_pos, display_contours)
    ax1 = axis;
    colorbar
    axis image
    
    raw = [x_pos' y_pos' spike_rates'];
    means = Z;
      
    %fit data here
    pars = gauss2Dfit(means,raw)
    
    %create interpolated arrays for data display
    x_interp = unique_x_pos(1): 1 : unique_x_pos(length(unique_x_pos));
    y_interp = unique_y_pos(1): 1 : unique_y_pos(length(unique_y_pos));
    z_gauss = zeros(length(x_interp), length(y_interp));
      
    %obtain fitted data for interpolated arrays
    for i=1:length(x_interp)
        for j = 1:length(y_interp)
            z_gauss(i,j) =  gauss2Dfunc(x_interp(i),y_interp(j), pars);
        end
    end
    
    z_gauss = rot90(z_gauss);
    z_gauss = rot90(z_gauss);
    z_gauss = rot90(z_gauss);
    z_gauss = fliplr(z_gauss);
    
        subplot (2,1,1);

    %print out the parameters of the data 
	axis([0 100 0 100]);
	axis('off');
    font_size = 9;
    bump_size = 7;
	
    xpos = -10;   
    ypos = 20;
   
    line = sprintf('Raw RF Map:');
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
    line = sprintf('Spontaneous Rate = %1.3f', null_rate);
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      

    
    
    
    curr_figure2 = figure;
    %now, print out some useful information in the upper subplot
    subplot(2, 1, 1);
    PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', '2D-Gaussian-Fitted Receptive Field Map');
    subplot (2,1,2);
    contourf(x_interp, y_interp, z_gauss)
    colorbar
    axis image
       
    subplot (2,1,1);

    %print out the parameters of the fit 
	axis([0 100 0 100]);
	axis('off');
    font_size = 9;
    bump_size = 7;
	
    xpos = -10;   
    ypos = 20;
   
    base_rate = pars(1);
    amplitude = pars(2);
    x_diam = pars(4) * 2;
    y_diam = pars(6) * 2;
    fit_xctr = pars(3);
    fit_yctr = pars(5);
    
    line = sprintf('Fitted RF:');
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
    line = sprintf('Base Rate = %1.3f', base_rate);
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
    line = sprintf('Amplitude = %1.3f', amplitude);
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
    line = sprintf('X Diameter = %1.3f', x_diam);
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
    line = sprintf('Y Diameter = %1.3f', y_diam);
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
    line = sprintf('Fitted CTR = (%1.3f, %1.3f)', fit_xctr, fit_yctr);
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
    line = sprintf('RFPLOT CTR = (%1.3f, %1.3f)', rf_xctr, rf_yctr);

        %output tuning curve metrics
    if (output == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\RF_Mapping\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'map'];
        eval(['save ' PATHOUT FILEOUT  ' base_rate amplitude x_diam y_diam fit_xctr fit_yctr rf_xctr rf_yctr SpikeChan StartCode StopCode BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
        
        %save figure
        FILEOUT = [FILE(1:i-1) '_raw_' num2str(SpikeChan) '.fig'];        
        saveas (curr_figure1, [PATHOUT FILEOUT]);
        FILEOUT = [FILE(1:i-1) '_fit_' num2str(SpikeChan) '.fig'];        
        saveas (curr_figure2, [PATHOUT FILEOUT]);
    end
