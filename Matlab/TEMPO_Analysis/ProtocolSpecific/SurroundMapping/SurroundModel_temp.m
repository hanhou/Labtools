function SurroundModel(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
   TEMPO_Defs;
   
   symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
   lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};
   
   %get the column of disparity values for the surround patches
   disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);
   
	%get indices of any NULL conditions (for measuring spontaneous activity)
   null_trials = logical( (disp == data.one_time_params(NULL_VALUE)) );
   patch_off_condition = logical((disp == data.one_time_params(PATCH_OFF)));
   
   unique_disp = munique(disp(~null_trials)');
   unique_disp = unique_disp(2:length(unique_disp));
   
   %get the column of values of offset angle of the surround patches
   ang = data.dots_params(DOTS_AP_OFF_ANG,BegTrial:EndTrial,PATCH1);
   unique_ang = munique(ang(~null_trials)');
   unique_ang = unique_ang(2:length(unique_ang));
   
   %get the column of different aperture sizes
   ap_size = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH1);
   surround_size = munique(ap_size(~null_trials)');
   surround_size = surround_size(2:length(surround_size));
   
   %now, get the firing rates for all the trials 
   spike_rates = data.spike_rates(SpikeChan, :);
   
   %get indices of monoc. and uncorrelated controls
   control_trials = logical( (disp == LEYE_CONTROL) | (disp == REYE_CONTROL) | (disp == UNCORR_CONTROL) );
   
   %now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:length(disp);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   disp_tuning_curves = zeros(length(unique_disp), length(unique_ang));

   %grab tuning curves for each of the surround patches
   for i = 1:length(unique_ang)
        hold on;
   end
   
   res = 0.5;
   
   %now create a grid based upon the size of the surround patches
   %fill that grid with a gaussian
   Surround_x_grid = -surround_size/2:res:surround_size/2;
   Surround_y_grid = -surround_size/2:res:surround_size/2;
   
   %calculate the parameters for a 2D gaussian, where grid_size/2 equals the center of the grid
   Surround_gauss_x = surround_size/2;
   Surround_gauss_y = surround_size/2;


   a = sqrt(-(Surround_gauss_x^2)/log(.1))*2;
   b = sqrt(-(Surround_gauss_y^2)/log(.1))*2;
   
   Ksurr = -.2;
   
   for i = 1:length(unique_ang)
      Positions(i+1).Gauss_grid = zeros(length(Surround_x_grid), length(Surround_y_grid));
      for j = 1:length(Positions(i+1).Gauss_grid)
         %get the x coordinate
         Surround_x_ind = Surround_x_grid(j);
         for k = 1:length(Positions(i+1).Gauss_grid)
            %get the y coordinate
            Surround_y_ind = Surround_y_grid(k);
            %calculate the gaussian value for this point
            Positions(i+1).Gauss_grid(j, k) = Ksurr*(exp(-((Surround_x_ind^2)/a + (Surround_y_ind^2)/b)));
        end
    end
   end

   %now create a grid that is the response of the surround
   %this is done by fitting the surround disparity tuning data with a gabor
   %then using the resulting function to create a response grid
  
   gaborfig = figure;
   tuningfig = figure;
   for i = 1:length(unique_ang)
      ang_select = logical(ang == unique_ang(i));

      plot_x = disp(ang_select & ~null_trials & ~patch_off_condition & ~control_trials & select_trials);
      plot_y = spike_rates(ang_select & ~null_trials & ~patch_off_condition & ~control_trials & select_trials);
      
      %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
      figure(tuningfig)
      hold on
      [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
      
      means = [px py];
      raw = [plot_x' plot_y'];
      [pars,freq] = gaborfit(means,raw);
      x_interp = (px(1)): .01 : (px(length(px)));
      y_interp =  pars(1) + pars(2)*exp(-0.5*((x_interp - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(x_interp - pars(3))+pars(6) );
      y_interp(y_interp < 0) = 0;
      
      figure(gaborfig)
      hold on
      plot(x_interp, y_interp,lines{i})
      
      % Note: this func MUST match that in gaborfit.m
      %x_gauss = (px(1)): .01 : (px(length(px)));
      %y_gauss =  pars(1) + pars(2)*exp(-0.5*((x_gauss - pars(3))/ pars(4)).^2);
      %x_sine = (px(1)): .01 : (px(length(px)));
      %y_sine =  pars(1) + 0.5*pars(2)*cos(2*pi*pars(5)*(x_sine - pars(3))+pars(6) );
   end
   
  
   %multiply the gabor response grid and the gaussian grid to create 
   %the response of the neuron
   