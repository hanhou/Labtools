function HGradModel(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
global distance_grid disparity_grid response_grid gauss_grid response_grid_after_gaussian tot_response

DEBUG = 0;

TEMPO_Defs;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--'};

%get position of peak response within the size tuning curve
%use this as the size of the receptive field (will be updated later for more sophisticated
%size analyses.

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);

%get the column of values of horiz. disparity magnitude in the dots_params matrix
mag_disp = data.dots_params(DOTS_HGRAD_MAG,BegTrial:EndTrial,PATCH1);
   
%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical((mag_disp == data.one_time_params(NULL_VALUE)) );
   
unique_mag_disp = munique(mag_disp(~null_trials)');	
   
%get the column of values of horiz. disparity angle of orientation in the dots_params matrix
disp_ang = data.dots_params(DOTS_HGRAD_ANGLE,BegTrial:EndTrial,PATCH1);
unique_disp_ang = munique(disp_ang(~null_trials)');

%get the column of mean disparity values
mean_disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);
unique_mean_disp = munique(mean_disp(~null_trials)');

%get the column of values of directions in the dots_params matrix
size = data.dots_params(DOTS_AP_XSIZ,:,PATCH1);
   
%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (size == data.one_time_params(NULL_VALUE)) );
   
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(size);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%plot a size tuning curve to determine what the CRF size will be
plot_x = size(~null_trials & select_trials);
plot_y = spike_rates(~null_trials & select_trials);
   
%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
[px, py, perr] = PlotTuningCurve(plot_x', plot_y', 'bo', 'b-', 0, 0);

size_range = [px py];
max_size_vector = sortrows(size_range);

%this will be the diameter of the grid
grid_size = max_size_vector(1);

%the resolution of the grid
res = .5; %degrees

x_grid = -grid_size/2:res:grid_size/2;
y_grid = -grid_size/2:res:grid_size/2;

%if(DEBUG)
%   unique_mean_disp = [-0.4 0 0.4]';
%end

%print out a legend graph that lists each mean disparity and the line type with which they are drawn
leg_x = [0 1];
leg_y = [0 0];
legend_figure = figure;
set(legend_figure, 'Position', [750 873 100 100], 'Name', 'Line Legend', 'MenuBar', 'none');
for i = 1:length(unique_mean_disp)
   axis([0 1 0 4]);
   axis('off')
   hold on
   leg_y = leg_y + 1;
   plot(leg_x, leg_y, char(lines(i)));
	string = sprintf('M. Disp = %1.3f', unique_mean_disp(i));
   text(0, leg_y(1)+.25, string, 'FontSize', 8);
   hold off
end

if(DEBUG)
    %define figure handle for the offsetgraph
	offsetgraph = figure;
   plot(0,0,'*');
end

if(DEBUG)
    %define figure handle for the gauss graph
    gauss_fig = figure;
end

if(1)
    %create a figure handle for every unique disparity magnitude
   for k = 1:length(unique_mag_disp)
      plot_figures(k) = figure;
      string = sprintf('Mag Disp = %4.1f', unique_mag_disp(k));
      set(gcf, 'Name', string);
   end
end

%these two main loops specify the offsets that will be simulated
offset = [1]; %in terms of degrees  
					 

for num_off = 1:length(offset)
    for offloopx = -1:1
        xoffset = (offset(num_off)/res) * offloopx;
        for offloopy = -1:1
	        if(DEBUG)
    		   figure(offsetgraph)
    		   hold on
    		   plot(offloopx, offloopy, '*');
    		   hold off
            end

    		yoffset = (offset(num_off)/res) * offloopy;

        	%load disparity tuning parameters for this cell
      	    load disptuning2
    	    sigma = perr;
    	    mu = py;
      
            %add in extra values on both ends of the tuning curve
            if(1)
                uncorr_resp = save_uncorr_resp;
		        for i = 1:3
                    new_x = min(px(:)) - 0.2;
                    px = [new_x; px];
                    py = [uncorr_resp; py];
                    sigma = [0;sigma];
                    mu = [0;mu];
                end
      
                for i = 1:3
                    new_x = max(px(:)) + 0.2;
                    px = [px; new_x ];
                    py = [py; uncorr_resp];
                    sigma = [sigma;0];
                    mu = [mu;0];
                end
            end
   
            pystart = py;
   
            %alter disparity tuning values here for debugging
            if (DEBUG)
    			unique_mean_disp = [-0.4 0 0.4]';
			    px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
                pystart = [100 100 100 120 130 150 130 120 100 100 100]';
            end
      
            %the (x,y) grid coordinates for the various graphs
            x_grid = -grid_size/2:res:grid_size/2;
            y_grid = -grid_size/2:res:grid_size/2;
      
            %the (x,y) coordinates of the center 
            %will be translated from grid coordinates into
            %array indices
            xctr = round(length(x_grid)/2)+xoffset;
            xctr = x_grid(xctr);

            yctr = round(length(y_grid)/2)+yoffset;
            yctr = y_grid(yctr);

            %calculate the parameters for a 2D gaussian, where grid_size/2 equals the center of the grid
            gauss_x = grid_size/2;
            gauss_y = grid_size/2;
            a = sqrt(-(gauss_x^2)/log(.1));
            b = sqrt(-(gauss_y^2)/log(.1));

            %initialize a matrix to hold the distance of the stimulus from the CRF center
            distance_grid = zeros(length(x_grid));
            gauss_grid = zeros(length(x_grid));
      
            %create a grid that contains the distance of each cell from the center of the receptive field
            %create a matrix that contains that values of a gaussian distribution centered
            %on (0,0) coordinates
    		for i = 1:length(distance_grid)
                %get the x coordinate
                x_ind = x_grid(i);
                for j = 1:length(distance_grid)
                    y_ind = y_grid(j);
                    %get the y coordinate
                    distance_grid(i, j) = sqrt((x_ind-xctr)^2 + (y_ind-yctr)^2);
	      
                    if distance_grid(i,j) > grid_size/2
                        distance_grid(i,j) = NaN;
                    end      
	      
                    %calculate the gaussian value for this point
                    gauss_grid(i,j) = exp(-((x_ind^2)/a + (y_ind^2)/b));
                end
            end
	
            if(DEBUG)
                figure(gauss_fig)
                %hold on;
                %surf(distance_grid);
                surf(gauss_grid)
                view(0, 90);
                hold off;
            end

            tot_response = {};
      
            %make a new set of (x,y) coordinates that are shifted to reflect
            %the offset of the experimental receptive field
            x_grid = (-grid_size/2)-xctr:res:(grid_size/2)-xctr;
            y_grid = (-grid_size/2)-yctr:res:(grid_size/2)-yctr;

            if(DEBUG)
                unique_mean_disp = unique_mean_disp(1);
                unique_disp_ang = unique_disp_ang(1:4);
                unique_mag_disp = unique_mag_disp(1);
            end
   
            %number of times to repeat the tuning curve
            num_reps = 5;

            for i=1:length(unique_mean_disp)  %for each mean disparity...
                for j =1:length(unique_mag_disp) %for each magnitude of disparity
                    for k =1:length(unique_disp_ang)
                        angle = unique_disp_ang(k) * (pi/180);
           
                        for reps = 1:num_reps
               
                            %create a matrix where each cell contains the disparity of that location
                            disparity_grid = zeros(length(x_grid)) + unique_mean_disp(i);
                            response_grid = zeros(length(x_grid));
                            %add in the error for the spikes here
                            randn('state',sum(100*clock))
                            perr = randn(1, length(pystart));
                            perr = perr' .* sigma;
                            py = perr + mu;
                            
                            for l = 1:length(disparity_grid) %for each x pt on the grid
                                %get the x coordinate of the grid
                                x_ind = x_grid(l);
                                for m = 1:length(disparity_grid) %for each y point on the grid
                                    %get the y coordinate of the grid
                                    y_ind = y_grid(m);
                                    if ~isnan(distance_grid(l,m))
                                        if(angle==0)
                                            disp_val = y_ind * unique_mag_disp(j);
                                            %is alpha != 0 then there will be a tilt to the gradient
                                        else
                                            len = distance_grid(l,m);
                                            if len ~= 0
                                                theta = asin(y_ind/len);
                                                if ~isreal(theta)
                                                    theta = real(theta);
                                                end
                                                %//sin only goes from -90 to 90 degrees.  need to fix for the other two quadrants.
                                                %//need to do quadrant corrections
    											%//cases (-, +), (-, -)
    											if ((x_ind <= 0) & (y_ind >= 0))
    											%//Quadrant III
    											%//theta is now greater than 90 degrees
                                                    temptheta = (pi/2)-theta;
    												theta = (pi/2) + temptheta;
    											elseif ((x_ind <= 0) & (y_ind <= 0))
    											%//Quadrant IV
    											%//theta is now less than -90 degrees
    												temptheta = (-pi/2)-theta;
                                                    theta = (-pi/2) + temptheta;
    							                end %end quad check if
                	
                                                %//calculate the distance of the point from the gradient tilt axis
                                                y = len * sin(theta-angle);
    		
                                                %//calculate the disparity difference
                                                disp_val = y * unique_mag_disp(j);
                                            else
                                                disp_val = 0;
                                            end %end len if
                                        end %end angle if
                                        disparity_grid(l,m) = disparity_grid(l,m) + disp_val;
                                        %using the values of the tuning curve, interpolate to find what the response of the cell
                                        %would have been to the disparity found in each cell
                                        response_grid(l,m) = interp1(px, py, disparity_grid(l,m));
                                    else %else isnan if
                                        disparity_grid(l,m) = NaN;
                                        response_grid(l,m) = 0;
                                    end %end isnan if
                                end %end m
                            end %end l
         
                            %the disparity matrix is now filled
                            %then multiply the response curve to get the response of the cell for each gradient
                            response_grid_after_gaussian = response_grid .* gauss_grid;
                            response(k,reps) = sum(response_grid_after_gaussian(:));
               
                            if(DEBUG)
                                [max_resp, max_i] = max(response_grid(:));
                                max_resp_gauss = max(response_grid_after_gaussian(:));
                                pref_disp = disparity_grid(max_i);
                                max_disp = max(disparity_grid(:));
                                min_disp = min(disparity_grid(:));
                            end
                        
                            if(DEBUG)
                                string = sprintf('Disparities: MN DSP = %4.1f, DSP Ang = %4.1f, Mag DSP = %4.1f', unique_mean_disp(i), unique_disp_ang(k), unique_mag_disp(j));
                                disp_fig = figure;
                                surf(disparity_grid)
                                set(disp_fig, 'Name', string);
                                view(0, 90);
                            end
         
                            if(DEBUG)
                                height = axis;
                                yheight = height(4);
                                string = sprintf('Summed Response = %4.1f', sum(response_grid_after_gaussian(:)));
                                text(height(1)+2, 0.95*yheight, string, 'FontSize', 8);
                                string = sprintf('Max Response = %4.1f', max_resp_gauss);
                                text(height(1)+2, 0.90*yheight, string, 'FontSize', 8);
                                string = sprintf('Max Disp = %4.1f', max_disp);
                                text(height(1)+2, 0.85*yheight, string, 'FontSize', 8);
                                string = sprintf('Min Disp = %4.1f', min_disp);
                                text(height(1)+2, 0.8*yheight, string, 'FontSize', 8);
                                string = sprintf('Pref Disp = %4.1f', pref_disp);
                                text(height(1)+2, 0.75*yheight, string, 'FontSize', 8);
                            end
                  
                            if(DEBUG)
                                resp_fig = figure;
                                surf(response_grid)
                                string = sprintf('Raw Response: MN DSP = %4.1f, DSP Ang = %4.1f, Mag DSP = %4.1f', unique_mean_disp(i), unique_disp_ang(k), unique_mag_disp(j));
                                set(resp_fig, 'Name', string);
                            end
      	   
                            if(DEBUG)
                                resp_fig_gauss = figure;
                                surf(response_grid_after_gaussian)
                                string = sprintf('Response After Gaussian: MN DSP = %4.1f, DSP Ang = %4.1f, Mag DSP = %4.1f', unique_mean_disp(i), unique_disp_ang(k), unique_mag_disp(j));
                                set(resp_fig_gauss, 'Name', string);
                            end
               
                        end %end reps loop
                    end %end (k) disp_ang loop
                    temp_resp(:, j) = response(:);
                    %temp_resp{j} = response;
                end %end (j) mag_disp loop
                tot_response{i} = temp_resp;
            end %end (i) mean_disp loop
   
            if(DEBUG)
                celldisp(tot_response);
            end
   
            %normalization of the values
            for i = 1:length(tot_response)
                temp = tot_response{i};
                maxtemp(i) = max(temp(:));
            end

            maxval = max(maxtemp(:));
   
            for i = 1:length(tot_response)
                temp = tot_response{i};
                tot_response{i} = temp/maxval;
            end

            if(DEBUG)
                celldisp(tot_response);
            end

            %plot the graph values
            if(DEBUG)
                graph = figure;
                string = sprintf('H Disp Angle of Rotation Tuning Curve (Xoff = %2d, Yoff = %2d)', offloopx, offloopy);
                set(graph,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 487 500 873/2], 'Name', string);
            end
   
            for i=1:length(unique_mag_disp)
                if(DEBUG)
                    subplot(length(unique_mag_disp), 1, i);
                end
      
                if(1)
                    figure(plot_figures(i));
                end
      
                if(1)
                    if(offloopx == -1)
                        if(offloopy == -1)
                            plot_here = 7;
                        elseif(offloopy == 0)
                            plot_here = 4;
                        elseif (offloopy == 1)
                            plot_here = 1;
                        end
                    elseif(offloopx == 0)
                        if(offloopy == -1)
                            plot_here = 8;
                        elseif(offloopy == 0)
                            plot_here = 5;
                        elseif (offloopy == 1)
                            plot_here = 2;
                        end
                    elseif(offloopx == 1)
                        if(offloopy == -1)
                            plot_here = 9;
                        elseif(offloopy == 0)
                            plot_here = 6;
                        elseif (offloopy == 1)
                            plot_here = 3;
                        end
                    end
         
                    subplot(3, 3, plot_here);
                end
      
                for j=1:length(tot_response)
                    temp = tot_response{j};
                    hold on
            
                    plot_x=[];
                    for ij=1:num_reps
                        plot_x = [plot_x; unique_disp_ang];
                    end
                    plot_y = temp(:,i);
                  	      
                    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
                    [px, py, perr] = PlotTuningCurve(plot_x, plot_y, symbols{j}, lines{j}, 1, 1);
                    hold off;
         
                    hold on;
                    errorbar(px, py, perr, perr, symbols{j});
                    hold off;
            
                    spike(:,1) = plot_x;
                    spike(:,2) = plot_y;
                    sortedspikes = sortrows(spike, [1]);
                    for ik = 1:length(unique_disp_ang)
                        spikes_arranged_by_angles(:,ik) = sortedspikes((num_reps*(ik-1)+1):ik*num_reps,2);
                    end
                    spikes_by_angles_by_mean_disp((num_reps*(j-1)+1):j*num_reps, :) = spikes_arranged_by_angles;
            
                end %end j loop tot_response loop
                p(i,:) = anova2(spikes_by_angles_by_mean_disp, num_reps);
                close(gcf);
                figure(plot_figures(i));
                height = axis;
                string = sprintf('P Values = %1.4f %1.4f %1.4f', p(i,:));
                text(height(1)+2, 1, string, 'FontSize', 8);
         
            end %end i loop mag_disp loop
      
        end %end offloopy
    end %end offloopx
end %end num_off loop