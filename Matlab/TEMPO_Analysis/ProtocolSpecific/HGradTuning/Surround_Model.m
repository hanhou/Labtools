function Surround_Model(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);

global done;

%first thing to do is to create a series of matrices for each surround area.
%1) gaussian receptive field 
%2) response of the surround using a the inverted version of the tuning curve of the CRF
%3) multiply the response by the gaussian
%4) now create a matrix large enough to hold the CRF in the center and have space for the surround fields
%5) add the response matrices of each of the surrounds and the CRF into the large matrix
%6) to calculate the entire spike rate of the neuron for this stimulus, make all negative values zero
%   and then add together all the cells in the matrix.


TEMPO_Defs;

%for h. disp i need to obtain the following info:
   %perr
   %px
   %py
   %uncorr_resp
   %plot_x
   %plot_y

%for surround mapping i need to obtain the following info:
   %surround size
   %center size

   %and for each surround patch:
      %perr
      %px
      %py
      %uncorr_resp
      %plot_x
      %plot_y

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
CRF_Size = max_size_vector(1);

Total_RF_Size = max_size_vector(2)
Surround_Size = (Total_RF_Size-CRF_Size)/2;

%caclulate the distance between the centers of the surround and the CRF
Distance_CRF_Sur_ctrs = round(Surround_Size/2+CRF_Size/2)


%resolution of matrices (in degrees)
res = .5;

%this loads the disparity tuning of the CRF
load disptuning1
sigma = perr;
mu = py;

uncorr_resp = save_uncorr_resp;
for i = 1:5
    new_x = min(px(:)) - 0.2;
    px = [new_x; px];
    py = [uncorr_resp; py];
    sigma = [0;sigma];
    mu = [0;mu];
end
      
for i = 1:5
    new_x = max(px(:)) + 0.2;
    px = [px; new_x ];
    py = [py; uncorr_resp];
    sigma = [sigma;0];
    mu = [mu;0];
end

CRF_x_grid = -CRF_Size/2:res:CRF_Size/2;
CRF_y_grid = -CRF_Size/2:res:CRF_Size/2;

%calculate the parameters for a 2D gaussian, where grid_size/2 equals the center of the grid
CRF_gauss_x = CRF_Size/2;
CRF_gauss_y = CRF_Size/2;

Kctr = 1.0;
a = sqrt(-(CRF_gauss_x^2)/log(.1));
b = sqrt(-(CRF_gauss_y^2)/log(.1));

%initialize a matrix to hold the distance of the stimulus from the CRF center
Disparity(1).Gauss_grid = zeros(length(CRF_x_grid));
      
%create a grid that contains the distance from the center of the receptive field
%create a matrix that contains that values of a gaussian distribution centered
%on (0,0) coordinates
for i = 1:length(Disparity(1).Gauss_grid)
    %get the x coordinate
    CRF_x_ind = CRF_x_grid(i);
    for j = 1:length(Disparity(1).Gauss_grid)
        %get the y coordinate
        CRF_y_ind = CRF_y_grid(j);
        %calculate the gaussian value for this point
        Disparity(1).Gauss_grid(i,j) = Kctr*exp(-((CRF_x_ind^2)/a + (CRF_y_ind^2)/b));
    end
end

Surround_x_grid = -Surround_Size/2:res:Surround_Size/2;
Surround_y_grid = -Surround_Size/2:res:Surround_Size/2;

%calculate the parameters for a 2D gaussian, where grid_size/2 equals the center of the grid
Surround_gauss_x = Surround_Size/2;
Surround_gauss_y = Surround_Size/2;


a = sqrt(-(Surround_gauss_x^2)/log(.1))*2;
b = sqrt(-(Surround_gauss_y^2)/log(.1))*2;

num_surround = 6;
%weights of the surround fields 
%can adjust to make a spatially heterogeneous surround
%1 - left field
%2 - bottom left
%3 - bottom right
%4 - right
%5 - top right
%6 - top left

                        %1 - right field     %2 - bottom right    %3 - bottom left   %4 - left  %5 - top left  %6 - top right
Ksurr(1:num_surround) = [-.2                   -.2                      -.2               -.2        -.2            -.2];
%initialize a matrix to hold the distance of the stimulus from the CRF center

      
%create a grid that contains the distance from the center of the receptive field
%create a matrix that contains that values of a gaussian distribution centered
%on (0,0) coordinates

for i = 1:num_surround
    Disparity(i+1).Gauss_grid = zeros(length(Surround_x_grid), length(Surround_y_grid));
    for j = 1:length(Disparity(i+1).Gauss_grid)
        %get the x coordinate
        Surround_x_ind = Surround_x_grid(j);
        for k = 1:length(Disparity(i+1).Gauss_grid)
            %get the y coordinate
            Surround_y_ind = Surround_y_grid(k);
            %calculate the gaussian value for this point
            Disparity(i+1).Gauss_grid(j, k) = Ksurr(i)*(exp(-((Surround_x_ind^2)/a + (Surround_y_ind^2)/b)));
        end
    end
end

%fit a gabor function to the CRF disparity responses
%	pars(1) is base rate
%   pars(2) is amplitude
%   pars(3) is center
%   pars(4) is size (sigma)
%	pars(5) is the frequency
%   pars(6) is phase of the sinusoid

load plot_values2
means = [px py];
raw = [plot_x' plot_y'];
[pars,freq] = gaborfit(means,raw);
x_interp = (px(1)): .01 : (px(length(px)));
y_interp =  pars(1) + pars(2)*exp(-0.5*((x_interp - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(x_interp - pars(3))+pars(6) );
y_interp(y_interp < 0) = 0;
% Note: this func MUST match that in gaborfit.m
x_gauss = (px(1)): .01 : (px(length(px)));
y_gauss =  pars(1) + pars(2)*exp(-0.5*((x_gauss - pars(3))/ pars(4)).^2);
x_sine = (px(1)): .01 : (px(length(px)));
y_sine =  pars(1) + 0.5*pars(2)*cos(2*pi*pars(5)*(x_sine - pars(3))+pars(6) );

no_shift = pars(3);
shift_far = pars(3) + 1;
shift_near = pars(3) - 1;

widen = 2*pars(4);
phase_shift = 0;
lower_freq = pars(5)/4;

%move baseline to 0;
pars(1) = 0;

y_surr=zeros(num_surround, length(y_interp));
for i=1:num_surround
    pars(4) = widen;
    pars(5) = lower_freq;
    pars(6) = phase_shift
    if (i==2) | (i == 3)
        pars(3) = shift_near;
        y_interp2 =  (pars(1) + pars(2)*exp(-0.5*((x_interp - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(x_interp - pars(3))+pars(6) ));
    elseif (i==5) | (i==6)
        pars(3) = shift_far;
        y_interp2 =  (pars(1) + pars(2)*exp(-0.5*((x_interp - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(x_interp - pars(3))+pars(6) ));
    else
        pars(3) = no_shift;
        y_interp2 =  (pars(1) + pars(2)*exp(-0.5*((x_interp - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(x_interp - pars(3))+pars(6) ));
    end
            
    y_surr(i,:)=y_interp2;
end



figure
hold on
plot(x_interp, y_interp, 'k--');
for i=1:num_surround
    plot(x_interp, y_surr(i,:), lines{i});
end
%plot(x_gauss, y_gauss, 'm--');
%plot(x_sine, y_sine, 'g:');     
hold off


Sur_size = length(Disparity(2).Gauss_grid);
CRF_size = length(Disparity(1).Gauss_grid);
length_grid = (CRF_size + (2*Sur_size));

%create a matrix that will contain the CRF and all surround fields
total_RF = zeros(length_grid);

%add CRF to the center of the RF
total_RF_x_ctr = round(length(total_RF)/2);
total_RF_y_ctr = round(length(total_RF)/2);

CRF_x_start = total_RF_x_ctr-round(CRF_size/2);
CRF_y_start = total_RF_y_ctr-round(CRF_size/2);
total_RF(CRF_x_start:CRF_x_start-1+CRF_size, CRF_y_start:CRF_y_start-1+CRF_size) = Disparity(1).Gauss_grid;

Disparity(1).px = x_interp;
Disparity(1).py = y_interp;
Disparity(1).present = zeros(length_grid);
Disparity(1).present(CRF_x_start:CRF_x_start-1+CRF_size, CRF_y_start:CRF_y_start-1+CRF_size) = 1;
Disparity(1).Gauss_in_RF = zeros(length_grid);
Disparity(1).Gauss_in_RF(CRF_x_start:CRF_x_start-1+CRF_size, CRF_y_start:CRF_y_start-1+CRF_size) = Disparity(1).Gauss_grid;

shift_val = round(.1*Sur_size);

sur_ang = [0 60 120 180 240 300];
for i = 2:num_surround+1
    sur_ang_rad = sur_ang(i-1) * (pi/180);
    Sur_ctr_x = Distance_CRF_Sur_ctrs * cos(sur_ang_rad) + total_RF_x_ctr;
    Sur_ctr_y = Distance_CRF_Sur_ctrs * sin(sur_ang_rad) + total_RF_y_ctr;
    
    Sur_x_start = round(Sur_ctr_x - Sur_size/2);
    Sur_y_start = round(Sur_ctr_y - Sur_size/2);

    total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size)+Disparity(i).Gauss_grid;
    
    Disparity(i).px = x_interp;
    Disparity(i).py = y_surr(i-1, :);
    Disparity(i).present = zeros(length_grid);
    Disparity(i).present(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = 1;
    Disparity(i).Gauss_in_RF = zeros(length_grid);
    Disparity(i).Gauss_in_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = Disparity(i).Gauss_grid;
    
end


figure
surf(total_RF)

distance_grid = total_RF;

x = zeros(length(distance_grid), 1);
y = zeros(length(distance_grid), 1);


%the next four loops change the indices of total_RF into meaningful coordinates in degrees
j = 1;
for i = total_RF_x_ctr-1:-1:1;
    x(i) = -(res*j);
    j = j + 1;
end

j = 1;
for i = total_RF_x_ctr+1:length(x)
    x(i) = (res*j);
    j = j+1;
end

j = 1;
for i = total_RF_y_ctr-1:-1:1;
    y(i) = (res*j);
    j = j + 1;
end

j = 1;
for i = total_RF_y_ctr+1:length(y)
    y(i) = -(res*j);
    j = j+1;
end



for i = 1:length(x)
    for j = 1:length(y)
        distance_grid(i, j) = sqrt((x(i)-0)^2 + (y(j)-0)^2);
        if distance_grid(i,j) > Total_RF_Size/2
            distance_grid(i,j) = NaN;
        end

    end
end

%debugging
%unique_disp_ang = [0 90 180 270]';

for mdisp=1:length(unique_mean_disp)  %for each mean disparity...
    resp_mag_ang = zeros(length(unique_mag_disp), length(unique_disp_ang));
    for magdisp =1:length(unique_mag_disp) %for each magnitude of disparity
        for dispang =1:length(unique_disp_ang)
            %create a matrix where each cell contains the disparity of that location
            disparity_grid = zeros(length(x)) + unique_mean_disp(mdisp);
            response_grid = zeros(length(x));

            for each_field=1:num_surround+1
                Disparity(each_field).response = zeros(length_grid);
            end

            %add in the error for the spikes here
            randn('state',sum(100*clock))
            perr = randn(1, length(py));
            perr = perr' .* sigma;
            py = perr + mu;

            angle = unique_disp_ang(dispang) * (pi/180);

            for l = 1:length(disparity_grid) %for each x pt on the grid
                %get the x coordinate of the grid
                x_ind = x(l);
                for m = 1:length(disparity_grid) %for each y point on the grid
                    %get the y coordinate of the grid
                    y_ind = y(m);
                    if ~isnan(distance_grid(l,m))
                        if(angle==0)
                            disp_val = y_ind * unique_mag_disp(magdisp);
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
                                y_tilt = len * sin(theta-angle);
                		
                                %//calculate the disparity difference
                                disp_val = y_tilt * unique_mag_disp(magdisp);
                            else
                                disp_val = 0;
                            end %end len if
                        end %end angle if
                        disparity_grid(l,m) = disparity_grid(l,m) + disp_val;
                        %using the values of the tuning curves for each field, interpolate to find what the response of the field
                        %would have been to the disparity found in each cell of the matrix
                        for each_field=1:num_surround+1
                            if Disparity(each_field).present(l,m) == 1
                                Disparity(each_field).response(l, m) = interp1(Disparity(each_field).px, Disparity(each_field).py, disparity_grid(l,m));
                            else
                                Disparity(each_field).response(l,m) = 0;
                            end
                            
                        end
                    else %else isnan if
                        for each_field=1:num_surround+1
                            Disparity(each_field).response(l, m) = 0;
                            disparity_grid(l,m) = NaN;
                        end
                    end %end isnan if
                end %end m
            end %end l

            %now multiply the response of each field by their respective gaussians
            for each_field = 1:num_surround+1
                Disparity(each_field).res_mult_gauss = zeros(length_grid);
                Disparity(each_field).res_mult_gauss = Disparity(each_field).response .* Disparity(each_field).Gauss_in_RF;
            end
            %figure
            %surf(Disparity(2).res_mult_gauss);
            
            %sum together the responses for each individual surround
            Surround_response  = 0;
            for each_field = 2:num_surround+1
                Surround_response = Surround_response + sum(Disparity(each_field).res_mult_gauss(:));
            end
            
%            if (unique_mean_disp(mdisp) == -.35) %| (unique_mean_disp(mdisp) == -.05)
%                if unique_mag_disp(magdisp) == .2
%                    if (unique_disp_ang(dispang) == 90) | (unique_disp_ang(dispang) == 270)
%                        response_total_RF = zeros(length_grid);
%                        res_mult_gauss_Total = zeros(length_grid);
%                        for each_field = 2:num_surround+1
%                            response_total_RF = response_total_RF+Disparity(each_field).response;
%                            res_mult_gauss_Total = res_mult_gauss_Total + Disparity(each_field).res_mult_gauss;
%                        end
                                                   
%                        figure
%                        surf(response_total_RF)
%                        figure
%                        surf(res_mult_gauss_Total)
%                        figure
%                        surf(disparity_grid)
%                    end
%                end
%            end

            Center_response = sum(Disparity(1).res_mult_gauss(:));
            Total_Response = Center_response + Surround_response;
            
            resp_angle(magdisp, dispang) = Total_Response;
        end %end for each angle
    end %end for each magnitude of disparity
    tot_response{mdisp} = resp_angle;
end %end mean disparity

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

num_reps = 1;
final_fig = figure;

for i=1:length(unique_mag_disp)
    subplot(2, 1, i);
    for j=1:length(tot_response)
        temp = tot_response{j};
        hold on
            
        plot_x=[];
        for ij=1:num_reps
            plot_x = [plot_x; unique_disp_ang];
        end
        plot_y = temp(i,:)';
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
    close(gcf)
    height = axis;
    string = sprintf('P Values = %1.4f %1.4f %1.4f', p(i,:));
    text(height(1)+2, 1, string, 'FontSize', 8);
         
end %end i loop mag_disp loop

%print out a legend graph that lists each mean disparity and the line type with which they are drawn
leg_x = [0 1];
leg_y = [0 0];
legend_figure = figure;
set(legend_figure, 'Position', [950 773 100 100], 'Name', 'Line Legend', 'MenuBar', 'none');
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