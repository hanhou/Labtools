function surround_test(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);

%first thing to do is to create a series of matrices for each surround area.
%1) gaussian receptive field 
%2) response of the surround using a the inverted version of the tuning curve of the CRF
%3) multiply the response by the gaussian
%4) now create a matrix large enough to hold the CRF in the center and have space for the surround fields
%5) add the response matrices of each of the surrounds and the CRF into the large matrix
%6) to calculate the entire spike rate of the neuron for this stimulus, make all negative values zero
%   and then add together all the cells in the matrix.


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
Total_RF_size = max_size_vector(1);

CRF_Size = max_size_vector(2);
Surround_Size = 8;
%Total_RF_Size = CRF_Size + (2*Surround_Size);
res = .5;

unique_mean_disp = [-0.4 0 0.4]';
px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';         %disparity values
pystart = [100 100 100 120 130 150 130 120 100 100 100]';   %neuronal responses to the disparities listed above

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

Ksurr = -.11;
a = sqrt(-(Surround_gauss_x^2)/log(.1))*2;
b = sqrt(-(Surround_gauss_y^2)/log(.1))*2;

num_surround = 6;
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
            Disparity(i+1).Gauss_grid(j, k) = Ksurr*(exp(-((Surround_x_ind^2)/a + (Surround_y_ind^2)/b)));
        end
    end
end


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

Disparity(1).px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
Disparity(1).py = [100 100 100 120 130 150 130 120 100 100 100]';
Disparity(1).present = zeros(length_grid);
Disparity(1).present(CRF_x_start:CRF_x_start-1+CRF_size, CRF_y_start:CRF_y_start-1+CRF_size) = 1;
Disparity(1).Gauss_in_RF = zeros(length_grid);
Disparity(1).Gauss_in_RF(CRF_x_start:CRF_x_start-1+CRF_size, CRF_y_start:CRF_y_start-1+CRF_size) = Disparity(1).Gauss_grid;

shift_val = round(.2*Sur_size);

%top left surround field
Sur_x_start = total_RF_x_ctr-round(Sur_size/2)-shift_val;
Sur_y_start = total_RF_y_ctr-round(CRF_size/2)-round(Sur_size/2)+1+shift_val;
total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size)+Disparity(2).Gauss_grid;

Disparity(2).px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
Disparity(2).py = [150 130 120 105 100 90 100 105 120 130 150]';
Disparity(2).present = zeros(length_grid);
Disparity(2).present(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = 1;
Disparity(2).Gauss_in_RF = zeros(length_grid);
Disparity(2).Gauss_in_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = Disparity(2).Gauss_grid;

%top right surround field
Sur_x_start = total_RF_x_ctr-round(Sur_size/2)+shift_val;
Sur_y_start = total_RF_y_ctr-round(CRF_size/2)-round(Sur_size/2)+1+shift_val;
total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size)+Disparity(3).Gauss_grid;

Disparity(3).px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
Disparity(3).py = [150 130 120 105 100 90 100 105 120 130 150]';
Disparity(3).present = zeros(length_grid);
Disparity(3).present(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = 1;
Disparity(3).Gauss_in_RF = zeros(length_grid);
Disparity(3).Gauss_in_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = Disparity(3).Gauss_grid;

%bottom left surround field
Sur_x_start = total_RF_x_ctr-round(Sur_size/2)-shift_val;
Sur_y_start = total_RF_y_ctr+round(CRF_size/2)-round(Sur_size/2)-1 -shift_val;
total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size)+Disparity(4).Gauss_grid;

Disparity(4).px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
Disparity(4).py = [150 130 120 105 100 90 100 105 120 130 150]';
Disparity(4).present = zeros(length_grid);
Disparity(4).present(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = 1;
Disparity(4).Gauss_in_RF = zeros(length_grid);
Disparity(4).Gauss_in_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = Disparity(4).Gauss_grid;

%bottom right surround field
Sur_x_start = total_RF_x_ctr-round(Sur_size/2)+shift_val;
Sur_y_start = total_RF_y_ctr+round(CRF_size/2)-round(Sur_size/2)-1-shift_val;
total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size)+Disparity(5).Gauss_grid;

Disparity(5).px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
Disparity(5).py = [150 130 120 105 100 90 100 105 120 130 150]';
Disparity(5).present = zeros(length_grid);
Disparity(5).present(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = 1;
Disparity(5).Gauss_in_RF = zeros(length_grid);
Disparity(5).Gauss_in_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = Disparity(5).Gauss_grid;

%right side surround field
Sur_x_start = total_RF_x_ctr+round(CRF_size/2)-round(Sur_size/2)-2;
Sur_y_start = total_RF_y_ctr-round(Sur_size/2);
total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size)+Disparity(6).Gauss_grid;

Disparity(6).px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
Disparity(6).py = [150 130 120 105 100 90 100 105 120 130 150]';
Disparity(6).present = zeros(length_grid);
Disparity(6).present(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = 1;
Disparity(6).Gauss_in_RF = zeros(length_grid);
Disparity(6).Gauss_in_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = Disparity(6).Gauss_grid;

%left side surround field
Sur_x_start = total_RF_x_ctr-round(CRF_size/2)-round(Sur_size/2)+2;
Sur_y_start = total_RF_y_ctr-round(Sur_size/2);
total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = total_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size)+Disparity(7).Gauss_grid;

Disparity(7).px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
Disparity(7).py = [150 130 120 105 100 90 100 105 120 130 150]';
Disparity(7).present = zeros(length_grid);
Disparity(7).present(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = 1;
Disparity(7).Gauss_in_RF = zeros(length_grid);
Disparity(7).Gauss_in_RF(Sur_x_start:Sur_x_start-1+Sur_size, Sur_y_start:Sur_y_start-1+Sur_size) = Disparity(7).Gauss_grid;


distance_grid = total_RF;

x = zeros(length(distance_grid), 1);
y = zeros(length(distance_grid), 1);

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
        if distance_grid(i,j) > 8
            distance_grid(i,j) = NaN;
        end

    end
end

unique_mean_disp = 0.2;
unique_disp_ang = 90;
unique_mag_disp = 0.2;

px = [-1.8 -1.6 -1.2 -.8 -.4 0 .4 .8 1.2 1.6 1.8]';
pystart = [100 100 100 120 130 150 130 120 100 100 100]';

py = pystart;

%create a matrix where each cell contains the disparity of that location
disparity_grid = zeros(length(x)) + unique_mean_disp;
response_grid = zeros(length(x));

for each_field=1:num_surround+1
    Disparity(each_field).response = zeros(length_grid);
end


%%add in the error for the spikes here
%randn('state',sum(100*clock))
%perr = randn(1, length(pystart));
%perr = perr' .* sigma;
%py = perr + mu;

angle = unique_disp_ang * (pi/180);

for l = 1:length(disparity_grid) %for each x pt on the grid
    %get the x coordinate of the grid
    x_ind = x(l);
    for m = 1:length(disparity_grid) %for each y point on the grid
        %get the y coordinate of the grid
        y_ind = y(m);
        if ~isnan(distance_grid(l,m))
            if(angle==0)
                disp_val = y_ind * unique_mag_disp;
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
                    disp_val = y_tilt * unique_mag_disp;
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
                disparity_grid(l,m) = 0;
            end
        end %end isnan if
    end %end m
end %end l

figure
surf(disparity_grid)

%not multiply the response of each field by their respective gaussians
for each_field = 1:num_surround+1
    Disparity(each_field).res_mult_gauss = zeros(length_grid);
    Disparity(each_field).res_mult_gauss = Disparity(each_field).response .* Disparity(each_field).Gauss_in_RF;
end

%sum together the responses for each individual surround
Surround_response  = 0;
for each_field = 2:num_surround+1
    Surround_response = Surround_response + sum(Disparity(each_field).res_mult_gauss(:));
end

Surround_response
Center_response = sum(Disparity(1).res_mult_gauss(:))

Total_Response = Center_response + Surround_response
figure
surf(total_RF)