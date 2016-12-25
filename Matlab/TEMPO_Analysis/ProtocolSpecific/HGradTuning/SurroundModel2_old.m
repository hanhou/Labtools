function SurroundModel2(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, HDisparity, SurrMap, UseSyncPulses);
TEMPO_Defs;

%following function checks offset times and outputs starting and ending analysis (in spike bins)
num_trials = size(data.event_data, 3);
[StartOffsetBin StopOffsetBin] = CheckTimeOffset(data, num_trials, StartCode, StopCode, StartOffset, StopOffset, UseSyncPulses);

if (~isempty(data.spike_data))
   %compute the firing rate over all trials during the period between StartCode and StopCode
   data.spike_rates = ComputeSpikeRates(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin);
end
   
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
ap_size = data.dots_params(DOTS_AP_XSIZ,:,PATCH1);
unique_size = munique(ap_size(~null_trials)');
   
%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (ap_size == data.one_time_params(NULL_VALUE)) );
   
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(ap_size);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%this will be the diameter of the grid
CRF_Size = SurrMap.ctrsize;

Total_RF_Size = SurrMap.ctrsize + (SurrMap.surrsize * 2);
Surround_Size = SurrMap.surrsize;
%number of surround patches
num_surr = length(SurrMap.angles);

%retrieve the tuning curves for the CRF and the surrounds
CRF.px = HDisparity.px;
CRF.py = HDisparity.py;
CRF.plot_x = HDisparity.plot_x;
CRF.plot_y = HDisparity.plot_y;
for i = 1:num_surr
   SurrPatch.px(:, i) = SurrMap.px(:, i);
   SurrPatch.py(:, i) = SurrMap.py(:, i);
   SurrPatch.plot_x(:, i) = SurrMap.plot_x(:, i);
   SurrPatch.plot_y(:, i) = SurrMap.plot_y(:, i);
end

%figure that displays tuning curves in appropriate locations
%surr_tuning_fig=figure;

%crf_axis = axes('position', [.35 .32 .3 .3]);
%surr_axis(1) = axes('position', [.68 .32 .3 .3]);
%surr_axis(2) = axes('position', [.52 .64 .3 .3]);
%surr_axis(3) = axes('position', [.18 .64 .3 .3]);
%surr_axis(4) = axes('position', [.02 .32 .3 .3]);
%surr_axis(5) = axes('position', [.18 0 .3 .3]);
%surr_axis(6) = axes('position', [.52 0 .3 .3]);


%get ctr only response from surround tuning data
CTR_only = SurrMap.ctronly;

dogabor = 0;
dosurroundgabor = 0;

if dogabor == 1
   means = [CRF.px CRF.py];
   raw = [CRF.plot_x' CRF.plot_y'];
   fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
   fixed_param_values = zeros(6,1); 
   [pars,freq] = gaborfit(means,raw, fixed_param_flags, fixed_param_values);
   CRF.x_interp = (CRF.px(1)): .01 : (CRF.px(length(CRF.px)));
   % Note: this func MUST match that in gaborfit.m
   CRF.y_interp =  pars(1) + pars(2)*exp(-0.5*((CRF.x_interp - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(CRF.x_interp - pars(3))+pars(6) );
   CRF.y_interp(CRF.y_interp < 0) = 0;
else
   CRF.x_interp = (CRF.px(1)): .01 : (CRF.px(length(CRF.px)));
   CRF.y_interp = interp1(CRF.px, CRF.py, CRF.x_interp, 'spline');  %spline fit
   CRF_max = max(CRF.y_interp);
   CRF_diff = CRF_max - CTR_only;
   CRF.y_interp = CRF.y_interp - CRF_diff;
end

tuning = figure;
hold on
plot(CRF.x_interp, CRF.y_interp, 'k--');
%plot crf tuning curve in center
%figure(surr_tuning_fig)
%axes(crf_axis)
%plot(CRF.x_interp, CRF.y_interp, 'k--');

for i=1:num_surr
   if dosurroundgabor == 1
      means = [SurrPatch.px(:, i) SurrPatch.py(:, i)];
      raw = [SurrPatch.plot_x(i, :)' SurrPatch.plot_y(i, :)'];
      fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
      fixed_param_values = zeros(6,1); 
      [pars,freq] = gaborfit(means,raw, fixed_param_flags, fixed_param_values);
      SurrPatch.x_interp(i,:) = (SurrPatch.px(1,i)): .01 : (SurrPatch.px(length(SurrPatch.px(:,i)), i));
      SurrPatch.y_interp(i,:) =  pars(1) + pars(2)*exp(-0.5*((SurrPatch.x_interp(i,:) - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(SurrPatch.x_interp(i,:) - pars(3))+pars(6) );
      SurrPatch.y_interp(SurrPatch.y_interp(i,:) < 0) = 0;
   else
      SurrPatch.x_interp(i,:) = (SurrPatch.px(1,i)): .01 : (SurrPatch.px(length(SurrPatch.px(:,i)), i));
      SurrPatch.y_interp(i,:) = interp1(SurrPatch.px(:, i), SurrPatch.py(:, i), SurrPatch.x_interp(i,:), 'spline') - CTR_only;  %spline fit
   end
   figure(tuning)
   hold on
   plot(SurrPatch.x_interp(i,:),SurrPatch.y_interp(i,:), lines{i});
end



figure(tuning)
angle1 = sprintf('%d', SurrMap.angles(1));
angle2 = sprintf('%d', SurrMap.angles(2));
angle3 = sprintf('%d', SurrMap.angles(3));
angle4 = sprintf('%d', SurrMap.angles(4));
angle5 = sprintf('%d', SurrMap.angles(5));
angle6 = sprintf('%d', SurrMap.angles(6));
legend('CRF', angle1, angle2, angle3, angle4, angle5, angle6, 0);

surr_max = max(SurrPatch.px(:,1));
surr_min = min(SurrPatch.px(:,1));

%draw a little schematic of what the surround and center look like
CRF_xdraw = 0-(CRF_Size/2);
CRF_ydraw = 0-(CRF_Size/2);
%circles = figure;
%hold on
%rectangle('Position', [-unique_size(1)/2 -unique_size(1)/2 unique_size(1) unique_size(1)], 'Curvature', [1 1], 'FaceColor', 'b');
%rectangle('Position', [CRF_xdraw CRF_ydraw CRF_Size CRF_Size], 'Curvature', [1 1], 'FaceColor','r');
%axis equal

Distance_between_Surr_and_Ctr_ctrs = (Surround_Size + CRF_Size)/2;

surr_xdraw = [];
surr_ydraw = [];
for i=1:num_surr
   ang = SurrMap.angles(i);
   ang_rad = ang * (pi/180);
   surr_xdraw(i) = Distance_between_Surr_and_Ctr_ctrs * cos(ang_rad) + 0 - Surround_Size/2;
   surr_ydraw(i) = Distance_between_Surr_and_Ctr_ctrs * sin(ang_rad) + 0 - Surround_Size/2;
 %  rectangle('Position', [surr_xdraw(i) surr_ydraw(i) Surround_Size Surround_Size], 'Curvature', [1 1]);
end
%resolution of grids in degrees
res = .1;

x_offset = abs(min(surr_xdraw)) + 1;
y_offset = abs(min(surr_ydraw)) + 1;

x_ind = fix((surr_xdraw + x_offset) * 1/res);
y_ind = fix((surr_ydraw + y_offset) * 1/res);

x_max = max(x_ind) + Surround_Size*(1/res);
y_max = max(y_ind) + Surround_Size*(1/res);

CRF_xind = fix((CRF_xdraw + x_offset) * 1/res);
CRF_yind = fix((CRF_ydraw + y_offset) * 1/res);

zero_xind = fix((0+x_offset) * 1/res);
zero_yind = fix((0+y_offset) * 1/res);


%big_circles = figure;
%hold on
%rectangle('Position', [-unique_size(1)/2 -unique_size(1)/2 unique_size(1) unique_size(1)], 'Curvature', [1 1], 'FaceColor', 'b');
%rectangle('Position', [CRF_xind CRF_yind CRF_Size*(1/res) CRF_Size*(1/res)], 'Curvature', [1 1], 'FaceColor','r');
%axis equal
%for i=1:num_surr
%   rectangle('Position', [x_ind(i) y_ind(i) Surround_Size*(1/res) Surround_Size*(1/res)], 'Curvature', [1 1]);
%end
%axis([0 x_max 0 y_max]);

%grid_size = max([x_max y_max]);
x_list = 1:x_max;
y_list = 1:y_max;
x_coord = (x_list-zero_xind) * res;
y_coord = (y_list-zero_yind) * res;

distance_grid = zeros(x_max, y_max);

for i = 1:length(x_coord)
   for j = 1:length(y_coord)
      distance_grid(i, j) = sqrt((x_coord(i)-0)^2 + (y_coord(j)-0)^2);
   end
end

grad_fig = figure;

disparity_and_angles = zeros(8, length(unique_disp_ang));

for magdisp =1:length(unique_mag_disp) %for each magnitude of disparity
   resp_mag_ang = zeros(length(unique_mag_disp), length(unique_disp_ang));
   for mdisp=1:length(unique_mean_disp)  %for each mean disparity...
      number_disp_outside_range_per_meandisp = 0;
      angle_and_response_data = [];
      for dispang =1:length(unique_disp_ang)
         %create a matrix where each cell contains the disparity of that location
         disparity_grid = zeros(x_max, y_max) + unique_mean_disp(mdisp);

         angle = unique_disp_ang(dispang) * (pi/180);

         for l = 1:x_max %for each x pt on the grid
            %get the x coordinate of the grid
            x_val = x_coord(l);
            for m = 1:y_max %for each y point on the grid
               %get the y coordinate of the grid
               y_val = y_coord(m);
               if(angle==0)
                  disp_val = y_val * unique_mag_disp(magdisp);
                  %is alpha != 0 then there will be a tilt to the gradient
               else
                  len = distance_grid(l,m);
                  if len ~= 0
                     theta = asin(y_val/len);
                     if ~isreal(theta)
                        theta = real(theta);
                     end
                     %//sin only goes from -90 to 90 degrees.  need to fix for the other two quadrants.
                     %//need to do quadrant corrections
                     %//cases (-, +), (-, -)
                     if ((x_val <= 0) & (y_val >= 0))
                        %//Quadrant III
                        %//theta is now greater than 90 degrees
                        temptheta = (pi/2)-theta;
                        theta = (pi/2) + temptheta;
                     elseif ((x_val <= 0) & (y_val <= 0))
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
            end %end m
         end %end l
         %now that we have the disparity gradient we can calculate the mean disparity
         %in each of the surround patches
         mean_disp_in_surr = [];
         resp = [];
         response = 0;
         
         %switches fo choosing which interpolation scheme to use
         dolinear = 1;
         docubic = 0;
         dospline = 0;
         
         number_disp_outside_range_per_angle = 0;
         
         %calculate mean disparity for each surround patch
         %then interpolate using the surround tuning curves to calculate the expected response of the surround
         %to that disparity
         %num_surr = 6;
         for i=1:num_surr
            disp_in_surr = disparity_grid(x_ind(i):x_ind(i)+Surround_Size*(1/res), y_ind(i):y_ind(i)+Surround_Size*(1/res));
            mean_disp_in_surr(i) = mean(disp_in_surr(:));
            if mean_disp_in_surr(i)>surr_max | mean_disp_in_surr(i)<surr_min
               number_disp_outside_range_per_angle = number_disp_outside_range_per_angle + 1;
            end
            if dolinear == 1
               resp(i) = interp1(SurrPatch.x_interp(i,:), SurrPatch.y_interp(i,:), mean_disp_in_surr(i), 'linear', 'extrap');
            elseif docubic == 1
               resp(i) = interp1(SurrPatch.x_interp(i,:), SurrPatch.y_interp(i,:), mean_disp_in_surr(i), 'cubic', 'extrap');
            elseif dospline == 1
               resp(i) = interp1(SurrPatch.x_interp(i,:), SurrPatch.y_interp(i,:), mean_disp_in_surr(i), 'spline', 'extrap');
            end
         end
         
         %calculate the mean disparity in the CRF
         %then using the disparity tuning curve of the center
         %inpteroplate to find the expected response
         mean_disp_and_res = [mean_disp_in_surr' resp'];
         disp_in_CRF = disparity_grid(CRF_xind:CRF_xind+CRF_Size*(1/res) ,  CRF_yind:CRF_yind + CRF_Size*(1/res));
         mean_disp_in_CRF = mean(disp_in_CRF(:));
         if dolinear == 1
            CRF_response = interp1(CRF.x_interp, CRF.y_interp, mean_disp_in_CRF, 'linear', 'extrap');
         elseif docubic == 1
            CRF_response = interp1(CRF.x_interp, CRF.y_interp, mean_disp_in_CRF, 'cubic', 'extrap');
         elseif dospline == 1
            CRF_response = interp1(CRF.x_interp, CRF.y_interp, mean_disp_in_CRF, 'spline', 'extrap');
         end
         
         disparities_and_angle(1,dispang) = unique_disp_ang(dispang);
         disparities_and_angle(2,dispang) = mean_disp_in_CRF;
         disparities_and_angle(3:8,dispang) = mean_disp_in_surr';

         avg_response = CRF_response + sum(resp);%resp(1) + resp(2) + resp(3) + resp(4) + resp(5) + resp(6) + ;
         angle_and_response_data(dispang,:) = [unique_disp_ang(dispang) avg_response];

         
         number_disp_outside_range_per_meandisp = number_disp_outside_range_per_meandisp + number_disp_outside_range_per_angle;
      end %end for each angle
      
      figure(grad_fig)
      hold on
      PlotTuningCurve(angle_and_response_data(:, 1), angle_and_response_data(:, 2), symbols{mdisp}, lines{mdisp}, 1, 1);
      leg{mdisp} = sprintf('%1.3f, %d pts out of range', unique_mean_disp(mdisp), number_disp_outside_range_per_meandisp);
   end %end mean disparity
   axhandles = get(gca, 'Children');
   numlines = length(unique_mean_disp);
   leghandles = [];
   legindex = 1;
   for legendloop = numlines:-1:1
      leghandles(legindex) = axhandles((legendloop*3));
      legindex = legindex + 1;
   end
   
   legend(leghandles, leg, 0);
end %end for each magnitude of disparity
disparities_and_angle