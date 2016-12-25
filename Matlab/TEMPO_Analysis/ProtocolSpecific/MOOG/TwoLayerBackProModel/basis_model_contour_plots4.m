
function basis_model_contour_plots2(unique_azimuth, unique_elevation, ...
    unique_point_azimuth, unique_point_elevation, r1_resp, r2_resp, r3_resp, r4_resp, r5_resp,...
    r6_resp, r7_resp, r8_resp, r9_resp, r10_resp, fig_num, tizzy)

plot_clabels = 0;   % whether to plot data labels in contour plots
plot_clabels2 = 0;   % whether to plot data labels in contour plots
plot_str_len = 70;  % maximum string length for special statistic text

% convert the azimuth and elevation into a format amenable to plotting
unique_azimuth = sort(unique_azimuth)' - 90;  % is the sort redundant?
unique_elevation = sort(unique_elevation)';
unique_azimuth = [unique_azimuth unique_azimuth(1) + 360];  % wrap around sphere
[plot_azimuth plot_elevation] = meshgrid(unique_azimuth, unique_elevation);
num_grid_azimuth = length(unique_azimuth);
num_grid_elevation = length(unique_elevation);
plot_mean_data = zeros(num_grid_elevation,num_grid_azimuth);
plot_fit_data = zeros(num_grid_elevation,num_grid_azimuth);
plot_residual_data = zeros(num_grid_elevation,num_grid_azimuth);

% set the bounds and ticks for the axis on the contour plots, 
% depending on how we formatted the data for plotting above.
x_min = -90; x_max = 270; y_min = -90; y_max = 90;
x_tick = [-90:45:270]; y_tick = [-90:45:90];

% get the min, max, and range of the data and the fitted data
% over repititions, gaze angles, and points
max_data = max( [max(r1_resp) max(r2_resp) max(r3_resp) max(r4_resp) max(r5_resp)] );
min_data = min( [min(r1_resp) min(r2_resp) min(r3_resp) min(r4_resp) min(r5_resp)] );
range_data = max_data - min_data;

% specify the elevations to draw contour lines at for data and residual
% plots.  specifying the number of contour lines to draw does not take
% the min and max range of all the data into account, so it is
% difficult to compare between contours.
data_clines = [min_data:range_data/9:max_data];  % 0.99 to see max and min points
residual_clines = [-range_data/2:range_data/9:range_data/2];        

% create a figure 
allfigs = get(0,'children');
fig_exists = any(allfigs == 10+fig_num);
hmain = figure(10 + fig_num); h = hmain;
clf reset;  % clear the figure
if ~fig_exists
    ss = get(0,'ScreenSize');
    pos = [0.2*ss(3) 0.1*ss(4) 0.8*ss(3) 0.8*ss(4)];
    set(h,'Position',pos);
end

% convert the mean and fit data into a format amenable to plotting
for m=1:num_grid_elevation
    for n=1:num_grid_azimuth
        % i can not think of a cleaner way to do this stupid
        % all azimuths at the pole thing.
        e = 1e-10;  % floating point tolerance
        select = ((abs(unique_point_azimuth - mod(unique_azimuth(n),360)) < e | ...
            unique_elevation(m) == 90 | unique_elevation(m) == -90) & ...
            abs(unique_point_elevation - unique_elevation(m)) < e);
        plot_r1_data(m,n) = r1_resp(select);
        plot_r2_data(m,n) = r2_resp(select);
        plot_r3_data(m,n) = r3_resp(select);
        plot_r4_data(m,n) = r4_resp(select);
        plot_r5_data(m,n) = r5_resp(select);
        plot_r6_data(m,n) = r6_resp(select);
        plot_r7_data(m,n) = r7_resp(select);
        plot_r8_data(m,n) = r8_resp(select);
        plot_r9_data(m,n) = r9_resp(select);
        plot_r10_data(m,n) = r10_resp(select);
%         plot_residual_data(m,n) = (plot_mean_data(m,n) - plot_fit_data(m,n));
    end
end

% if isempty(data_clines)
%     return;
% end

num_rows = 2;
num_cols = 5;

% plot the r1 data
subplot(num_rows,num_cols,1);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r1_data,data_clines);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;    
%xlabel('azimuth'); ylabel('elevation');
subplot(num_rows,num_cols,6);
%[C h] = contourf('v6',plot_azimuth,plot_elevation,plot_r6_data);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r6_data);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); if plot_clabels2, clabel(C,h), end;    
% %xlabel('azimuth'); ylabel('elevation');

subplot(num_rows,num_cols,2);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r2_data,data_clines);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;
%xlabel('azimuth'); ylabel('elevation'); 
subplot(num_rows,num_cols,7);
%[C h] = contourf('v6',plot_azimuth,plot_elevation,plot_r7_data);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r7_data);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); if plot_clabels2, clabel(C,h), end;    

subplot(num_rows,num_cols,3);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r3_data,data_clines);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;
title(tizzy);
%xlabel('azimuth'); ylabel('elevation'); 
subplot(num_rows,num_cols,8);
%[C h] = contourf('v6',plot_azimuth,plot_elevation,plot_r8_data);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r8_data);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); if plot_clabels2, clabel(C,h), end;    

subplot(num_rows,num_cols,4);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r4_data,data_clines);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;
%xlabel('azimuth'); ylabel('elevation'); 
subplot(num_rows,num_cols,9);
%[C h] = contourf('v6',plot_azimuth,plot_elevation,plot_r9_data);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r9_data);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); if plot_clabels2, clabel(C,h), end;    

subplot(num_rows,num_cols,5);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r5_data,data_clines);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;
%xlabel('azimuth'); ylabel('elevation'); 
subplot(num_rows,num_cols,10);
%[C h] = contourf('v6',plot_azimuth,plot_elevation,plot_r10_data);
[C h] = contourf(plot_azimuth,plot_elevation,plot_r10_data);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); if plot_clabels2, clabel(C,h), end;    

% % plot the first residual
% subplot(1,3,3);
% [C h] = contourf(plot_azimuth,plot_elevation,plot_residual_data,residual_clines);
% set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
% axis([x_min,x_max,y_min,y_max]); caxis([-range_data/2 range_data/2]); if plot_clabels, clabel(C,h), end;
% xlabel('azimuth'); ylabel('elevation');
