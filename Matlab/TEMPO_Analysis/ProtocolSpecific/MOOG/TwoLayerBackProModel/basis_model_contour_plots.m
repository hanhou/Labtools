
function basis_model_contour_plots(unique_azimuth, unique_elevation, ...
    unique_point_azimuth, unique_point_elevation, expected_hc_resp, actual_hc_resp)

plot_clabels = 0;   % whether to plot data labels in contour plots
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
max_data = max( [max(expected_hc_resp) max(actual_hc_resp)] );
min_data = min( [min(expected_hc_resp) min(actual_hc_resp)] );
range_data = max_data - min_data;

% specify the elevations to draw contour lines at for data and residual
% plots.  specifying the number of contour lines to draw does not take
% the min and max range of all the data into account, so it is
% difficult to compare between contours.
data_clines = [min_data:range_data/9:max_data];  % 0.99 to see max and min points
residual_clines = [-range_data/2:range_data/9:range_data/2];        

% create a figure 
hmain = figure(10); h = hmain;
clf reset;  % clear the figure
ss = get(0,'ScreenSize');
pos = [0.2*ss(3) 0.1*ss(4) 0.8*ss(3) 0.8*ss(4)];
set(h,'Position',pos);

% convert the mean and fit data into a format amenable to plotting
for m=1:num_grid_elevation
    for n=1:num_grid_azimuth
        % i can not think of a cleaner way to do this stupid
        % all azimuths at the pole thing.
        e = 1e-10;  % floating point tolerance
        select = ((abs(unique_point_azimuth - mod(unique_azimuth(n),360)) < e | ...
            unique_elevation(m) == 90 | unique_elevation(m) == -90) & ...
            abs(unique_point_elevation - unique_elevation(m)) < e);
        if sum(select) == 0;
            eat_shit_and_die = 1;
        end
        plot_mean_data(m,n) = expected_hc_resp(select);
        plot_fit_data(m,n) = actual_hc_resp(select);
        plot_residual_data(m,n) = (plot_mean_data(m,n) - plot_fit_data(m,n));
    end
end

% plot the mean data
% subplot(1,3,1);
% [C h] = contourf(plot_azimuth,plot_elevation,plot_mean_data);
% set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
% axis([x_min,x_max,y_min,y_max]); %caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;    
% xlabel('azimuth'); ylabel('elevation');
subplot(1,3,1);
[C h] = contourf(plot_azimuth,plot_elevation,plot_mean_data,data_clines);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;    
xlabel('azimuth'); ylabel('elevation');
%                 str = {sprintf('%s', FILE) sprintf('gaze_angle=%d vary_fixation_%s', ...
%                         unique_gaze_angle(gas), fixation_type_str)}; 
%                 h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); 

% plot the first fit
subplot(1,3,2);
[C h] = contourf(plot_azimuth,plot_elevation,plot_fit_data,data_clines);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;
xlabel('azimuth'); ylabel('elevation'); 
%                 str = {sprintf('%s', func2str(fitting_models(j)))};
%                 if gas == 1, h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); end

% plot the first residual
% subplot(1,3,3);
% [C h] = contourf(plot_azimuth,plot_elevation,plot_residual_data);
% set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
% axis([x_min,x_max,y_min,y_max]); %caxis([-range_data/2 range_data/2]); if plot_clabels, clabel(C,h), end;
% xlabel('azimuth'); ylabel('elevation');
subplot(1,3,3);
[C h] = contourf(plot_azimuth,plot_elevation,plot_residual_data,residual_clines);
set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
axis([x_min,x_max,y_min,y_max]); caxis([-range_data/2 range_data/2]); if plot_clabels, clabel(C,h), end;
xlabel('azimuth'); ylabel('elevation');
%                 str = [sprintf('%s SSE=%.5f ', func2str(fitting_models(j)), fit_sse(i,j,gas)) ];
%                 str = [str repmat([' '], 1, plot_str_len - size(str,2))];
%                 pstr{gas} = [pstr{gas};str];
