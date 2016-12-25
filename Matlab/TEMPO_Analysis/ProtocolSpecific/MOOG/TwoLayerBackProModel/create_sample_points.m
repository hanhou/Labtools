
% create data points about the sphere where we will sample the error used
% for training the network.
global num_unique_elevation;  % this is a hack to set this param from a top level looping script
global num_unique_azimuth;  % this is a hack to set this param from a top level looping script
unique_azimuth = (0:1:num_unique_azimuth-1)'/num_unique_azimuth*360;
unique_elevation = (-num_unique_elevation+1:2:num_unique_elevation-1)'/(num_unique_elevation-1)*90;
[el az] = meshgrid(unique_elevation(2:end-1), unique_azimuth);  % exclude the poles
unique_point_azimuth = [0 az(:)' 0];  % put the poles back
unique_point_elevation = [-90 el(:)' 90];
num_unique_points = length(unique_point_elevation);  % or azimuth
unique_point_azimuth_r = unique_point_azimuth/180*pi;
unique_point_elevation_r = unique_point_elevation/180*pi;

num_gaze_samples = 5;  % make this independent of the number of gaze angle inputs
gas = [range_gaze_az(1):(range_gaze_az(2)-range_gaze_az(1))/(num_gaze_samples-1):range_gaze_az(2)];
%gas = [range_gaze_az(1):(range_gaze_az(2)-range_gaze_az(1))/(num_eye_neurons-1):range_gaze_az(2)];
num_gaze_angles = length(gas);
unique_gaze_angles_xr = [gas zeros(1,num_gaze_angles)];
unique_gaze_angles_yr = [zeros(1,num_gaze_angles) gas];
num_gaze_angles = length(unique_gaze_angles_xr);
