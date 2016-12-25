
num_unique_elevation = 17;  % number of latitudes to generate (including poles)
num_unique_azimuth = 32;    % number of longitutdes to generate
unique_azimuth = (0:1:num_unique_azimuth-1)'/num_unique_azimuth*360;
unique_elevation = (-num_unique_elevation+1:2:num_unique_elevation-1)'/(num_unique_elevation-1)*90;
[el az] = meshgrid(unique_elevation(2:end-1), unique_azimuth);  % exclude the poles
unique_point_azimuth = [0 az(:)' 0];  % put the poles back
unique_point_elevation = [-90 el(:)' 90];
num_unique_points = length(unique_point_elevation);  % or azimuth
unique_point_azimuth_r = unique_point_azimuth/180*pi;
unique_point_elevation_r = unique_point_elevation/180*pi;

az = unique_point_azimuth_r;
el = unique_point_elevation_r;
mean_az = pi;
mean_el = -pi/4+.04;
std_az = pi/4;
std_el = pi/6;

t0 = 0;
p0 = pi/3;

sel = (az >= 3*pi/2);
az(sel) = az(sel) - 2*pi;

[az el] = Curvefit_rotate_coords(az,el,t0,p0);

%g = exp( -(el-mean_el).^2./(2*std_el.^2) );
%g = exp( -(az-mean_az).^2./(2*std_az.^2) );
%g = exp( -(az-mean_az).^2./(2*std_az.^2) - (el-mean_el).^2./(2*std_el.^2) );
g = exp( -az.^2./(2*std_az.^2) - el.^2./(2*std_el.^2) );
%g = exp( -(cos(az).*cos(el)-mean_az).^2./(2*std_az.^2) - (cos(az).*sin(el)-mean_el).^2./(2*std_el.^2) );
h = cos(az).*cos(el); h(h < 0) = 0;

%g = az;
%h = el;

% sel = (g > 0.01);
% g(sel) = 0.01;

basis_model_contour_plots(unique_azimuth, unique_elevation, ...
    unique_point_azimuth, unique_point_elevation, g, h)

