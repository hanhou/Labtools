% This is to calculate shift ratio for hidden units under vestibular,
% visual condition, horizontally and vertically  -GY 12/05/05
% only consider the largest angle, temporally, 40 degree horizontal and
% vertical
for k = 1 : 3
	for n = 1 : 150 % num of hidden units
    gg(1:26) = hidden_outputs(k,3,:,n);
    gg =gg+1;
	[x, y, z] = sph2cart(unique_point_azimuth_r, unique_point_elevation_r, gg );
	x = sum(x); y = sum(y); z = sum(z);
	[azz eee magg] = cart2sph(x,y,z);
	HTI_condition(k,n) = magg / sum(abs( gg ));  % norm of vector sums divided by sum of vector norms    
    maxx(k,n) = max(hidden_outputs(k,3,:,n));
    minn(k,n) = min(hidden_outputs(k,3,:,n));
	end
end
% dlmwrite('HTI.txt',HTI_condition');
% dlmwrite('maxx.txt',maxx');
% dlmwrite('minn.txt',minn');

% visual
pref_n45_dir_visual_az = pref_n45_dir_visual_az*180/pi;
pref_dir_visual_az = pref_dir_visual_az*180/pi;
pref_n45_dir_visual_az(find(pref_n45_dir_visual_az < 0)) = 360 + pref_n45_dir_visual_az(find(pref_n45_dir_visual_az < 0));
shift_az = abs(pref_n45_dir_visual_az - pref_dir_visual_az);
shift_az(find(shift_az>180)) = 360-shift_az(find(shift_az>180));
shift_vi = shift_az;

% vestibular
pref_n45_dir_vestib_az = pref_n45_dir_vestib_az*180/pi;
pref_dir_vestib_az = pref_dir_vestib_az*180/pi;
pref_n45_dir_vestib_az(find(pref_n45_dir_vestib_az < 0)) = 360 + pref_n45_dir_vestib_az(find(pref_n45_dir_vestib_az < 0));
shift_az = abs(pref_n45_dir_vestib_az - pref_dir_vestib_az);
shift_az(find(shift_az>180)) = 360-shift_az(find(shift_az>180));
shift_ve = shift_az;

shift_ratio(:,1) = shift_ve';
shift_ratio(:,2) = shift_vi';
shift_ratio

