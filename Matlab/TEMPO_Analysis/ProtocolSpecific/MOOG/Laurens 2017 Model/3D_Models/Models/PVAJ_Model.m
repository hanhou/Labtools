% Position-Velocity-acceleration-jerk model for 3D tuning 20170603LBY
% a is aeter set
% u_azi is unique azimuth ( 0,45,90,135,180,225,270,315 )
% u_ele is unique elevation ( 0, -+45, -+90 )
% t is PSTH time points

function r = PVAJ_Model(a,st_data)

u_ele = st_data(1:5);
u_azi = st_data(6:13);
t = st_data(14:end);

stim_sig = sqrt(sqrt(2))/6;

%position model
% time profile
pos_time = vel_func(a(3:4), t);
% spatial profiles
ele_azi_p = cos_tuning(a(4:6), [u_ele; u_azi]);
ele_azi_p = reshape(ele_azi_p, length(u_azi), length(u_ele));

% velocity model
% time profile
vel_time = vel_func(a(3:4), t);
% spatial profiles
ele_azi_v = cos_tuning(a(8:10), [u_ele; u_azi]);
ele_azi_v = reshape(ele_azi_v, length(u_azi), length(u_ele));

% acceleration model
%time profile
acc_time = acc_func(a(3:4), t);
%spatial profiles
ele_azi_a = cos_tuning(a(12:14), [u_ele; u_azi]);
ele_azi_a = reshape(ele_azi_a, length(u_azi), length(u_ele));

% jerk model
%time profile
jerk_time = jerk_func(a(3:4), t);
%spatial profiles
ele_azi_j = cos_tuning(a(16:18), [u_ele; u_azi]);
ele_azi_j = reshape(ele_azi_j, length(u_azi), length(u_ele));


%compute results
r = zeros(size(ele_azi_v,1), size(ele_azi_v,2), length(vel_time));
for i=1:size(r,1),
    for j=1:size(r,2),
        rr =a(1)*( ...
            (1-a(22))*( ...
            (1-a(21))*( ...
            a(20)*ele_azi_p(i,j)*pos_time + ...
            (1-a(20))*ele_azi_v(i,j)*vel_time) + ...
            a(21)*ele_azi_a(i,j)*acc_time) + ...
            a(22)*ele_azi_j(i,j)*jerk_time) + ...
            a(2);
        r(i,j,:) = rr;
    end
end

end
