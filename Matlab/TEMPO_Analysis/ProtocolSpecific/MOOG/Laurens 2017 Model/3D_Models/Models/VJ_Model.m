% Velocity-jerk model for 3D tuning 20170602LBY
% a is parameter set
% u_azi is unique azimuth ( 0,45,90,135,180,225,270,315 )
% u_ele is unique elevation ( 0, -+45, -+90 )
% t is PSTH time points

function r = VJ_Model(a,st_data)

u_ele = st_data(1:5);
u_azi = st_data(6:13);
t = st_data(14:end);

stim_sig = sqrt(sqrt(2))/6;

% velocity model
% time profile
vel_time = vel_func(a(3:4), t);
% spatial profiles
ele_azi_v = d_cos_tuning(a(5:12), [u_ele; u_azi]);
ele_azi_v = reshape(ele_azi_v, length(u_azi), length(u_ele));


% jerk model
%time profile
jerk_time = jerk_func(a(3:4), t);
%spatial profiles
ele_azi_j = d_cos_tuning(a(13:20), [u_ele; u_azi]);
ele_azi_j = reshape(ele_azi_j, length(u_azi), length(u_ele));


%compute results
r = zeros(size(ele_azi_v,1), size(ele_azi_v,2), length(vel_time));
for i=1:size(r,1),
    for j=1:size(r,2),
        rr =a(1)*(a(21)*ele_azi_v(i,j)*vel_time + (1-a(21))*ele_azi_j(i,j)*jerk_time)+ a(2);
        r(i,j,:) = rr;
    end
end

end
