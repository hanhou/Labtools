function [exp_plus_sph, exp_minus_sph, shift_ratio_plus, shift_ratio_minus, shift_ratio_plusminus, noise_ratio_plus, noise_ratio_minus, noise_ratio_plusminus] = rotation(zero, actual_plus, actual_minus, alph, gaze_dir, fignum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROTATION (zero, actual_plus, actual_minus, alph, gaze_dir, fignum) ** but now without plots, for bootstrap and net_data **
% Given the calculated vectorsum in a 3D_Vary_Fixation experiment, 
% computes rotated vectors for the 'eye-centered' prediction and finds the 
% difference angle between expected and actual vectors, then plots all five 
% vectors.   CRF 12/29/03
% Added shift ratios.  CRF 9/04

% for us, a 'positive' elevation is actually downward, so:
zero(2) = -zero(2);
actual_plus(2) = -actual_plus(2);
actual_minus(2) = -actual_minus(2);

% change angle to radians and define rotation matrices
alph = alph * pi/180;
if gaze_dir == 0       % horizontal gaze
  rotate_plus = [cos(alph),sin(alph),0; -sin(alph),cos(alph),0; 0,0,1];
  rotate_minus = [cos(-alph),sin(-alph),0; -sin(-alph),cos(-alph),0; 0,0,1];
elseif gaze_dir == 1   % vertical gaze
  rotate_plus = [1,0,0; 0,cos(alph),-sin(alph); 0,sin(alph),cos(alph)];
  rotate_minus = [1,0,0; 0,cos(-alph),-sin(-alph); 0,sin(-alph),cos(-alph)];
else 
  disp('Error: Value of GAZE_DIR must be 0 or 1');
  return;
end

% change vectors to cartesian coords (to radians first), and rotate
[zero_x,zero_y,zero_z] = sph2cart(zero(1)*pi/180,zero(2)*pi/180,zero(3));
[actual_plus_x,actual_plus_y,actual_plus_z] = sph2cart(actual_plus(1)*pi/180,actual_plus(2)*pi/180,actual_plus(3));
[actual_minus_x,actual_minus_y,actual_minus_z] = sph2cart(actual_minus(1)*pi/180,actual_minus(2)*pi/180,actual_minus(3));

expected_plus = rotate_plus * [zero_x;zero_y;zero_z];
expected_plus = expected_plus';
expected_minus = rotate_minus * [zero_x;zero_y;zero_z];
expected_minus = expected_minus';

% % add a point at 0 and plot the five vectors
% expected_plus_0 = [expected_plus(1,:);0,0,0];
% expected_minus_0 = [expected_minus(1,:);0,0,0];
% zero_x_0 = [zero_x;0];
% zero_y_0 = [zero_y;0];
% zero_z_0 = [zero_z;0];
% actual_plus_x_0 = [actual_plus_x;0];
% actual_plus_y_0 = [actual_plus_y;0];
% actual_plus_z_0 = [actual_plus_z;0];
% actual_minus_x_0 = [actual_minus_x;0];
% actual_minus_y_0 = [actual_minus_y;0];
% actual_minus_z_0 = [actual_minus_z;0];
% 
% figure(fignum);
% axes('position',[0.05,0.05, 0.8,0.8] );
% set(fignum,'Position', [800*(fignum-4),0, 800,650], 'Name', 'Actual and expected rotations of vectorsum');
% plot3(expected_plus_0(:,1), expected_plus_0(:,2), expected_plus_0(:,3), ...             %blue
%       zero_x_0, zero_y_0, zero_z_0, ...                                                 %green
%       expected_minus_0(:,1), expected_minus_0(:,2), expected_minus_0(:,3), ...          %red
%       actual_plus_x_0, actual_plus_y_0, actual_plus_z_0, ...                            %cyan
%       actual_minus_x_0, actual_minus_y_0, actual_minus_z_0, ...                         %magenta
%       [-zero(3);zero(3)],[0;0],[0;0], 'k', [0;0], [-zero(3);zero(3)], [0;0], 'k', ...   
%       [0;0], [0;0], [-zero(3);zero(3)], 'k');                                           %black (pseudo-axes)
%   
% xlim([-zero(3),zero(3)]);
% ylim([-zero(3),zero(3)]);
% zlim([-zero(3),zero(3)]);

% convert back to spherical coords
[exp_plus_sph(1),exp_plus_sph(2),exp_plus_sph(3)] = cart2sph(expected_plus(1),expected_plus(2),expected_plus(3));
[exp_minus_sph(1), exp_minus_sph(2), exp_minus_sph(3)] = cart2sph(expected_minus(1),expected_minus(2),expected_minus(3));

%____________________________________________________
% ******* Begin shift ratio calculation *******

% first take the normal (cross product) of each pair of vectors, and use them (and their lengths) 
% in the formula from the website to get the plane projection of actuals
normal_plusminus = cross(expected_plus, expected_minus);  %(not needed yet)
normal_plus = cross(expected_plus, [zero_x,zero_y,zero_z]);
normal_minus = cross([zero_x,zero_y,zero_z], expected_minus);
 
normal_pm_length = sqrt(normal_plusminus(1)^2+normal_plusminus(2)^2+normal_plusminus(3)^2);
normal_p_length = sqrt(normal_plus(1)^2+normal_plus(2)^2+normal_plus(3)^2);
normal_m_length = sqrt(normal_minus(1)^2+normal_minus(2)^2+normal_minus(3)^2);

proj_plus = cross(normal_plus,(cross([actual_plus_x,actual_plus_y,actual_plus_z],normal_plus)/normal_p_length))/normal_p_length;
[proj_plus_sph(1),proj_plus_sph(2),proj_plus_sph(3)] = cart2sph(proj_plus(1),proj_plus(2),proj_plus(3));
proj_minus = cross(normal_minus,(cross([actual_minus_x,actual_minus_y,actual_minus_z],normal_minus)/normal_m_length))/normal_m_length;
[proj_minus_sph(1),proj_minus_sph(2),proj_minus_sph(3)] = cart2sph(proj_minus(1),proj_minus(2),proj_minus(3));

% calculate difference angles using the formula:
% cos(x) = sin(Ea)sin(Eb) + cos(Eb)sin(Ab)cos(Ea)sin(Aa) + cos(Eb)cos(Ab)cos(Ea)cos(Aa)
% where Ea and Eb is the elevation of vectors a and b, and Aa and Ab = azimuth of a and b

% (changing back and forth between deg and rad is a huge pain; I'll use temp vars)
zero_temp(1) = zero(1) * pi/180;
zero_temp(2) = zero(2) * pi/180;
zero_temp(3) = zero(3);
actual_plus_temp(1) = actual_plus(1) * pi/180;
actual_plus_temp(2) = actual_plus(2) * pi/180;
actual_plus_temp(3) = actual_plus(3);
actual_minus_temp(1) = actual_minus(1) * pi/180;
actual_minus_temp(2) = actual_minus(2) * pi/180;
actual_minus_temp(3) = actual_minus(3);

% weird Matlab azimuth convention
if zero_temp(1) > pi
    zero_temp(1) = zero_temp(1) - 2*pi;
end
if actual_plus_temp(1) > pi
    actual_plus_temp(1) = actual_plus_temp(1) - 2*pi;
end
if actual_minus_temp(1) > pi
    actual_minus_temp(1) = actual_minus_temp(1) - 2*pi;
end

% Diff between the two plane-projections of actuals, and between each of them and zero,
diff_proj_plusminus = acos(sin(proj_plus_sph(2))*sin(proj_minus_sph(2)) + cos(proj_minus_sph(2))*sin(proj_minus_sph(1))*...
                cos(proj_plus_sph(2))*sin(proj_plus_sph(1)) + cos(proj_minus_sph(2))*cos(proj_minus_sph(1))*...
                cos(proj_plus_sph(2))*cos(proj_plus_sph(1)));
diff_proj_plus = acos(sin(proj_plus_sph(2))*sin(zero_temp(2)) + cos(zero_temp(2))*sin(zero_temp(1))*...
                cos(proj_plus_sph(2))*sin(proj_plus_sph(1)) + cos(zero_temp(2))*cos(zero_temp(1))*...
                cos(proj_plus_sph(2))*cos(proj_plus_sph(1)));
diff_proj_minus = acos(sin(proj_minus_sph(2))*sin(zero_temp(2)) + cos(zero_temp(2))*sin(zero_temp(1))*...
                cos(proj_minus_sph(2))*sin(proj_minus_sph(1)) + cos(zero_temp(2))*cos(zero_temp(1))*...
                cos(proj_minus_sph(2))*cos(proj_minus_sph(1)));

% ...between the two expecteds, and between each of them and zero (not always alph! [see m3c148r1]),
diff_exp_plusminus = acos(sin(exp_plus_sph(2))*sin(exp_minus_sph(2)) + cos(exp_minus_sph(2))*sin(exp_minus_sph(1))*...
                cos(exp_plus_sph(2))*sin(exp_plus_sph(1)) + cos(exp_minus_sph(2))*cos(exp_minus_sph(1))*...
                cos(exp_plus_sph(2))*cos(exp_plus_sph(1)));
diff_exp_plus = acos(sin(exp_plus_sph(2))*sin(zero_temp(2)) + cos(zero_temp(2))*sin(zero_temp(1))*...
                cos(exp_plus_sph(2))*sin(exp_plus_sph(1)) + cos(zero_temp(2))*cos(zero_temp(1))*...
                cos(exp_plus_sph(2))*cos(exp_plus_sph(1)));
diff_exp_minus = acos(sin(exp_minus_sph(2))*sin(zero_temp(2)) + cos(zero_temp(2))*sin(zero_temp(1))*...
                cos(exp_minus_sph(2))*sin(exp_minus_sph(1)) + cos(zero_temp(2))*cos(zero_temp(1))*...
                cos(exp_minus_sph(2))*cos(exp_minus_sph(1)));

% ...and between actuals (projs) and expecteds.
diff_plus = acos(sin(proj_plus_sph(2))*sin(exp_plus_sph(2)) + cos(exp_plus_sph(2))*sin(exp_plus_sph(1))*...
            cos(proj_plus_sph(2))*sin(proj_plus_sph(1)) + cos(exp_plus_sph(2))*cos(exp_plus_sph(1))*...
            cos(proj_plus_sph(2))*cos(proj_plus_sph(1)));
diff_minus = acos(sin(proj_minus_sph(2))*sin(exp_minus_sph(2)) + cos(exp_minus_sph(2))*sin(exp_minus_sph(1))*...
             cos(proj_minus_sph(2))*sin(proj_minus_sph(1)) + cos(exp_minus_sph(2))*cos(exp_minus_sph(1))*...
             cos(proj_minus_sph(2))*cos(proj_minus_sph(1)));

% Also need angle between proj_plus and exp_minus, and vice versa, for +/- assignment
diff_plusminus = acos(sin(proj_plus_sph(2))*sin(exp_minus_sph(2)) + cos(exp_minus_sph(2))*sin(exp_minus_sph(1))*...
            cos(proj_plus_sph(2))*sin(proj_plus_sph(1)) + cos(exp_minus_sph(2))*cos(exp_minus_sph(1))*...
            cos(proj_plus_sph(2))*cos(proj_plus_sph(1)));
diff_minusplus = acos(sin(proj_minus_sph(2))*sin(exp_plus_sph(2)) + cos(exp_plus_sph(2))*sin(exp_plus_sph(1))*...
             cos(proj_minus_sph(2))*sin(proj_minus_sph(1)) + cos(exp_plus_sph(2))*cos(exp_plus_sph(1))*...
             cos(proj_minus_sph(2))*cos(proj_minus_sph(1)));

% +/- assignment: if the projs are closer to the *opposite*
% expected than to the correct one, then it's a negative shift ratio
if diff_plusminus < diff_plus
    diff_proj_plus = -diff_proj_plus;
end
if diff_minusplus < diff_minus
    diff_proj_minus = -diff_proj_minus;
end
% for the 'plusminus' ratio the logic is more complicated (see notes)
if (diff_plusminus < diff_plus & diff_minusplus < diff_minus)
    diff_proj_plusminus = -diff_proj_plusminus;
end
if (diff_plusminus > diff_plus & diff_minusplus < diff_minus & diff_proj_plus < diff_proj_minus)
    diff_proj_plusminus = -diff_proj_plusminus;
end
if (diff_plusminus < diff_plus & diff_minusplus > diff_minus & diff_proj_minus < diff_proj_plus)
    diff_proj_plusminus = -diff_proj_plusminus;
end
    
% Now to calculate 'noise ratio' (vertical shift ratio for horiz, 
% and vice versa), need diff angles between actual and proj:
diff_noise_plus = acos(sin(actual_plus_temp(2))*sin(proj_plus_sph(2)) + cos(proj_plus_sph(2))*sin(proj_plus_sph(1))*...
            cos(actual_plus_temp(2))*sin(actual_plus_temp(1)) + cos(proj_plus_sph(2))*cos(proj_plus_sph(1))*...
            cos(actual_plus_temp(2))*cos(actual_plus_temp(1)));
diff_noise_minus = acos(sin(actual_minus_temp(2))*sin(proj_minus_sph(2)) + cos(proj_minus_sph(2))*sin(proj_minus_sph(1))*...
             cos(actual_minus_temp(2))*sin(actual_minus_temp(1)) + cos(proj_minus_sph(2))*cos(proj_minus_sph(1))*...
             cos(actual_minus_temp(2))*cos(actual_minus_temp(1)));

% For plusminus, let's try the difference angle between the PLANE defined by the actuals
% and the plane defined by the expecteds (same as diff between NORMALS)
[normal_pm_sph(1),normal_pm_sph(2),normal_pm_sph(3)] = cart2sph(normal_plusminus(1),normal_plusminus(2),normal_plusminus(3));
normal_actuals = cross([actual_plus_x,actual_plus_y,actual_plus_z],[actual_minus_x,actual_minus_y,actual_minus_z]);
[normal_act_sph(1),normal_act_sph(2),normal_act_sph(3)] = cart2sph(normal_actuals(1),normal_actuals(2),normal_actuals(3));

diff_noise_plusminus = acos(sin(normal_pm_sph(2))*sin(normal_act_sph(2)) + cos(normal_act_sph(2))*sin(normal_act_sph(1))*...
            cos(normal_pm_sph(2))*sin(normal_pm_sph(1)) + cos(normal_act_sph(2))*cos(normal_act_sph(1))*...
            cos(normal_pm_sph(2))*cos(normal_pm_sph(1)));

% assigning +/- for noise ratio is arbitrary, yet even harder than for SR
% At least for horiz it's easy: compare the elevations of the actual to
% its proj, and if the actual is lower then the NR is negative
if gaze_dir == 0
    if actual_plus_temp(2) < proj_plus_sph(2)
        diff_noise_plus = -diff_noise_plus;
    end
    if actual_minus_temp(2) < proj_minus_sph(2)
        diff_noise_minus = -diff_noise_minus;
    end
    % and for plusminus, just compare the elevations of actuals
    if actual_plus_temp(2) < actual_minus_temp(2)
        diff_noise_plusminus = -diff_noise_plusminus;
    end
else % for vert, need to rotate vectors to align with the Y axis
     % (so I first need to rotate the zero vector around the Z axis by the amount that its
     % azimuth differs from 90, then down (around the X axis) so that it aligns with the Y)
    if zero(1) < 90
        azi_rot = (zero(1) + 270) * pi/180;
    else
        azi_rot = (zero(1) - 90) * pi/180;
    end
    ele_rot = zero(2) * pi/180;
    
    rotate_Z = [cos(azi_rot),sin(azi_rot),0; -sin(azi_rot),cos(azi_rot),0; 0,0,1];
    rotate_X = [1,0,0; 0,cos(ele_rot),-sin(ele_rot); 0,sin(ele_rot),cos(ele_rot)];

    rotated_zero = rotate_X * (rotate_Z * [zero_x;zero_y;zero_z]);
    rotated_actual_plus = rotate_X * (rotate_Z * [actual_plus_x;actual_plus_y;actual_plus_z]);
    rotated_actual_minus = rotate_X * (rotate_Z * [actual_minus_x;actual_minus_y;actual_minus_z]);
    rotated_exp_plus = rotate_X * (rotate_Z * expected_plus');
    rotated_exp_minus = rotate_X * (rotate_Z * expected_minus');
    
    % Now the 'direction' of the noise should be the component of the
    % actual along the Z axis (arbitrary, but consistent within cells)
    
    % (turns out this works better - not sure why)
    if rotated_exp_plus(3) - rotated_actual_plus(3) < 0
        diff_noise_plus = -diff_noise_plus;
    end
    if rotated_exp_minus(3) - rotated_actual_minus(3) < 0
        diff_noise_minus = -diff_noise_minus;
    end
    % and for plusminus (should work, but who knows):
    if rotated_actual_plus(3) < rotated_actual_minus(3)
        diff_noise_plusminus = -diff_noise_plusminus;
    end    
end

% convert back to degrees
exp_plus_sph(1) = exp_plus_sph(1) * 180/pi;
exp_plus_sph(2) = exp_plus_sph(2) * 180/pi;
exp_minus_sph(1) = exp_minus_sph(1) * 180/pi;
exp_minus_sph(2) = exp_minus_sph(2) * 180/pi;
proj_plus_sph(1) = proj_plus_sph(1) * 180/pi;
proj_plus_sph(2) = proj_plus_sph(2) * 180/pi;
proj_minus_sph(1) = proj_minus_sph(1) * 180/pi;
proj_minus_sph(2) = proj_minus_sph(2) * 180/pi;
diff_plus = diff_plus * 180/pi;
diff_minus = diff_minus * 180/pi;
diff_plusminus = diff_plusminus * 180/pi;
diff_minusplus = diff_minusplus * 180/pi;
diff_proj_plusminus = diff_proj_plusminus * 180/pi;
diff_proj_plus = diff_proj_plus * 180/pi;
diff_proj_minus = diff_proj_minus * 180/pi;
diff_exp_plusminus = diff_exp_plusminus * 180/pi;
diff_exp_plus = diff_exp_plus * 180/pi;
diff_exp_minus = diff_exp_minus * 180/pi;
diff_noise_plusminus = diff_noise_plusminus * 180/pi;
diff_noise_plus = diff_noise_plus * 180/pi;
diff_noise_minus = diff_noise_minus * 180/pi;

% (if negative, add 360 so the range of azimuths is 0 to 360)
if exp_plus_sph(1) < 0
    exp_plus_sph(1) = exp_plus_sph(1) + 360;
end
if exp_minus_sph(1) < 0
    exp_minus_sph(1) = exp_minus_sph(1) + 360;
end
if proj_plus_sph(1) < 0
    proj_plus_sph(1) = proj_plus_sph(1) + 360;
end
if proj_minus_sph(1) < 0
    proj_minus_sph(1) = proj_minus_sph(1) + 360;
end

% and finally, the shift(noise) ratio is the diff_proj(noise) / diff_exp
shift_ratio_plusminus = diff_proj_plusminus / diff_exp_plusminus;
shift_ratio_plus = diff_proj_plus / diff_exp_plus;
shift_ratio_minus = diff_proj_minus / diff_exp_minus;

noise_ratio_plusminus = diff_noise_plusminus / diff_exp_plusminus;
noise_ratio_plus = diff_noise_plus / diff_exp_plus;
noise_ratio_minus = diff_noise_minus / diff_exp_plus;

%***May want to change the noise ratio for plusminus:
%currently it's the difference angle between the planes, but
%for head-centered cells this can become very large, with only minor
%deviations of the actuals from the plane of expecteds (see m3c139r1)****

% ***** End shift ratio calculation *****
%______________________________________________________________


% (reverse the elevations again)
zero(2) = -zero(2);
exp_plus_sph(2) = -exp_plus_sph(2);
exp_minus_sph(2) = -exp_minus_sph(2);
actual_plus(2) = -actual_plus(2);
actual_minus(2) = -actual_minus(2);
proj_plus_sph(2) = -proj_plus_sph(2);
proj_minus_sph(2) = -proj_minus_sph(2);
