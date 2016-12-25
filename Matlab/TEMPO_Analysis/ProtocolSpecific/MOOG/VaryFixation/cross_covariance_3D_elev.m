% CROSS_COVARIANCE_3D_elev.m -- compute 'orthogonal DI' by shifting 3D tuning 
% curves in elevation instead of azimuth.

function DI = cross_covariance_3D_elev(resp_mat, unique_stim_type, unique_condition_num, bin, method, Vec_sum, FILE);
tic

shift = 40;  % will rotate tuning surface in elevation by this amount, both positive and negative
x = [0 45 90 135 180 225 270 315 360];
xi = 0 : bin : 360-bin;
y = [90; 45; 0; -45; -90];
yi = flipud((-90:bin:90)');

for i = 1:length(xi)
    ele(:,i) = yi;
end

for h = 1:length(yi)
    azi(h,:) = xi;
end

for n = 1:length(unique_stim_type)
    for k = 1:length(unique_condition_num)

        resp_mat_360{k,n} = squeeze(resp_mat(k+3*(n-1),:,:));   % kluge, to make the resp_mat more intuitive
        resp_mat_360{k,n}(:,end+1) = resp_mat_360{k,n}(:,1);
        resp_mat_360{k,n}(1,:) = resp_mat_360{k,n}(1,1);
        resp_mat_360{k,n}(end,:) = resp_mat_360{k,n}(end,1);
        
        resp_mat_360{k,n} = flipud(resp_mat_360{k,n});  % because of our weird ele sign convention
        
        for j = 1:length(y)
            clear temp_y;
            temp_y = resp_mat_360{k,n}(j,:);
            resp_mat_temp{k,n}(j,:) = interp1(x, temp_y, xi, method); % interpolate in azimuth
        end
        
        for i = 1:length(xi)
            clear temp_x;
            temp_x = resp_mat_temp{k,n}(:,i);
            resp_mat_interp{k,n}(:,i) = interp1(y, temp_x, yi, method); % interpolate in elevation
        end
        
        
        for t = 0:shift/bin                                                                 % for each rotation step:
%         for t = 0:2  %TEMP% - to graph it, for testing
            FILE
            n
            k
            t
            clear vecsum rot_vect xx yy zz rotate_plus rotate_minus rot_plus* rot_minus*;
            alph = t*bin*pi/180;
%             alph = 15*t*pi/180;  %TEMP% - to graph it, for testing

            % rotate_plus = [1,0,0; 0,cos(alph),-sin(alph); 0,sin(alph),cos(alph)];           % 1) define rotation matrices
%   WRONG!  % rotate_minus = [1,0,0; 0,cos(-alph),-sin(-alph); 0,sin(-alph),cos(-alph)];

%---------------------------------------------------------------------------------------------            
%       Instead, rotate around an arbitrary axis defined as the vector in the horizontal 
%       plane that is normal to the preferred direction vector.  This should rotate 
%       the tuning in elevation only (at least for well-behaved tuning functions).**
%---------------------------------------------------------------------------------------------
            % First find the correct rotation axis (unit vector) by taking the cross product of 
            % the preferred heading vector and the unit vector normal to the horizontal plane:
            [vecsum(1) vecsum(2) vecsum(3)] = sph2cart(Vec_sum{k+3*(n-1)}(1)*pi/180, -Vec_sum{k+3*(n-1)}(2)*pi/180, 1);
            rot_vect = cross(vecsum, [0 0 1]);
            xx = rot_vect(1); yy = rot_vect(2); zz = rot_vect(3);
            % Then the new rotation matrix is: (see http://www.euclideanspace.com/maths/algebra/matrix/orthogonal/rotation/):
            rotate_plus = [1+(1-cos(alph))*(xx*xx-1), -zz*sin(alph)+(1-cos(alph))*xx*yy, yy*sin(alph)+(1-cos(alph))*xx*zz ;
                           zz*sin(alph)+(1-cos(alph))*xx*yy, 1+(1-cos(alph))*(yy*yy-1), -xx*sin(alph)+(1-cos(alph))*yy*zz ;
                          -yy*sin(alph)+(1-cos(alph))*xx*zz, xx*sin(alph)+(1-cos(alph))*yy*zz, 1+(1-cos(alph))*(zz*zz-1) ];
            xx = -xx; yy = -yy; zz = -zz;
            rotate_minus = [1+(1-cos(alph))*(xx*xx-1), -zz*sin(alph)+(1-cos(alph))*xx*yy, yy*sin(alph)+(1-cos(alph))*xx*zz ;
                            zz*sin(alph)+(1-cos(alph))*xx*yy, 1+(1-cos(alph))*(yy*yy-1), -xx*sin(alph)+(1-cos(alph))*yy*zz ;
                           -yy*sin(alph)+(1-cos(alph))*xx*zz, xx*sin(alph)+(1-cos(alph))*yy*zz, 1+(1-cos(alph))*(zz*zz-1) ];
%---------------------------------------------------------------------------------------------

            [X Y Z] = sph2cart(azi*pi/180,ele*pi/180,1);                                    % 2) convert to cartesian
            for i = 1:length(xi)
                for h = 1:length(yi)
                    clear temp_plus temp_minus;
                    temp_plus = rotate_plus * [X(h,i);Y(h,i);Z(h,i)];                       % 3) rotate
                    temp_minus = rotate_minus * [X(h,i);Y(h,i);Z(h,i)];
                    rot_plus_X(h,i) = temp_plus(1); rot_minus_X(h,i) = temp_minus(1);
                    rot_plus_Y(h,i) = temp_plus(2); rot_minus_Y(h,i) = temp_minus(2);
                    rot_plus_Z(h,i) = temp_plus(3); rot_minus_Z(h,i) = temp_minus(3);       
                end
            end
                                                                                            % 4) convert back to spherical:
            [rot_plus_azi, rot_plus_ele, rot_plus_amp] = cart2sph(rot_plus_X, rot_plus_Y, rot_plus_Z);
            [rot_minus_azi, rot_minus_ele, rot_minus_amp] = cart2sph(rot_minus_X, rot_minus_Y, rot_minus_Z);

            rot_plus_azi = bin * round((rot_plus_azi*180/pi)/bin);                          % 5) round to nearest binth degree
            rot_plus_azi(rot_plus_azi < 0) = rot_plus_azi(rot_plus_azi < 0) + 360;          % (and add 360 to negative azimuths)
            rot_plus_ele = bin * round((rot_plus_ele*180/pi)/bin);
            rot_minus_azi = bin * round((rot_minus_azi*180/pi)/bin);
            rot_minus_azi(rot_minus_azi < 0) = rot_minus_azi(rot_minus_azi < 0) + 360;
            rot_minus_ele = bin * round((rot_minus_ele*180/pi)/bin);
            for i = 1:length(xi)
                for h = 1:length(yi)
                    select_plus = logical(azi==rot_plus_azi(h,i) & ele==rot_plus_ele(h,i)); % 6) find where that azi/ele pair exists
                    select_minus = logical(azi==rot_minus_azi(h,i) & ele==rot_minus_ele(h,i));
                    if sum(sum(select_plus)) == 1 & sum(sum(select_minus)) == 1
                        resp_mat_rot_plus{t+1,k,n}(h,i) = resp_mat_interp{k,n}(select_plus);% 7) and assign the firing rate value from that location
                        resp_mat_rot_minus{t+1,k,n}(h,i) = resp_mat_interp{k,n}(select_minus);
                    elseif sum(sum(select_plus)) == 0 | sum(sum(select_minus)) == 0
                        disp('wtf? sum = 0');
                        return;
                    else
                        disp('wtf? sum > 1');
                        return;
                    end
                end
            end
%             figure(t+10); contourf(resp_mat_rot_plus{t+1,k,n}); %TEMP% - to graph it, for testing
%             figure(t+20); contourf(resp_mat_rot_minus{t+1,k,n});
        end
        
    end
    
    % Now find covariance** between each pair of tuning surfaces (measured at different gaze angles)
    % **(CORRECTION: here I need to use corrcoef instead of cov, because the overall column variances have 
    % changed due to the rotation.  Corrcoef normalizes by the column variance, so it gives the correct answer.)
    %______________________________________________________________________
    % Zero with Itself as a test:
    clear cov_self;
    for t = 0:shift/bin
        clear cov_temp;
        cov_temp = corrcoef(resp_mat_rot_plus{t+1,2,n}, resp_mat_interp{2,n});
        cov_self(shift/bin+1+t) = cov_temp(1,2);
        
        clear cov_temp;
        cov_temp = corrcoef(resp_mat_rot_minus{t+1,2,n}, resp_mat_interp{2,n});
        cov_self(shift/bin+1-t) = cov_temp(1,2);
    end
    DI_index = find(cov_self == max(cov_self));
    DI_self(n) = (DI_index - (shift/bin) - 1) * bin / unique_condition_num(3)
    %______________________________________________________________________
    % Minus with Zero:
    clear cov_MZ;
    for t = 0:shift/bin
        clear cov_temp;
        cov_temp = corrcoef(resp_mat_rot_plus{t+1,1,n}, resp_mat_interp{2,n});
        cov_MZ(shift/bin+1+t) = cov_temp(1,2);
        
        clear cov_temp;
        cov_temp = corrcoef(resp_mat_rot_minus{t+1,1,n}, resp_mat_interp{2,n});
        cov_MZ(shift/bin+1-t) = cov_temp(1,2);
    end
    DI_index = find(cov_MZ == max(cov_MZ));
    DI(1,n) = (DI_index - (shift/bin) - 1) * bin / unique_condition_num(3);
    %______________________________________________________________________
    % Minus with Plus
    clear cov_MP;
    for t = 0:shift/bin
        clear cov_temp;
        cov_temp = corrcoef(resp_mat_rot_plus{t+1,1,n}, resp_mat_interp{3,n});
        cov_MP(shift/bin+1+t) = cov_temp(1,2);
        
        clear cov_temp;
        cov_temp = corrcoef(resp_mat_rot_minus{t+1,1,n}, resp_mat_interp{3,n});
        cov_MP(shift/bin+1-t) = cov_temp(1,2);
    end
    DI_index = find(cov_MP == max(cov_MP));
    DI(2,n) = (DI_index - (shift/bin) - 1) * bin / (unique_condition_num(3)*2);
    %______________________________________________________________________
    % Zero with Plus    
    clear cov_ZP;
    for t = 0:shift/bin
        clear cov_temp;
        cov_temp = corrcoef(resp_mat_rot_plus{t+1,2,n}, resp_mat_interp{3,n});
        cov_ZP(shift/bin+1+t) = cov_temp(1,2);
        
        clear cov_temp;
        cov_temp = corrcoef(resp_mat_rot_minus{t+1,2,n}, resp_mat_interp{3,n});
        cov_ZP(shift/bin+1-t) = cov_temp(1,2);
    end
    DI_index = find(cov_ZP == max(cov_ZP));
    DI(3,n) = (DI_index - (shift/bin) - 1) * bin / unique_condition_num(3);
    
end
toc
return;

% **Here's the reasoning behind this: a rotation about the vertical
% axis (i.e., for regular DI) only changes the azimuth of each point.
% However, a rotation about the horizontal (interaural or interocular)
% axis changes elevation AND azimuth in a complex way.  So what happens is
% you get a 'diagonal' shift of the tuning peak, as viewed from a 2D
% projection.  This produces an artificially high orthogonal DI for the visual
% condition, because the FALSE shift in azimuth (due to this 'vertical'
% rotation) causes the tuning to be best matched with TRUE shift in azimuth
% (the eye-centeredness).
% Instead, we do our best to rotate the business end of the tuning curve 
% PRIMARILY in elevation, by choosing the rotation axis to be orthogonal to 
% the preferred direction vector (strictly speaking, a rotation that only 
% changes the elevation of a vector is impossible for all vectors other than 
% the preferred direction)**