% analyze degree of overlap between two clusters
% each cluster has azimuth and elevation
% check whether the overlapped area is <0.05
% GY 02/02/2007

aa = dlmread('Rotation_fixation.dat','',1,1);  % data during fixation
bb = dlmread('Rotation_darkness.dat','',1,1);  % data in darkness
aa_pre = dlmread('Rotation_fixation.txt'); % preferred direction during fication
bb_pre = dlmread('Rotation_darkness.txt'); % preferred direction during darkness
dim = size(aa);

for i = 1:dim(1)    
    fix_azi(i,:) = aa(i,1:1000)*pi/180;
    fix_ele(i,:) = aa(i,1001:2000)*pi/180;
    fix_azi_pre(i) = aa_pre(i,1)*pi/180;
    fix_ele_pre(i) = aa_pre(i,2)*pi/180;
    
    dark_azi(i,:) = bb(i,1:1000)*pi/180;
    dark_ele(i,:) = bb(i,1001:2000)*pi/180;
    dark_azi_pre(i) = bb_pre(i,1)*pi/180;
    dark_ele_pre(i) = bb_pre(i,2)*pi/180;  
    
    for j = 1:1000
        diff_fix(i,j) = (180/3.14159) * acos( sin(fix_ele_pre(i)) * sin(fix_ele(i,j)) + cos(fix_ele(i,j)) * sin(fix_azi(i,j)) * cos(fix_ele_pre(i)) * sin(fix_azi_pre(i)) + cos(fix_ele(i,j)) * cos(fix_azi(i,j)) * cos(fix_ele_pre(i)) * cos(fix_azi_pre(i)) );
        diff_dark(i,j) = (180/3.14159) * acos( sin(dark_ele_pre(i)) * sin(dark_ele(i,j)) + cos(dark_ele(i,j)) * sin(dark_azi(i,j)) * cos(dark_ele_pre(i)) * sin(dark_azi_pre(i)) + cos(dark_ele(i,j)) * cos(dark_azi(i,j)) * cos(dark_ele_pre(i)) * cos(dark_azi_pre(i)) );
    end
    diff_fix_dark(i) = (180/3.14159) * acos( sin(fix_ele_pre(i)) * sin(dark_ele_pre(i)) + cos(dark_ele_pre(i)) * sin(dark_azi_pre(i)) * cos(fix_ele_pre(i)) * sin(fix_azi_pre(i)) + cos(dark_ele_pre(i)) * cos(dark_azi_pre(i)) * cos(fix_ele_pre(i)) * cos(fix_azi_pre(i)) );

    % judge how much data points are overlapped
   diffsum = diff_fix(i,:)+diff_dark(i,:);
   p(i) = length(find(diffsum>diff_fix_dark(i) )) / 1000; % the total number is 1000 data points   

end


