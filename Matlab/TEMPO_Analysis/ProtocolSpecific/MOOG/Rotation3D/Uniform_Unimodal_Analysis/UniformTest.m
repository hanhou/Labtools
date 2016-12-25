% uniform distibution, based on permutation of uniform
%GY -09/29/2006
% notice: data and range need to be set properly first for your data
clear all

data = 1; % if 3D, set 3 (2 vectors with azimuth and elevation), lambert transform will be forced (cosine transformed)
          % if 2D, set 2 (2 vectors with azimuth), no lambert transform
          % if 1D, set 1 (1 vector)
if data ==1
    % change this according to your data
     range  = 0 : 360;    % circular data
%     range  = 0 : 180;  % uncircular data
%    range  = -1 : 1;   % lambert transform
end

bin_num = 12;  % asign how many bins should be there, this is a test of the bin, results would be a little varied according to how you bin your data

aa = dlmread('ModeTestData.txt');
dim = size(aa);
index_n = 1:1:dim(1);
num_perm =1000;

if data == 3
    V1azi = aa(:,1)*pi/180;
    V1ele = aa(:,2)*pi/180;
    V2azi = aa(:,3)*pi/180;
    V2ele = aa(:,4)*pi/180;

    for i = 1:dim(1)
        diff_angle(i) = (180/3.14159) * acos( sin(V1ele(i)) * sin(V2ele(i)) + cos(V2ele(i)) * sin(V2azi(i)) * cos(V1ele(i)) * sin(V1azi(i)) + cos(V2ele(i)) * cos(V2azi(i)) * cos(V1ele(i)) * cos(V1azi(i)) );
    end
    hist_raw = hist(cos(diff_angle*pi/180),bin_num); % lambert transform    

    R = sum( (hist_raw - mean(hist_raw)).^2 );

    % do permutation now, bin first, and then permute the bin
    for b = 1 : num_perm
        index_perm = index_n( randperm(dim(1)) );
        for i = 1:dim(1)
            diff_angle_perm(i) = (180/3.14159) * acos( sin(V1ele(index_perm(i))) * sin(V2ele(i)) + cos(V2ele(i)) * sin(V2azi(i)) * cos(V1ele(index_perm(i))) * sin(V1azi(index_perm(i))) + cos(V2ele(i)) * cos(V2azi(i)) * cos(V1ele(index_perm(i))) * cos(V1azi(index_perm(i))) );
        end
        hist_raw_perm = hist(cos(diff_angle_perm*pi/180),bin_num); % lambert transform
        R_perm(b) = sum( (hist_raw_perm - mean(hist_raw)).^2 );
        plot(hist_raw_perm,'b-');
        hold on;
    end
    plot(hist_raw,'r-');
    plot([1 bin_num],[mean(hist_raw),mean(hist_raw)],'r-');
    
elseif data ==2 % 2D data
    V1 = aa(:,1);
    V2 = aa(:,2);

    for i = 1:dim(1)
        diff_angle(i) = V1(i)-V2(i);
        if diff_angle(i)>180
            diff_angle(i) = 360- diff_angle(i);
        end
        if diff_angle(i)<-180
            diff_angle(i) = 360+ diff_angle(i);
        end
    end
    hist_raw = hist(diff_angle,bin_num); 

    R = sum( (hist_raw - mean(hist_raw)).^2 );

    % do permutation now, bin first, and then permute the bin
    for b = 1 : num_perm
        index_perm = index_n( randperm(dim(1)) );
        for i = 1:dim(1)
            diff_angle_perm(i) = V1(index_perm(i))-V2(i);
            if diff_angle_perm(i)>180
                diff_angle_perm(i) = 360- diff_angle_perm(i);
            end
            if diff_angle_perm(i)<-180
                diff_angle_perm(i) = 360+ diff_angle_perm(i);
            end
        end
        hist_raw_perm = hist(diff_angle_perm,bin_num); 
        R_perm(b) = sum( (hist_raw_perm - mean(hist_raw)).^2 );
        plot(hist_raw_perm,'b-');
        hold on;
    end
    plot(hist_raw,'r-');
    plot([1 bin_num],[mean(hist_raw),mean(hist_raw)],'r-');
    
else % 1D data
    V1 = aa;   
    hist_raw = hist(V1,bin_num); 

    R = sum( (hist_raw - mean(hist_raw)).^2 );
    
    for b = 1:num_perm
        if min(range) == 0
           diff_angle_perm = rand(1,dim(1))*max(range);        
        elseif min(range) == -1
           diff_angle_perm = rand(1,dim(1))*2-1;
        else
           warning = 'reset data range properly';
        end
        hist_raw_perm = hist(diff_angle_perm,bin_num); 
        R_perm(b) = sum( (hist_raw_perm - mean(hist_raw)).^2 );
        plot(hist_raw_perm,'b-');
        hold on;
    end
    plot(hist_raw,'r-');
    plot([1 bin_num],[mean(hist_raw),mean(hist_raw)],'r-');
    
end

% now calculate p value or significant test
p =length(find(R_perm>R)) / num_perm;
p

