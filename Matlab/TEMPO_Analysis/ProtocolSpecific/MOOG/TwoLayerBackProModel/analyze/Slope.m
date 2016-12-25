% calculate slope of gain field  ---GY, 05/12/05
%--------------------------------------------------------------
global maxx;
resp_mat=[];
figure(2);
set(2,'Position', [5,15 1000,650], 'Name', 'gain field');
for s = 1:2
	for n = 1 : 150
		for k = 1 : 5
			for j = 1 : 5  % elevation arraged as [-90,-45,0,45,90]
                for i = 1 : 9 % azimuth arraged as [270,225,180,135,90,45,0,315,270]
                    if (270-(i-1)*45) >= 0 
                        select = find (unique_point_azimuth == (270-(i-1)*45) & unique_point_elevation == (-90+(j-1)*45) );
                    else 
                        select = find (unique_point_azimuth == (360+270-(i-1)*45) & unique_point_elevation == (-90+(j-1)*45) );
                    end
                    if select > 0
                       resp_mat{k}(j, i) = hidden_outputs(s,k,select,n);
                   else
                       select = find (unique_point_azimuth == 0 & unique_point_elevation == (-90+(j-1)*45) );
                       resp_mat{k}(j, i) = hidden_outputs(s,k,select,n);
                   end
                end
			end
            resp_mat{k}(:, :) = resp_mat{k}(:, :) - min(min(resp_mat{k}(:,:))) ; % make output range from 0~2
        end
		% plot gain fields for all units, presented as line
		minn1= min(min(resp_mat{1}(:,:)));
		minn2= min(min(resp_mat{2}(:,:)));
		minn3= min(min(resp_mat{3}(:,:)));
        minn4= min(min(resp_mat{4}(:,:)));
		minn5= min(min(resp_mat{5}(:,:)));
		minn=[minn1,minn2,minn3,minn4,minn5];
		maxx1= max(max(resp_mat{1}(:,:)));
		maxx2= max(max(resp_mat{2}(:,:)));
		maxx3= max(max(resp_mat{3}(:,:)));
        maxx4= max(max(resp_mat{4}(:,:)));
		maxx5= max(max(resp_mat{5}(:,:)));
		maxx=[maxx1,maxx2,maxx3,maxx4,maxx5];
        subplot(1,2,s);
        plot(maxx,'-');
        set(gca,'xticklabel',[-40,-20,0,20,40]);
        ylim([0,2]);
        hold on;
        % now calculate slope in two ways
        % first, calculate the difference between points directly
        diff1 = abs(maxx3-maxx1)/maxx1; 
        diff2 = abs(maxx5-maxx1)/maxx1;
        slope_d(s,n)= (diff1+diff2)/2;        
        % second, fit a linear line and calculate the slope
        gaze_x = [-40,-20,0,20,40];
        % estimate first 
        estimate(1,1) = (maxx4-maxx2)/40;
        estimate(1,2) = (maxx4+maxx2)/2;
        quick = fminsearch('slope_func',estimate);
        slope_f(s,n) = quick(1);
        % this is to seperate three types of slope: flat, concave/convex
        % and monotonically change        
        slope_left = regress([maxx1,maxx2,maxx3]',[-20,0,20]'); % here only use the y error to fit linear line
        slope_right = regress([maxx3,maxx4,maxx5]',[-20,0,20]'); % here only use the y error to fit linear line
        slope_g(s,n) = slope_left*slope_right; % data are seperated into 3 groups 
%         if s==1
%             subplot(1,2,1);
%             plot([1,2], [maxx1-minn1,maxx2-minn2],'-');
%             hold on;
%             plot([2,3], [maxx2-minn2,maxx3-minn3],'-');
%             hold on;
%             plot([3,4], [maxx3-minn3,maxx4-minn4],'-');
%             hold on;
%             plot([4,5], [maxx4-minn4,maxx5-minn5],'-');
%             hold on;
%             set(gca,'xticklabel',[-40,-20,0,20,40]);
%             ylim([0,2]);
%         else
%             subplot(1,2,2);
%             plot([1,2], [maxx1-minn1,maxx2-minn2],'-');
%             hold on;
%             plot([2,3], [maxx2-minn2,maxx3-minn3],'-');
%             hold on;
%             plot([3,4], [maxx3-minn3,maxx4-minn4],'-');
%             hold on;
%             plot([4,5], [maxx4-minn4,maxx5-minn5],'-');
%             hold on;
%             set(gca,'xticklabel',[-40,-20,0,20,40]);
%             ylim([0,2]);
%         end
    end
end
% dlmwrite('slope_d.txt',slope_d');
dlmwrite('slope_f.txt',slope_f');
%dlmwrite('slope_g.txt',slope_g');

