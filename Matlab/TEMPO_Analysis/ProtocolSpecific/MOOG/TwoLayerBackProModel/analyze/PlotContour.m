% plot contour map from raw output of hidden units ---GY, 05/12/05
%--------------------------------------------------------------

stim_type = 1;  % range 1,2,3, representing ves,vis,comb, respectively
gaze_angle_h = 1; % check unique_gaze_angle_xr and check unique_gaze_angle_yr
unit_num = 127; % can select any unit from 1~150
resp_mat=[];

for n = 127 : 127
	figure(2);
	set(2,'Position', [5,15 1000,650], 'Name', num2str(n));
    
	for s = 1:2
	%	for n = 1 : 150
			for k = 1 : 3
				for j = 1 : 5  % elevation arraged as [-90,-45,0,45,90]
                    for i = 1 : 9 % azimuth arraged as [270,225,180,135,90,45,0,315,270]
                        if (270-(i-1)*45) >= 0 
                            select = find (unique_point_azimuth == (270-(i-1)*45) & unique_point_elevation == (-90+(j-1)*45) );
                        else 
                            select = find (unique_point_azimuth == (360+270-(i-1)*45) & unique_point_elevation == (-90+(j-1)*45) );
                        end
                        if select > 0
                           resp_mat{s,k}(j, i) = hidden_outputs(s,k*2-1,select,n);
                       else
                           select = find (unique_point_azimuth == 0 & unique_point_elevation == (-90+(j-1)*45) );
                           resp_mat{s,k}(j, i) = hidden_outputs(s,k*2-1,select,n);
                       end
                    end
				end
                resp_mat{s,k}(:, :) = resp_mat{s,k}(:, :) - min(min(resp_mat{s,k}(:,:))) ; % make output range from 0~2
            end
            min_resp(1,s) = min([min(min(resp_mat{s,1}(:,:))), min(min(resp_mat{s,2}(:,:))),min(min(resp_mat{s,3}(:,:)))]); 
            max_resp(1,s) = max([max(max(resp_mat{s,1}(:,:))), max(max(resp_mat{s,2}(:,:))),max(max(resp_mat{s,3}(:,:)))]); 
			% plot gain fields for all units, presented as line
	% 		minn1= min(min(resp_mat{1}(:,:)));
	% 		minn2= min(min(resp_mat{2}(:,:)));
	% 		minn3= min(min(resp_mat{3}(:,:)));
	%         minn4= min(min(resp_mat{4}(:,:)));
	% 		minn5= min(min(resp_mat{5}(:,:)));
	% 		minn=min([minn1,minn2,minn3,minn4,minn5]);
	% 		maxx1= max(max(resp_mat{1}(:,:)));
	% 		maxx2= max(max(resp_mat{2}(:,:)));
	% 		maxx3= max(max(resp_mat{3}(:,:)));
	%         maxx4= max(max(resp_mat{4}(:,:)));
	% 		maxx5= max(max(resp_mat{5}(:,:)));
	% 		maxx=max([maxx1,maxx2,maxx3,maxx4,maxx5]);
	%         subplot(1,2,1);
	%         plot([1,2], [maxx1-minn1,maxx2-minn2],'-');
	%         hold on;
	%         plot([2,3], [maxx2-minn2,maxx3-minn3],'-');
	%         hold on;
	%         plot([3,4], [maxx3-minn3,maxx4-minn4],'-');
	%         hold on;
	%         plot([4,5], [maxx4-minn4,maxx5-minn5],'-');
	%         hold on;
	%         set(gca,'xticklabel',[-40,-20,0,20,40]);
	%         ylim([0,2]);
	%         subplot(1,2,2);
	%         plot([1,2], [maxx1,maxx2],'-');
	%         hold on;
	%         plot([2,3], [maxx2,maxx3],'-');
	%         hold on;
	%         plot([3,4], [maxx3,maxx4],'-');
	%         hold on;
	%         plot([4,5], [maxx4,maxx5],'-');
	%         hold on;
	%         set(gca,'xticklabel',[-40,-20,0,20,40]);
	%         ylim([0,2]);    
	%         slope_am(s,n)= abs( maxx5-minn5 - (maxx1-minn1) );
	%         slope(s,n)= abs( maxx5-maxx1 );
	%	end
	end
    azi_cos = [1,2,3,4,5,6,7,8,9];
    ele_sin = [-1,-0.707,0,0.707,1];
	for s = 1 : 2
        for k = 1 : 3
            axes('position',[0.05+(s-1)*0.45 0.05+(k-1)*0.305 0.4 0.29]);
            contourf(azi_cos, ele_sin, resp_mat{s,k});        
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            caxis( [ min(min_resp), max(max_resp) ] );
            colorbar;    
        end
	end
%     dlmwrite('ve-40.txt',resp_mat{1,1});
%     dlmwrite('ve0.txt',resp_mat{1,2});
%     dlmwrite('ve40.txt',resp_mat{1,3});
%     dlmwrite('vi-40.txt',resp_mat{2,1});
%     dlmwrite('vi0.txt',resp_mat{2,2});
%     dlmwrite('vi40.txt',resp_mat{2,3});
%	pause;
%	close;
end




