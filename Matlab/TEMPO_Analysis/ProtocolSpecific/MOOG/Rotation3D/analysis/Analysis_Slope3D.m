% analyze slope for the 3D tuning curve in space

% choose protocol first
%analyze_protocol = 'rotation';
analyze_protocol = 'translat';

if analyze_protocol == 'rotation';
	aa = dlmread('mean3D_rot.dat','',0,1);  % load data
	bb = dlmread('std3D_rot.txt');
	cc = dlmread('p_rot.txt');
elseif analyze_protocol == 'translat';
	aa = dlmread('mean3D_tran.dat','',0,1);  % load data
	bb = dlmread('std3D_tran.txt');
	cc = dlmread('p_tran.txt');
end
p_ves = cc(:,1);
p_vis = cc(:,2);
dim = size(aa); %dim(1) should be the number of cells

%transform data into 3D matrix, notice that since data is circular in both
%azimuth and elevation direction, we need to unwrap the data
% extend map show as bellow (azimuth/elevation)

% 135/-45, 180/-45, 225/-45, 270/-45, 315/-45, 0/-45,   45/-45,  90/-45,  135/-45, 180/-45, 225/-45
% 0/-90,   0/-90,   0/-90,   0/-90,   0/-90,   0/-90,   0/-90,   0/-90,   0/-90,   0/-90,   0/-90
% 315/-45, 0/-45,   45/-45,  90/-45,  135/-45, 180/-45, 225/-45, 270/-45, 315/-45, 0/-45,   45/-45
% 315/0,   0/0,     45/0,    90/0,    135/0,   180/0,   225/0,   270/0,   315/0,   0/0,     45/0
% 315/45,  0/45,    45/45,   90/45,   135/45,  180/45,  225/45,  270/45,  315/45,  0/45,    45/45
% 0/90,    0/90,    0/90,    0/90,    0/90,    0/90,    0/90,    0/90,    0/90,    0/90,    0/90
% 135/45,  180/45,  225/45,  270/45,  315/45,  0/45,    45/45,   90/45,   135/45,  180/45,  225/45

x=[-45,0,45,90,135,180,225,270,315,360,405];  % -45,360 and 405 were added to extend data
y=[-135,-90,-45,0,45,90,135]; %-135 and 135 were added to extend data
[xx,yy]=meshgrid(x,y);
xs = -45:1:405;
ys = -135:1:135;
[xxs,yys]=meshgrid(xs,ys);

step = 5; % 5 degree
slide_size = 1;
fisher_ves_azi = zeros(271,360);
fisher_ves_ele = zeros(181,451);
fisher_vis_azi = zeros(271,360);
fisher_vis_ele = zeros(181,451);
n_ves = 0;
n_vis = 0;

for i = 84 : 84
    % extract data first from regular matrix
    ves(1,1:8) = aa(i,1);
    ves(2,:) = aa(i,2:9);
    ves(3,:) = aa(i,10:17);
    ves(4,:) = aa(i,18:25);
    ves(5,:) = aa(i,26);
    vis(1,1:8) = aa(i,27);
    vis(2,:) = aa(i,28:35);
    vis(3,:) = aa(i,36:43);
    vis(4,:) = aa(i,44:51);
    vis(5,:) = aa(i,52);
    % now extend data for later spline fit use
    ves_ext{i}(1,:) = [ves(2,4:8),ves(2,1:5), ves(2,6)];
    ves_ext{i}(2,:) = ves(1,1);
    ves_ext{i}(3,:) = [ves(2,8), ves(2,1:8), ves(2,1:2)];
    ves_ext{i}(4,:) = [ves(3,8), ves(3,1:8), ves(3,1:2)];
    ves_ext{i}(5,:) = [ves(4,8), ves(4,1:8), ves(4,1:2)];
    ves_ext{i}(6,:) = ves(5,1);
    ves_ext{i}(7,:) = [ves(4,4:8),ves(4,1:5), ves(4,6)];
    
    vis_ext{i}(1,:) = [vis(2,4:8),vis(2,1:5), ves(2,6)];
    vis_ext{i}(2,:) = vis(1,1);
    vis_ext{i}(3,:) = [vis(2,8), vis(2,1:8), vis(2,1:2)];
    vis_ext{i}(4,:) = [vis(3,8), vis(3,1:8), vis(3,1:2)];
    vis_ext{i}(5,:) = [vis(4,8), vis(4,1:8), vis(4,1:2)];
    vis_ext{i}(6,:) = vis(5,1);
    vis_ext{i}(7,:) = [vis(4,4:8),vis(4,1:5), vis(4,6)];
    % translate data into 3D
    ves_fit{i} = interp2(xx,yy,ves_ext{i},xxs,yys,'spline');
    vis_fit{i} = interp2(xx,yy,vis_ext{i},xxs,yys,'spline');
    % make responses larger than 1 (should be 0, but for log use later)
    ves_fit{i}(find(ves_fit{i}<1)) = 1;
    vis_fit{i}(find(vis_fit{i}<1)) = 1;
    i    
    % calculate varience/mean ratio for each cell
    aa(find(aa<1)) = 1; % to avoid Log10 to be 0;
    bb(find(bb<1)) = 1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Tunning Width    
%     % cropping the extra values off and subtract minimum
%     ves_fit_ori = ves_fit(46:226, 46:405)-min(min(ves_fit(46:226, 46:405)));    
%     [max_ele_temp, max_azi_temp] = find( ves_fit_ori==max(max(ves_fit_ori)) );
%     max_ele_ves = max_ele_temp(1);
%     max_azi_ves = max_azi_temp(1);
%     ves_fit_ori_azi = ves_fit_ori(max_ele_ves,:) - min(ves_fit_ori(max_ele_ves,:)); % for azimuth width searching
%     ves_fit_ori_ele = ves_fit_ori(:,max_azi_ves) - min(ves_fit_ori(:,max_azi_ves)); % for elevation width searching
%     half_height_ves_azi = (max(max(ves_fit_ori))-min(ves_fit_ori(max_ele_ves,:)))/2;
%     half_height_ves_ele = (max(max(ves_fit_ori))-min(ves_fit_ori(:,max_azi_ves)))/2;
%      
%     % vestibular 
%     % azimuth direction first
%     ves_azi_left = ves_fit_ori_azi(max_azi_ves);
%     ves_azi_right = ves_fit_ori_azi(max_azi_ves);
%     cc_azi_left_ves = 1;
%     while (ves_azi_left - half_height_ves_azi) > 0  & cc_azi_left_ves < 360 % left side of azimuth
%         if max_azi_ves-cc_azi_left_ves >= 1 
%            ves_azi_left = ves_fit_ori_azi(max_azi_ves-cc_azi_left_ves);
%         else
%            ves_azi_left = ves_fit_ori_azi(max_azi_ves-cc_azi_left_ves+360); 
%         end
%         cc_azi_left_ves = cc_azi_left_ves + 1;
%     end
%     cc_azi_right_ves = 1;
%     while (ves_azi_right - half_height_ves_azi) > 0  & cc_azi_right_ves < 360% right side of azimuth
%         if max_azi_ves+ cc_azi_right_ves <=360 
%            ves_azi_right = ves_fit_ori_azi(max_azi_ves+cc_azi_right_ves);
%         else
%            ves_azi_right = ves_fit_ori_azi(max_azi_ves+cc_azi_right_ves - 360);  
%         end
%         cc_azi_right_ves = cc_azi_right_ves + 1;
%     end
%     if cc_azi_left_ves >=360 & cc_azi_right_ves < 360  % only use one side
%        ves_azi_width(i) = 2 * cc_azi_right_ves;
%    elseif cc_azi_left_ves <360 & cc_azi_right_ves >= 360 % only use one side
%        ves_azi_width(i) = 2 * cc_azi_left_ves;
%    else
%        ves_azi_width(i) = cc_azi_left_ves + cc_azi_right_ves;
%    end
%     % elevation, elevation is not circular, it stops at -90 and 90
%     ves_ele_up = ves_fit_ori_ele(max_ele_ves); % up is -90 deg
%     ves_ele_down = ves_fit_ori_ele(max_ele_ves); % down is 90 deg
%     cc_ele_up_ves = 1;
%     while (ves_ele_up - half_height_ves_ele) > 0  & (max_ele_ves - cc_ele_up_ves) > 0 % up side of elevation
%         ves_ele_up = ves_fit_ori_ele(max_ele_ves-cc_ele_up_ves);      
%         cc_ele_up_ves = cc_ele_up_ves + 1;
%     end
%     cc_ele_down_ves = 1;
%     while (ves_ele_down - half_height_ves_ele) > 0  & (max_ele_ves + cc_ele_down_ves) < 181 % down side of elevation
%         ves_ele_down = ves_fit_ori_ele(max_ele_ves+cc_ele_down_ves);
%         cc_ele_down_ves = cc_ele_down_ves + 1;
%     end
%     ves_ele_width(i) = cc_ele_up_ves + cc_ele_down_ves;
%     
%     % visual 
%     vis_fit_ori = vis_fit(46:226, 46:405)-min(min(vis_fit(46:226, 46:405)));
%     [max_ele_temp, max_azi_temp] = find( vis_fit_ori==max(max(vis_fit_ori)) );
%     max_ele_vis = max_ele_temp(1);
%     max_azi_vis = max_azi_temp(1);
%     vis_fit_ori_azi = vis_fit_ori(max_ele_vis,:) - min(vis_fit_ori(max_ele_vis,:)); % for azimuth width searching
%     vis_fit_ori_ele = vis_fit_ori(:,max_azi_vis) - min(vis_fit_ori(:,max_azi_vis)); % for elevation width searching
%     half_height_vis_azi = (max(max(vis_fit_ori))-min(vis_fit_ori(max_ele_vis,:)))/2;
%     half_height_vis_ele = (max(max(vis_fit_ori))-min(vis_fit_ori(:,max_azi_vis)))/2;
%      
%     % visual 
%     % azimuth direction first
%     vis_azi_left = vis_fit_ori_azi(max_azi_vis);
%     vis_azi_right = vis_fit_ori_azi(max_azi_vis);
%     cc_azi_left_vis = 1;
%     while (vis_azi_left - half_height_vis_azi) > 0  & cc_azi_left_vis < 360 % left side of azimuth
%         if max_azi_vis-cc_azi_left_vis >= 1 
%            vis_azi_left = vis_fit_ori_azi(max_azi_vis-cc_azi_left_vis);
%         else
%            vis_azi_left = vis_fit_ori_azi(max_azi_vis-cc_azi_left_vis+360); 
%         end
%         cc_azi_left_vis = cc_azi_left_vis + 1;
%     end
%     cc_azi_right_vis = 1;
%     while (vis_azi_right - half_height_vis_azi) > 0  & cc_azi_right_vis < 360% right side of azimuth
%         if max_azi_vis+ cc_azi_right_vis <=360 
%            vis_azi_right = vis_fit_ori_azi(max_azi_vis+cc_azi_right_vis);
%         else
%            vis_azi_right = vis_fit_ori_azi(max_azi_vis+cc_azi_right_vis - 360);  
%         end
%         cc_azi_right_vis = cc_azi_right_vis + 1;
%     end
%     if cc_azi_left_vis >=360 & cc_azi_right_vis < 360  % only use one side
%        vis_azi_width(i) = 2 * cc_azi_right_vis;
%    elseif cc_azi_left_vis <360 & cc_azi_right_vis >= 360 % only use one side
%        vis_azi_width(i) = 2 * cc_azi_left_vis;
%    else
%        vis_azi_width(i) = cc_azi_left_vis + cc_azi_right_vis;
%    end
%     % elevation, elevation is not circular, it stops at -90 and 90
%     vis_ele_up = vis_fit_ori_ele(max_ele_vis); % up is -90 deg
%     vis_ele_down = vis_fit_ori_ele(max_ele_vis); % down is 90 deg
%     cc_ele_up_vis = 1;
%     while (vis_ele_up - half_height_vis_ele) > 0  & (max_ele_vis - cc_ele_up_vis) > 0 % up side of elevation
%         vis_ele_up = vis_fit_ori_ele(max_ele_vis-cc_ele_up_vis);      
%         cc_ele_up_vis = cc_ele_up_vis + 1;
%     end
%     cc_ele_down_vis = 1;
%     while (vis_ele_down - half_height_vis_ele) > 0  & (max_ele_vis + cc_ele_down_vis) < 181 % down side of elevation
%         vis_ele_down = vis_fit_ori_ele(max_ele_vis+cc_ele_down_vis);
%         cc_ele_down_vis = cc_ele_down_vis + 1;
%     end
%     vis_ele_width(i) = cc_ele_up_vis + cc_ele_down_vis;
    

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Slope and Fisher information

    if max(aa(i,1:26)) == 1  % polyfit might have problem with such data      
       vmr_ves = [1,0];
    elseif max(aa(i,27:52)) == 1
       vmr_vis = [1,0];
    else
        vmr_ves = polyfit( log10(aa(i,1:26)), log10(bb(i,1:26).^2),1 );
        vmr_vis = polyfit( log10(aa(i,27:52)), log10(bb(i,27:52).^2),1 );   
    end
    varience_ves = 10.^ (log10(ves_fit{i})*vmr_ves(1)+vmr_ves(2) );
    varience_vis = 10.^ (log10(vis_fit{i})*vmr_vis(1)+vmr_vis(2) );
    ves_fit{i}(find(ves_fit{i}<1)) = 1;
    vis_fit{i}(find(vis_fit{i}<1)) = 1;
    
    % find steepest slope, for each point, calculate the slope based on +5 and
    % -5 degree
    for j = 0 : slide_size : 359  % azimuth
 %       slope_ves_azi_temp{i}(:, j+1) = abs( ( ves_fit{i}(:,j+45+1+step)-ves_fit{i}(:,j+45+1-step) ) / 2*step );
        slope_vis_azi_temp{i}(:, j+1) = vis_fit{i}(:,j+45+1+1)-vis_fit{i}(:,j+45+1-0);
        % Fisher information
  %      fisher_ves_azi_temp{i}(:, j+1) = (( ves_fit{i}(:,j+45+1+step)-ves_fit{i}(:,j+45+1-step) ) / (2*step)).^2 ./ varience_ves(:,j+45+1);
  %      fisher_vis_azi_temp{i}(:, j+1) = (( vis_fit{i}(:,j+45+1+step)-vis_fit{i}(:,j+45+1-step) ) / (2*step)).^2 ./ varience_vis(:,j+45+1);
        fisher_ves_azi_temp{i}(:, j+1) = (( ves_fit{i}(:,j+45+1+1)-ves_fit{i}(:,j+45+1-0) ) / (2*0.5)).^2 ./ varience_ves(:,j+45+1);
        fisher_vis_azi_temp{i}(:, j+1) = (( vis_fit{i}(:,j+45+1+1)-vis_fit{i}(:,j+45+1-0) ) / (2*0.5)).^2 ./ varience_vis(:,j+45+1);
    end
    for m = -90 : slide_size : 90 % elevation
%         slope_ves_ele_temp{i}(m+91, :) = abs( ( ves_fit{i}(m+135+1+step,:)-ves_fit{i}(m+135+1-step,:) ) / 2*step );  
%         slope_vis_ele_temp{i}(m+91, :) = abs( ( vis_fit{i}(m+135+1+step,:)-vis_fit{i}(m+135+1-step,:) ) / 2*step );  
        % fisher information
%        fisher_ves_ele_temp{i}(m+91, :) = (( ves_fit{i}(m+135+1+step,:)-ves_fit{i}(m+135+1-step,:) ) / (2*step)).^2 ./ varience_ves(m+135+1,:);  
%        fisher_vis_ele_temp{i}(m+91, :) = (( vis_fit{i}(m+135+1+step,:)-vis_fit{i}(m+135+1-step,:) ) / (2*step)).^2 ./ varience_vis(m+135+1,:);  
        fisher_ves_ele_temp{i}(m+91, :) = (( ves_fit{i}(m+135+1+1,:)-ves_fit{i}(m+135+1-0,:) ) / (2*0.5)).^2 ./ varience_ves(m+135+1,:);  
        fisher_vis_ele_temp{i}(m+91, :) = (( vis_fit{i}(m+135+1+1,:)-vis_fit{i}(m+135+1-0,:) ) / (2*0.5)).^2 ./ varience_vis(m+135+1,:);  

    end 
    
    % only choose significant ones 
    if p_ves(i)<=0.05
        fisher_ves_azi = fisher_ves_azi + fisher_ves_azi_temp{i};
        fisher_ves_ele = fisher_ves_ele + fisher_ves_ele_temp{i};
        n_ves = n_ves+1;
    else
        fisher_ves_azi = fisher_ves_azi + 0;
        fisher_ves_ele = fisher_ves_ele + 0;
    end
    if p_vis(i)<=0.05
        fisher_vis_azi = fisher_vis_azi + fisher_vis_azi_temp{i};
        fisher_vis_ele = fisher_vis_ele + fisher_vis_ele_temp{i};
        n_vis = n_vis+1;
    else
        fisher_vis_azi = fisher_vis_azi + 0;
        fisher_vis_ele = fisher_vis_ele + 0;
    end

    [max_azi_ves_ele_temp, max_azi_ves_azi_temp] = find(fisher_ves_azi_temp{i}(46:225,:) == max(max(fisher_ves_azi_temp{i}(46:225,:))) );
    max_azi_ves_ele(i) = max_azi_ves_ele_temp(1)-91;
    max_azi_ves_azi(i) = max_azi_ves_azi_temp(1)-1;
    [max_ele_ves_ele_temp, max_ele_ves_azi_temp] = find(fisher_ves_ele_temp{i}(:,46:404) == max(max(fisher_ves_ele_temp{i}(:,46:404))) );
    max_ele_ves_ele(i) = max_ele_ves_ele_temp(1)-91;
    max_ele_ves_azi(i) = max_ele_ves_azi_temp(1)-1;

    [max_azi_vis_ele_temp, max_azi_vis_azi_temp] = find(fisher_vis_azi_temp{i}(46:225,:) == max(max(fisher_vis_azi_temp{i}(46:225,:))) );
    max_azi_vis_ele(i) = max_azi_vis_ele_temp(1)-91;
    max_azi_vis_azi(i) = max_azi_vis_azi_temp(1)-1;
    [max_ele_vis_ele_temp, max_ele_vis_azi_temp] = find(fisher_vis_ele_temp{i}(:,46:404) == max(max(fisher_vis_ele_temp{i}(:,46:404))) );
    max_ele_vis_ele(i) = max_ele_vis_ele_temp(1)-91;
    max_ele_vis_azi(i) = max_ele_vis_azi_temp(1)-1;
end
% n_ves
% n_vis
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot distribution of maximum slope in azi and ele
% figure(1);
% set(1,'Position', [5,5 1000,680], 'Name', 'Envelope');
% orient landscape;
% subplot(2,2,1);
% plot(max_azi_ves_azi(1:200),max_azi_ves_ele(1:200),'o');
% xlim([0, 360]);
% ylim([-90, 90]);
% set(gca, 'ydir' , 'reverse');
% set(gca, 'ytick',[-90,-45,0,45,90]);
% set(gca, 'yticklabel','-90|-45|0|45|90');
% set(gca, 'xtick',[0,90,180,270,360]);
% set(gca, 'xticklabel','0|90|180|270|360'); 
% title('ves: max azi');
% 
% subplot(2,2,2);
% plot(max_ele_ves_azi(1:200),max_ele_ves_ele(1:200),'o');
% xlim([0, 360]);
% ylim([-90, 90]);
% set(gca, 'ydir' , 'reverse');
% set(gca, 'ytick',[-90,-45,0,45,90]);
% set(gca, 'yticklabel','-90|-45|0|45|90');
% set(gca, 'xtick',[0,90,180,270,360]);
% set(gca, 'xticklabel','0|90|180|270|360');
% title('ves: max ele');
% 
% subplot(2,2,3);
% plot(max_azi_vis_azi,max_azi_vis_ele,'o');
% xlim([0, 360]);
% ylim([-90, 90]);
% set(gca, 'ydir' , 'reverse');
% set(gca, 'ytick',[-90,-45,0,45,90]);
% set(gca, 'yticklabel','-90|-45|0|45|90');
% set(gca, 'xtick',[0,90,180,270,360]);
% set(gca, 'xticklabel','0|90|180|270|360');
% title('vis: max azi');
% 
% subplot(2,2,4);
% plot(max_ele_vis_azi,max_ele_vis_ele,'o');
% xlim([0, 360]);
% ylim([-90, 90]);
% set(gca, 'ydir' , 'reverse');
% set(gca, 'ytick',[-90,-45,0,45,90]);
% set(gca, 'yticklabel','-90|-45|0|45|90');
% set(gca, 'xtick',[0,90,180,270,360]);
% set(gca, 'xticklabel','0|90|180|270|360');
% title('vis: max ele');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
set(2,'Position', [5,5 1000,680], 'Name', 'Envelope');
orient landscape;
subplot(2,2,1);
contourf(sqrt(fisher_ves_azi(:,:)));
colorbar;
ylim([46,226]);
set(gca, 'ydir' , 'reverse');
set(gca, 'ytick',[46,91,136,181,226]);
set(gca, 'yticklabel','-90|-45|0|45|90');
set(gca, 'xtick',[1,90,180,270,360]);
set(gca, 'xticklabel','0|90|180|270|360'); 
title('Fisher in azimuth: Ves');

subplot(2,2,2);
contourf(sqrt(fisher_ves_ele(:,:)));
colorbar;
xlim([46,406]);
set(gca, 'ydir' , 'reverse');
set(gca, 'ytick',[1,46,91,136,181]);
set(gca, 'yticklabel','-90|-45|0|45|90');
set(gca, 'xtick',[46,136,226,316,406]);
set(gca, 'xticklabel','0|90|180|270|360'); 
title('Fisher in elevation: Ves');

subplot(2,2,3);
contourf(sqrt(fisher_vis_azi(:,:)));
colorbar;
ylim([46,226]);
set(gca, 'ydir' , 'reverse');
set(gca, 'ytick',[46,91,136,181,226]);
set(gca, 'yticklabel','-90|-45|0|45|90');
set(gca, 'xtick',[1,90,180,270,360]);
set(gca, 'xticklabel','0|90|180|270|360'); 
title('Fisher in azimuth: Vis');

subplot(2,2,4);
contourf(sqrt(fisher_vis_ele(:,:)));
colorbar;
xlim([46,406]);
set(gca, 'ydir' , 'reverse');
set(gca, 'ytick',[1,46,91,136,181]);
set(gca, 'yticklabel','-90|-45|0|45|90');
set(gca, 'xtick',[46,136,226,316,406]);
set(gca, 'xticklabel','0|90|180|270|360'); 
title('Fisher in elevation: Vis');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % Plot each cell with fitted curve
% 
% % vestibular
% % plot raw first
% for f = 1:5
%     figure(f+1);
%     set(f+1,'Position', [5,5 1000,680], 'Name', 'Envelope');
%     orient landscape;
% 	axis off;
% 	for j = 1 : 5     
%         for i = 1 : 4  
%             if (j-1)*4+(f-1)*20+i <= dim(1)
%                 axes('position',[0.18*(j-1)+0.02 (0.98-0.22*i) 0.16 0.18]);  
%                 % raw data
%                 contourf(ves_ext{(j-1)*4+(f-1)*20+i});  
%                 xlim([2,10]);
%                 ylim([2,6]);
%                 set(gca, 'ydir' , 'reverse'); % upward motion is negative, on top of the y axis
%                 if j==1 & i~=4
%                    set(gca, 'yticklabel','-90|-45|0|45|90');
%                    set(gca, 'xtick',[]);
%                 elseif i==4 & j~=1
%                    set(gca, 'xticklabel','0|90|180|270|360');  
%                    set(gca, 'ytick',[]);
%                 elseif i==4 & j==1
%                    set(gca, 'xticklabel','0|90|180|270|360');
%                    set(gca, 'yticklabel','-90|-45|0|45|90');
%                 else
%                    set(gca, 'xtick',[]);
%                    set(gca, 'ytick',[]);
%                 end
%                 title(num2str(p_ves((j-1)*4+(f-1)*20+i))); 
%             end
%         end 
% 	end
% end
% % now spline fitted data
% for f = 1:5
% 	figure(f+1+5);
% 	set(f+1+5,'Position', [5,5 1000,680], 'Name', 'Envelope');
% 	orient landscape;
% 	axis off;
% 	for j = 1 : 5     
%         for i = 1 : 4      
%             if (j-1)*4+(f-1)*20+i <= dim(1)
%                 axes('position',[0.18*(j-1)+0.02 (0.98-0.22*i) 0.16 0.18]);  
%                 % fit data
%                 contourf(ves_fit{(j-1)*4+(f-1)*20+i});  
%                 xlim([46,406]);
%                 ylim([46,226]);
%                 set(gca, 'ydir' , 'reverse'); % upward motion is negative, on top of the y axis
%                 if j==1 & i~=4
%                    set(gca, 'yticklabel','-90|-45|0|45|90');
%                    set(gca, 'xtick',[]);
%                 elseif i==4 & j~=1
%                    set(gca, 'xtick',[46,136,226,316,406]);
%                    set(gca, 'xticklabel','0|90|180|270|360'); 
%                    set(gca, 'ytick',[]);
%                 elseif i==4 & j==1
%                    set(gca, 'xtick',[46,136,226,316,406]);
%                    set(gca, 'xticklabel','0|90|180|270|360');
%                    set(gca, 'yticklabel','-90|-45|0|45|90');
%                 else
%                    set(gca, 'xtick',[]);
%                    set(gca, 'ytick',[]);
%                 end
%                 title(num2str(p_ves((j-1)*4+(f-1)*20+i))); 
%                 hold on;
%                 % superimpose the maximum slope 
%                 xa = max_azi_ves_azi((j-1)*4+(f-1)*20+i)+46;
%                 ya = max_azi_ves_ele((j-1)*4+(f-1)*20+i)+136;
%                 xe = max_ele_ves_azi((j-1)*4+(f-1)*20+i)+46;
%                 ye = max_ele_ves_ele((j-1)*4+(f-1)*20+i)+136;        
%                 plot(xa, ya, 'o','MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 0] );  % max azi slope
%                 plot(xe, ye, 'o','MarkerEdgeColor','w', 'MarkerFaceColor',[1 1 1]);  % max ele slope
%             end
%         end 
% 	end
% end
% 
% % visual
% % plot raw first
% for f = 1:5
%     figure(f+1+10);
%     set(f+1+10,'Position', [5,5 1000,680], 'Name', 'Envelope');
%     orient landscape;
% 	axis off;
% 	for j = 1 : 5     
%         for i = 1 : 4  
%             if (j-1)*4+(f-1)*20+i <= dim(1)
%                 axes('position',[0.18*(j-1)+0.02 (0.98-0.22*i) 0.16 0.18]);  
%                 % raw data
%                 contourf(vis_ext{(j-1)*4+(f-1)*20+i});  
%                 xlim([2,10]);
%                 ylim([2,6]);
%                 set(gca, 'ydir' , 'reverse'); % upward motion is negative, on top of the y axis
%                 if j==1 & i~=4
%                    set(gca, 'yticklabel','-90|-45|0|45|90');
%                    set(gca, 'xtick',[]);
%                 elseif i==4 & j~=1
%                    set(gca, 'xticklabel','0|90|180|270|360');  
%                    set(gca, 'ytick',[]);
%                 elseif i==4 & j==1
%                    set(gca, 'xticklabel','0|90|180|270|360');
%                    set(gca, 'yticklabel','-90|-45|0|45|90');
%                 else
%                    set(gca, 'xtick',[]);
%                    set(gca, 'ytick',[]);
%                 end
%                 title(num2str(p_vis((j-1)*4+(f-1)*20+i))); 
%             end
%         end 
% 	end
% end
% % now spline fitted data
% for f = 1:5
% 	figure(f+1+15);
% 	set(f+1+15,'Position', [5,5 1000,680], 'Name', 'Envelope');
% 	orient landscape;
% 	axis off;
% 	for j = 1 : 5     
%         for i = 4 : 4      
%             if (j-1)*4+(f-1)*20+i <= dim(1)
%                 axes('position',[0.18*(j-1)+0.02 (0.98-0.22*i) 0.16 0.18]);  
%                 % fit data
%                 contourf(vis_fit{(j-1)*4+(f-1)*20+i});  
%                 xlim([46,406]);
%                 ylim([46,226]);
%                 set(gca, 'ydir' , 'reverse'); % upward motion is negative, on top of the y axis
%                 if j==1 & i~=4
%                    set(gca, 'yticklabel','-90|-45|0|45|90');
%                    set(gca, 'xtick',[]);
%                 elseif i==4 & j~=1
%                    set(gca, 'xtick',[46,136,226,316,406]);
%                    set(gca, 'xticklabel','0|90|180|270|360'); 
%                    set(gca, 'ytick',[]);
%                 elseif i==4 & j==1
%                    set(gca, 'xtick',[46,136,226,316,406]);
%                    set(gca, 'xticklabel','0|90|180|270|360');
%                    set(gca, 'yticklabel','-90|-45|0|45|90');
%                 else
%                    set(gca, 'xtick',[]);
%                    set(gca, 'ytick',[]);
%                 end
%                 title(num2str(p_vis((j-1)*4+(f-1)*20+i))); 
%                 hold on;
%                 % superimpose the maximum slope 
%                 xa = max_azi_vis_azi((j-1)*4+(f-1)*20+i)+46;
%                 ya = max_azi_vis_ele((j-1)*4+(f-1)*20+i)+136;
%                 xe = max_ele_vis_azi((j-1)*4+(f-1)*20+i)+46;
%                 ye = max_ele_vis_ele((j-1)*4+(f-1)*20+i)+136;        
%                 plot(xa, ya, 'o','MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 0] );  % max azi slope
%                 plot(xe, ye, 'o','MarkerEdgeColor','w', 'MarkerFaceColor',[1 1 1]);  % max ele slope
%             end
%         end 
% 	end
% end

